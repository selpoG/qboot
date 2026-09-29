#ifndef QBOOT_TASK_QUEUE_HPP_
#define QBOOT_TASK_QUEUE_HPP_

#include <concepts>            // for constructible_from, default_initializable, invocable, same_as
#include <condition_variable>  // for condition_variable_any
#include <cstddef>             // for size_t
#include <cstdint>             // for uint32_t
#include <exception>           // for current_exception, exception_ptr, rethrow_exception
#include <functional>          // for function
#include <future>              // for future, packaged_task
#include <memory>              // for unique_ptr
#include <mutex>               // for mutex, lock_guard, unique_lock
#include <queue>               // for queue
#include <stdexcept>           // for logic_error
#include <stop_token>          // for stop_token
#include <string>              // for string
#include <string_view>         // for string_view
#include <thread>              // for jthread, thread
#include <type_traits>         // for decay_t, invoke_result_t
#include <utility>             // for forward, move
#include <vector>              // for vector

namespace qboot
{
	// MPFR retains per-thread caches unless they are freed before the worker exits.
	void _free_mpfr_cache();

	template <class T>
	std::vector<T> _seq_eval(const std::vector<std::function<T()>>& fs)
	{
		std::vector<T> ans;
		for (const auto& f : fs) ans.push_back(f());
		return ans;
	}
	inline void _seq_eval(const std::vector<std::function<void()>>& fs)
	{
		for (const auto& f : fs) f();
	}
	// evaluate fs in parallel and returns its value
	// if p <= 1, evaluate sequentially
	template <std::default_initializable T>
	std::vector<T> _parallel_evaluate(const std::vector<std::function<T()>>& fs,
	                                  uint32_t p = std::thread::hardware_concurrency())
	{
		if (p <= 1) return _seq_eval(fs);
		auto N = fs.size();
		std::vector<T> ans(N);
		std::mutex mtx{};
		std::size_t now = 0;
		std::exception_ptr error{};
		{
			std::vector<std::jthread> worker;
			for (uint32_t i = 0; i < p; ++i)
				worker.emplace_back([&mtx, &now, N, &fs, &ans, &error] {
					try
					{
						while (true)
						{
							std::size_t now_local;
							{
								std::lock_guard<std::mutex> lock(mtx);
								if (error || now >= N) break;
								now_local = now++;
							}
							if constexpr (std::same_as<T, bool>)
							{
								// vector<bool> packs distinct elements into shared storage.
								const bool value = fs[now_local]();
								std::lock_guard<std::mutex> lock(mtx);
								ans[now_local] = value;
							}
							else
								ans[now_local] = fs[now_local]();
						}
					}
					catch (...)
					{
						std::lock_guard<std::mutex> lock(mtx);
						if (!error) error = std::current_exception();
					}
					_free_mpfr_cache();
				});
		}
		if (error) std::rethrow_exception(error);
		return ans;
	}
	// evaluate fs in parallel
	// if p <= 1, evaluate sequentially
	inline void _parallel_evaluate(const std::vector<std::function<void()>>& fs,
	                               uint32_t p = std::thread::hardware_concurrency())
	{
		if (p <= 1) return _seq_eval(fs);
		std::vector<std::function<bool()>> fs_dummy;
		for (const auto& f : fs)
			fs_dummy.emplace_back([&f] {
				f();
				return true;
			});
		_parallel_evaluate(fs_dummy, p);
	}

	class _task_queue
	{
		std::mutex mtx_{};  // lock for q_ and killed_
		std::condition_variable_any cond_{};
		bool killed_ = false;
		std::queue<std::packaged_task<void()>> q_{};
		// Destroy workers before the state they access, including during construction failure.
		std::vector<std::jthread> ts_{};
		void work(std::stop_token stop)
		{
			while (true)
			{
				std::packaged_task<void()> task;
				{
					std::unique_lock<std::mutex> lk(mtx_);
					if (!cond_.wait(lk, stop, [this] { return killed_ || !q_.empty(); }) || killed_) return;
					task = std::move(q_.front());
					q_.pop();
				}
				task();
			}
		}

	public:
		_task_queue(uint32_t p = std::thread::hardware_concurrency())
		{
			if (p == 0) p = 1;
			for (uint32_t i = 0; i < p; ++i)
				ts_.emplace_back([this](std::stop_token stop) {
					work(stop);
					_free_mpfr_cache();
				});
		}
		~_task_queue() { signal_done(); }
		void signal_done() &
		{
			std::lock_guard<std::mutex> lk(mtx_);
			killed_ = true;
			cond_.notify_all();
		}
		template <class Func>
			requires (std::invocable<std::decay_t<Func>&> && std::constructible_from<std::decay_t<Func>, Func>)
		std::future<std::invoke_result_t<std::decay_t<Func>&>> push(Func&& task) &
		{
			std::lock_guard<std::mutex> lk(mtx_);
			if (killed_) throw std::logic_error("Task queue has stopped");
			std::packaged_task<std::invoke_result_t<std::decay_t<Func>&>()> packaged(std::forward<Func>(task));
			auto result = packaged.get_future();
			q_.emplace([queued = std::move(packaged)]() mutable { queued(); });
			cond_.notify_one();
			return result;
		}
	};

	class _event_base
	{
#ifndef NDEBUG
	public:
		_event_base() = default;
		_event_base(const _event_base&) = default;
		_event_base(_event_base&&) noexcept = default;
		_event_base& operator=(const _event_base&) = default;
		_event_base& operator=(_event_base&&) noexcept = default;
		virtual ~_event_base();
		virtual void on_begin(std::string_view tag) = 0;
		virtual void on_end(std::string_view tag) = 0;
#endif
	};

#ifndef NDEBUG
	class _scoped_event
	{
		std::string tag_;
		_event_base* event_;

	public:
		explicit _scoped_event(std::string_view tag, const std::unique_ptr<_event_base>& event = {})
		    : tag_(tag), event_(event.get())
		{
			if (event_ != nullptr) event_->on_begin(tag_);
		}
		_scoped_event(const _scoped_event&) = delete;
		_scoped_event(_scoped_event&&) = delete;
		_scoped_event& operator=(const _scoped_event&) = delete;
		_scoped_event& operator=(_scoped_event&&) = delete;
		~_scoped_event()
		{
			if (event_ != nullptr) event_->on_end(tag_);
		}
	};
#else
	class _scoped_event
	{
	public:
		explicit _scoped_event([[maybe_unused]] std::string_view tag,
		                       [[maybe_unused]] const std::unique_ptr<_event_base>& event = {})
		{
		}
	};
#endif
}  // namespace qboot

#endif  // QBOOT_TASK_QUEUE_HPP_
