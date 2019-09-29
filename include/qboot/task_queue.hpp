#ifndef QBOOT_TASK_QUEUE_HPP_
#define QBOOT_TASK_QUEUE_HPP_

#include <atomic>              // for atomic
#include <condition_variable>  // for condition_variable
#include <functional>          // for function
#include <future>              // for future, promise
#include <memory>              // for unique_ptr, make_unique
#include <mutex>               // for mutex, lock_guard, unique_lock
#include <queue>               // for queue
#include <string>              // for string
#include <string_view>         // for string_view
#include <thread>              // for thread
#include <type_traits>         // for is_default_constructible_v
#include <utility>             // for declval, move
#include <vector>              // for vector

namespace qboot
{
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
	template <class T>
	std::vector<T> parallel_evaluate(const std::vector<std::function<T()>>& fs,
	                                 uint32_t p = std::thread::hardware_concurrency())
	{
		static_assert(std::is_default_constructible_v<T>);
		if (p <= 1) return _seq_eval(fs);
		auto N = fs.size();
		std::vector<T> ans(N);
		std::mutex mtx{};
		uint32_t now = 0;
		std::vector<std::thread> worker;
		for (uint32_t i = 0; i < p; ++i)
			worker.emplace_back([&mtx, &now, N, &fs, &ans]() {
				while (true)
				{
					uint32_t now_local;
					{
						std::lock_guard<std::mutex> lock(mtx);
						now_local = now;
						++now;
					}
					if (now_local >= N) break;
					ans[now_local] = fs[now_local]();
				}
			});
		for (uint32_t i = 0; i < p; ++i) worker[i].join();
		return ans;
	}
	// evaluate fs in parallel
	// if p <= 1, evaluate sequentially
	inline void parallel_evaluate(const std::vector<std::function<void()>>& fs,
	                              uint32_t p = std::thread::hardware_concurrency())
	{
		if (p <= 1) return _seq_eval(fs);
		std::vector<std::function<bool()>> fs_dummy;
		for (const auto& f : fs)
			fs_dummy.emplace_back([&f]() {
				f();
				return true;
			});
		parallel_evaluate(fs_dummy, p);
	}

	class _task_base
	{
	public:
		_task_base() = default;
		_task_base(const _task_base&) = default;
		_task_base(_task_base&&) noexcept = default;
		_task_base& operator=(const _task_base&) = default;
		_task_base& operator=(_task_base&&) noexcept = default;
		virtual ~_task_base();
		virtual void run() = 0;
	};

	template <class T>
	class _task : public _task_base
	{
		std::promise<T> p_;
		std::function<T()> f_;

	public:
		_task(std::promise<T>&& p, std::function<T()> f) : _task_base(), p_(std::move(p)), f_(std::move(f)) {}
		_task(const _task&) = default;
		_task(_task&&) noexcept = default;
		_task& operator=(const _task&) = default;
		_task& operator=(_task&&) noexcept = default;
		~_task() override = default;
		void run() override { p_.set_value(f_()); }
	};

	class TaskQueue
	{
		std::mutex mtx_{};                // lock for q_ and killed_
		std::condition_variable cond_{};  // awake if task is queued or killed
		std::vector<std::thread> ts_{};   // worker threads
		std::atomic<bool> killed_ = false;
		std::queue<std::unique_ptr<_task_base>> q_{};
		void work()
		{
			std::unique_ptr<_task_base> task;
			while (true)
			{
				{
					std::unique_lock<std::mutex> lk(mtx_);
					cond_.wait(lk, [this] { return killed_ || !q_.empty(); });
					if (killed_) return;
					task = std::move(q_.front());
					q_.pop();
				}
				task->run();
				task.reset();
			}
		}

	public:
		TaskQueue(uint32_t p = std::thread::hardware_concurrency())
		{
			for (uint32_t i = 0; i < p; ++i) ts_.emplace_back([this] { work(); });
		}
		~TaskQueue()
		{
			signal_done();
			for (auto& t : ts_) t.join();
		}
		void signal_done() &
		{
			std::unique_lock<std::mutex> lk(mtx_);
			if (killed_.exchange(true)) return;
			cond_.notify_all();
		}
		template <class Func, class T = decltype(std::declval<Func>()())>
		std::future<T> push(Func task) &
		{
			std::unique_lock<std::mutex> lk(mtx_);
			std::promise<T> p_;
			auto f = p_.get_future();
			q_.push(std::make_unique<_task<T>>(std::move(p_), task));
			cond_.notify_one();
			return f;
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
		explicit _scoped_event(std::string_view tag, _event_base* event = nullptr) : tag_(tag), event_(event)
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
		template <class T>
		_scoped_event(std::string_view tag [[maybe_unused]], [[maybe_unused]] T&& event)
		{
		}
	};
#endif
}  // namespace qboot

#endif  // QBOOT_TASK_QUEUE_HPP_
