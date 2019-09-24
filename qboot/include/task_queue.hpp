#ifndef QBOOT_TASK_QUEUE_HPP_
#define QBOOT_TASK_QUEUE_HPP_

#include <condition_variable>  // for condition_variable
#include <functional>          // for function
#include <future>              // for future, promise
#include <memory>              // for unique_ptr
#include <mutex>               // for mutex, unique_lock
#include <queue>               // for queue
#include <string>              // for string
#include <string_view>         // for string_view
#include <thread>              // for thread
#include <utility>             // for declval, move
#include <vector>              // for vector

namespace qboot
{
	// producer-consumer pattern
	// T must be default constructive
	template <class T>
	class _bounded_queue
	{
		std::queue<T> q_{};
		std::mutex mtx_{};
		std::condition_variable not_empty_{};

	public:
		// ensure q.size() >= len
		void fill_dummy(uint32_t len) &
		{
			std::unique_lock<std::mutex> lk(mtx_);
			while (q_.size() < len) q_.push({});
			if (len > 0) not_empty_.notify_all();
		}
		void push(T val) &
		{
			std::unique_lock<std::mutex> lk(mtx_);
			q_.push(std::move(val));
			not_empty_.notify_all();
		}
		T pop()
		{
			std::unique_lock<std::mutex> lk(mtx_);
			not_empty_.wait(lk, [this] { return !q_.empty(); });
			T ret = std::move(q_.front());
			q_.pop();
			return ret;
		}
	};
	struct _task_base
	{
		_task_base() = default;
		_task_base(const _task_base&) = default;
		_task_base(_task_base&&) noexcept = default;
		_task_base& operator=(const _task_base&) = default;
		_task_base& operator=(_task_base&&) noexcept = default;
		virtual ~_task_base();
		virtual void run() = 0;
	};

	template <class T>
	struct _task : public _task_base
	{
		_task(std::promise<T>&& p, std::function<T()> f) : _task_base(), p_(std::move(p)), f_(std::move(f)) {}
		_task(const _task&) = default;
		_task(_task&&) noexcept = default;
		_task& operator=(const _task&) = default;
		_task& operator=(_task&&) noexcept = default;
		~_task() override = default;
		void run() override { p_.set_value(f_()); }

	private:
		std::promise<T> p_;
		std::function<T()> f_;
	};
	class TaskQueue
	{
		bool fin_ = false;
		_bounded_queue<std::unique_ptr<_task_base>> q_{};
		std::vector<std::thread> ts_{};
		void consume()
		{
			while (!fin_)
			{
				auto t = q_.pop();
				if (t) t->run();
			}
		}

	public:
		// cap: capacity for internal job queue
		// p: actual concurrency
		explicit TaskQueue(uint32_t p = std::thread::hardware_concurrency())
		{
			for (uint32_t i = 0; i < p; i++) ts_.emplace_back([this]() { consume(); });
		}
		TaskQueue(const TaskQueue&) = delete;
		TaskQueue(TaskQueue&&) = delete;
		TaskQueue& operator=(const TaskQueue&) = delete;
		TaskQueue& operator=(TaskQueue&&) = delete;
		~TaskQueue()
		{
			fin_ = true;
			q_.fill_dummy(uint32_t(ts_.size()));
			for (auto&& t : ts_)
				if (t.joinable()) t.detach();
		}
		template <class Func, class T = decltype(std::declval<Func>()())>
		std::future<T> push(Func f)
		{
			std::promise<T> p;
			auto fut = p.get_future();
			q_.push(std::make_unique<_task<T>>(std::move(p), f));
			return fut;
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
#define QBOOT_scope(varname, tag, event) _scoped_event varname(tag, event)
#else
#define QBOOT_scope(varname, tag, event) ((void)0)
#endif
}  // namespace qboot

#endif  // QBOOT_TASK_QUEUE_HPP_
