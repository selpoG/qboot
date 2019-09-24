#include "task_queue.hpp"

namespace qboot
{
	_task_base::~_task_base() = default;
#ifndef NDEBUG
	_event_base::~_event_base() = default;
#endif
}  // namespace qboot
