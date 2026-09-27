#include "qboot/task_queue.hpp"

#include "mpfr.h"  // for mpfr_free_cache

namespace qboot
{
	void _free_mpfr_cache() { mpfr_free_cache(); }

	_task_base::~_task_base() = default;
#ifndef NDEBUG
	_event_base::~_event_base() = default;
#endif
}  // namespace qboot
