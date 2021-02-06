#ifndef QBOOT_MY_CONCEPTS_HPP_
#define QBOOT_MY_CONCEPTS_HPP_

#include <type_traits>  // for is_constructible_v, is_integral_v, is_nothrow_destructible_v, is_signed_v

namespace qboot
{
	// define concepts not defined by clang
	template <class T>
	concept integral = std::is_integral_v<T>;
	template <class T>
	concept signed_integral = integral<T>&& std::is_signed_v<T>;
	template <class T>
	concept unsigned_integral = integral<T> && !signed_integral<T>;

	template <class T>
	concept _destructible = std::is_nothrow_destructible_v<T>;
	template <class T, class... Args>
	concept _constructible_from = _destructible<T>&& std::is_constructible_v<T, Args...>;
	template <class T>
	concept default_initializable = _constructible_from<T>&& requires
	{
		T{};
	}
	&&requires { ::new (static_cast<void*>(nullptr)) T; };
}  // namespace qboot

#endif  // QBOOT_MY_CONCEPTS_HPP_
