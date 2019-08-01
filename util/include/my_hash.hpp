#ifndef MY_HASH_HPP_
#define MY_HASH_HPP_

#include <array>        // for array
#include <cstddef>      // for size_t
#include <optional>     // for optional
#include <tuple>        // for tuple, apply
#include <type_traits>  // for declval
#include <utility>      // for pair
#include <variant>      // for variant

namespace util
{
	using std::size_t;
	constexpr inline void hash_combine(size_t* h, size_t k) noexcept
	{
		constexpr size_t m = 0xc6a4a7935bd1e995ULL;
		constexpr int r = 47;

		k *= m;
		k ^= k >> r;
		k *= m;

		*h ^= k;
		*h *= m;

		*h += 0xe6546b64;
	}
	union _double_and_size_t {
		double d;
		size_t s;
		constexpr _double_and_size_t(const double& x) noexcept : d(x) {}
	};
	struct MyHash
	{
		constexpr size_t operator()(const int64_t& r) const noexcept { return size_t(r); }
		constexpr size_t operator()(const size_t& r) const noexcept { return r; }
		constexpr size_t operator()(const double& r) const noexcept
		{
			// -0 and +0 must be equal
			if (-0.0 <= r && r <= 0.0) return 0;
			_double_and_size_t u(r);
			return (*this)(u.s);
		}
		template <class T,
		          class = std::void_t<std::enable_if<std::is_same_v<std::size_t, decltype(std::declval<T>().hash())>>>>
		constexpr size_t operator()(const T& r) const noexcept
		{
			return r.hash();
		}
		template <size_t n, typename T>
		constexpr size_t operator()(const std::array<T, n>& a) const noexcept
		{
			size_t hash = 0;
			for (auto& x : a) hash_combine(&hash, (*this)(x));
			return hash;
		}
		template <class... Types>
		constexpr size_t operator()(const std::tuple<Types...>& t) const noexcept
		{
			size_t hash = 0;
			std::apply([&](auto... args) constexpr { (hash_combine(&hash, (*this)(args)), ...); }, t);
			return hash;
		}
		template <class T1, class T2>
		constexpr size_t operator()(const std::pair<T1, T2>& p) const noexcept
		{
			size_t hash = (*this)(p.first);
			hash_combine(&hash, (*this)(p.second));
			return hash;
		}
		template <class... Types>
		constexpr size_t operator()(const std::variant<Types...>& t) const noexcept
		{
			size_t hash = std::visit([this](auto&& arg) { return (*this)(arg); }, t);
			hash_combine(&hash, (*this)(t.index()));
			return hash;
		}
		template <class T>
		constexpr size_t operator()(const std::optional<T>& p) const noexcept
		{
			if (!p.has_value()) return 0;
			size_t hash = 1;
			hash_combine(&hash, (*this)(p));
			return hash;
		}
	};
}  // namespace util

#endif  // MY_HASH_HPP_
