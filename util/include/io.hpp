#ifndef IO_HPP_
#define IO_HPP_

#include <cstddef>      // for size_t
#include <istream>      // for basic_istream
#include <optional>     // for optional, nullopt
#include <ostream>      // for basic_ostream
#include <string>       // for string
#include <type_traits>  // for declval, remove_reference
#include <utility>      // for declval

#include "multi_array.hpp"

namespace util
{
	template <class T>
	class _left_shift
	{
		using type = typename std::remove_reference<T>::type;

	public:
		template <class Char, class Traits>
		std::basic_ostream<Char, Traits>& operator()(std::basic_ostream<Char, Traits>& os, const type& v) const
		{
			return os << v;
		}
	};
	template <
	    class T, class E = decltype(*std::begin(std::declval<T>())), class Printer = _left_shift<E>,
	    class = std::void_t<decltype(std::declval<const Printer&>()(std::declval<std::ostream&>(), std::declval<E>()))>>
	class OutSeqWrap
	{
		const T& val;
		const std::string left, right, sep;
		const Printer& printer;
		const std::optional<std::string> on_empty;

	public:
		explicit OutSeqWrap(const T& val, const std::string& left = "{", const std::string& right = "}",
		                    const std::string& sep = ", ", const std::optional<std::string>& on_empty = std::nullopt,
		                    const Printer& printer = Printer())
		    : val(val), left(left), right(right), sep(sep), printer(printer), on_empty(on_empty)
		{
		}
		OutSeqWrap(const T& val, const Printer& printer, const std::string& left = "{", const std::string& right = "}",
		           const std::string& sep = ", ", const std::optional<std::string>& on_empty = std::nullopt)
		    : val(val), left(left), right(right), sep(sep), printer(printer), on_empty(on_empty)
		{
		}
		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os,
		                                                    const OutSeqWrap<T, E, Printer, void>& v)
		{
			auto b = std::begin(v.val);
			auto e = std::end(v.val);
			if (b == e && v.on_empty.has_value()) return os << v.on_empty.value();
			bool flag = false;
			os << v.left;
			for (; b != e; ++b)
			{
				const auto& x = *b;
				if (flag) os << v.sep;
				v.printer(os, x);
				flag = true;
			}
			return os << v.right;
		}
	};
	template <class T, size_t r, size_t fixed_, class Printer>
	class _out_array_wrap
	{
		multi_array_const_view<T, r, fixed_> val;
		const std::string left, right, sep;
		const Printer& printer;
		const std::optional<std::string> on_empty;

	public:
		_out_array_wrap(const multi_array_const_view<T, r, fixed_>& val, const std::string& left,
		                const std::string& right, const std::string& sep, const std::optional<std::string>& on_empty,
		                const Printer& printer)
		    : val(val), left(left), right(right), sep(sep), printer(printer), on_empty(on_empty)
		{
		}
		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os,
		                                                    const _out_array_wrap& v)
		{
			if constexpr (r == fixed_)
				return v.printer(os, *v.val);
			else
			{
				if (v.val.size() == 0 && v.on_empty.has_value()) return os << v.on_empty.value();
				auto flag = false;
				os << v.left;
				for (size_t x = 0; x < v.val.size(); ++x)
				{
					if (flag) os << v.sep;
					os << _out_array_wrap<T, r, fixed_ + 1, Printer>(v.val[x], v.left, v.right, v.sep, v.on_empty,
					                                                 v.printer);
					flag = true;
				}
				return os << v.right;
			}
		}
	};
	template <class T, size_t r, size_t fixed_, class Printer = _left_shift<T>,
	          class = std::void_t<decltype(std::declval<Printer>()(std::declval<std::ostream&>(), std::declval<T>()))>>
	auto OutArrayWrap(const multi_array_const_view<T, r, fixed_>& val, const std::string& left = "{",
	                  const std::string& right = "}", const std::string& sep = ", ",
	                  const std::optional<std::string>& on_empty = std::nullopt, const Printer& printer = Printer())
	{
		return _out_array_wrap<T, r, fixed_, Printer>(val, left, right, sep, on_empty, printer);
	}
	template <class T, size_t r, size_t fixed_, class Printer = _left_shift<T>,
	          class = std::void_t<decltype(std::declval<Printer>()(std::declval<std::ostream&>(), std::declval<T>()))>>
	auto OutArrayWrap(const multi_array_const_view<T, r, fixed_>& val, const Printer& printer,
	                  const std::string& left = "{", const std::string& right = "}", const std::string& sep = ", ",
	                  const std::optional<std::string>& on_empty = std::nullopt)
	{
		return _out_array_wrap<T, r, fixed_, Printer>(val, left, right, sep, on_empty, printer);
	}
	template <class T, size_t r, class Printer = _left_shift<T>,
	          class = std::void_t<decltype(std::declval<Printer>()(std::declval<std::ostream&>(), std::declval<T>()))>>
	auto OutArrayWrap(const multi_array<T, r>& val, const std::string& left = "{", const std::string& right = "}",
	                  const std::string& sep = ", ", const std::optional<std::string>& on_empty = std::nullopt,
	                  const Printer& printer = Printer())
	{
		return _out_array_wrap<T, r, 0, Printer>(val.get_view(), left, right, sep, on_empty, printer);
	}
	template <class T, size_t r, class Printer>
	auto OutArrayWrap(const multi_array<T, r>& val, const Printer& printer, const std::string& left = "{",
	                  const std::string& right = "}", const std::string& sep = ", ",
	                  const std::optional<std::string>& on_empty = std::nullopt)
	{
		return _out_array_wrap<T, r, 0, Printer>(val.get_view(), left, right, sep, on_empty, printer);
	}
	template <class T, size_t r, class Char, class Traits>
	std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os, const multi_array<T, r>& v)
	{
		return os << OutArrayWrap(v);
	}
}  // namespace util

#endif  // IO_HPP_
