#ifndef QBOOT_RATIONAL_HPP_
#define QBOOT_RATIONAL_HPP_

#include <istream>      // for basic_istream
#include <ostream>      // for basic_ostream
#include <stdexcept>    // for runtime_error
#include <string>       // for to_string, string_literals
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, is_integral_v, is_signed_v, is_unsigned_v
#include <utility>      // for move

#include "gmpxx.h"
#include "mpfr.h"

#include "integer.hpp"
#include "rational.hpp"

namespace mpfr
{
	template <class Tp>
	inline constexpr bool _mpq_is_other_operands = _mpz_is_other_operands<Tp> || std::is_same_v<Tp, integer>;

	inline integer floor(const rational& z);

	class rational
	{
		friend class integer;
		friend class real;
		friend mpq_srcptr _take(const rational& x) { return x._x; }

	public:
		mpq_t _x;  // NOLINT

		rational() { mpq_init(_x); }
		rational(const rational& o) : rational() { mpq_set(_x, o._x); }
		rational(rational&& o) noexcept : rational() { mpq_swap(_x, o._x); }
		rational& operator=(const rational& o)
		{
			if (this != &o) mpq_set(_x, o._x);
			return *this;
		}
		rational& operator=(rational&& o) noexcept
		{
			if (this != &o) mpq_swap(_x, o._x);
			return *this;
		}
		~rational() { mpq_clear(_x); }
		void swap(rational& o) & { mpq_swap(_x, o._x); }

		[[nodiscard]] rational clone() const { return *this; }
		[[nodiscard]] bool iszero() const { return mpq_sgn(_x) == 0; }
		void negate() & { mpq_neg(_x, _x); }
		void invert() & { mpq_inv(_x, _x); }

		[[nodiscard]] std::string str() const { return std::string(mpq_get_str(nullptr, 10, _x)); }

		template <class T, class = std::enable_if_t<_mpq_is_other_operands<T> || std::is_same_v<T, double>>>
		explicit rational(T o) : rational()
		{
			_mp_ops<T>::set(_x, o);
		}
		template <class T1, class T2,
		          class = std::enable_if_t<std::is_integral_v<T1> && std::is_integral_v<T2> && std::is_unsigned_v<T2>>>
		rational(T1 num, T2 den) : rational()
		{
			if constexpr (std::is_signed_v<T1>)
				mpq_set_si(_x, num, den);
			else
				mpq_set_ui(_x, num, den);
			mpq_canonicalize(_x);
		}

		explicit rational(std::string_view o) : rational()
		{
			using namespace std::string_literals;
			if (auto err = mpq_set_str(_x, o.data(), 10); err == -1)
			{
				throw std::runtime_error("in mpfr::rational(string_view):\n  invalid input format "s += o);
				mpq_clear(_x);
			}
			mpq_canonicalize(_x);
		}

		template <class T, class = std::enable_if_t<_mpq_is_other_operands<T> || std::is_same_v<T, double>>>
		rational& operator=(T o) &
		{
			_mp_ops<T>::set(_x, o);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
			return *this;
#pragma GCC diagnostic pop
		}
		rational& operator=(std::string_view o)
		{
			using namespace std::string_literals;
			if (auto err = mpq_set_str(_x, o.data(), 10); err == -1)
				throw std::runtime_error("in mpfr::rational(string_view):\n  invalid input format "s += o);
			mpq_canonicalize(_x);
			return *this;
		}

		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		explicit operator T() const
		{
			return _mp_ops<T>::get(_x);
		}

		explicit operator integer() const { return floor(*this); }

		// _cmp(a, b) returns the sign of a - b

		friend int _cmp(const rational& r1, const rational& r2) { return mpq_cmp(r1._x, r2._x); }
		template <class T, class = std::enable_if_t<_mpq_is_other_operands<T>>>
		friend int _cmp(const rational& r1, T r2)
		{
			return _mp_ops<T>::cmp(r1._x, r2);
		}
		template <class T, class = std::enable_if_t<_mpq_is_other_operands<T>>>
		friend int _cmp(T r1, const rational& r2)
		{
			return -_cmp(r2, r1);
		}

		template <class Tp>
		inline rational& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_add(_x, _x, o._x);
			else
				_mp_ops<Tp>::addmul(mpq_numref(_x), mpq_denref(_x), o);
			return *this;
		}

		template <class Tp>
		inline rational& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_sub(_x, _x, o._x);
			else
				_mp_ops<Tp>::submul(mpq_numref(_x), mpq_denref(_x), o);
			return *this;
		}

		template <class Tp>
		inline rational& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_mul(_x, _x, o._x);
			else
			{
				_mp_ops<Tp>::mul(mpq_numref(_x), mpq_numref(_x), o);
				mpq_canonicalize(_x);
			}
			return *this;
		}

		template <class Tp>
		inline rational& operator/=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_div(_x, _x, o._x);
			else
			{
				_mp_ops<Tp>::mul(mpq_denref(_x), mpq_denref(_x), o);
				mpq_canonicalize(_x);
			}
			return *this;
		}

		friend rational mul(const rational& r1, const rational& r2) { return r1 * r2; }
		friend rational mul(rational&& r1, const rational& r2) { return std::move(r1) * r2; }
		friend rational mul(const rational& r1, rational&& r2) { return r1 * std::move(r2); }
		friend rational mul(rational&& r1, rational&& r2) { return std::move(r1) * r2; }
		friend rational mul_scalar(const rational& r1, const rational& r2) { return r1 * r2; }
		friend rational mul_scalar(rational&& r1, const rational& r2) { return std::move(r1) * r2; }
		friend rational mul_scalar(const rational& r1, rational&& r2) { return r1 * std::move(r2); }
		friend rational mul_scalar(rational&& r1, rational&& r2) { return std::move(r1) * r2; }

		friend rational operator+(const rational& a, const rational& b)
		{
			rational z;
			mpq_add(z._x, a._x, b._x);
			return z;
		}
		friend rational operator+(rational&& a, const rational& b) { return a += b; }
		friend rational operator+(const rational& a, rational&& b) { return b += a; }
		friend rational operator+(rational&& a, rational&& b) { return a += b; }

		friend rational operator-(const rational& a, const rational& b)
		{
			rational z;
			mpq_sub(z._x, a._x, b._x);
			return z;
		}
		friend rational operator-(rational&& a, const rational& b) { return a -= b; }
		friend rational operator-(const rational& a, rational&& b)
		{
			mpq_sub(b._x, a._x, b._x);
			return std::move(b);
		}
		friend rational operator-(rational&& a, rational&& b) { return a -= b; }

		friend rational operator*(const rational& a, const rational& b)
		{
			rational z;
			mpq_mul(z._x, a._x, b._x);
			return z;
		}
		friend rational operator*(rational&& a, const rational& b) { return a *= b; }
		friend rational operator*(const rational& a, rational&& b) { return b *= a; }
		friend rational operator*(rational&& a, rational&& b) { return a *= b; }

		friend rational operator/(const rational& a, const rational& b)
		{
			rational z;
			mpq_div(z._x, a._x, b._x);
			return z;
		}
		friend rational operator/(rational&& a, const rational& b) { return a /= b; }
		friend rational operator/(const rational& a, rational&& b)
		{
			mpq_div(b._x, a._x, b._x);
			return std::move(b);
		}
		friend rational operator/(rational&& a, rational&& b) { return a /= b; }

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator+(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::add(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator+(rational&& r1, const Tp& r2)
		{
			return r1 += r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator+(const Tp& r1, const rational& r2)
		{
			return r2 + r1;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator+(const Tp& r1, rational&& r2)
		{
			return r2 += r1;
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator-(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::sub_a(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator-(rational&& r1, const Tp& r2)
		{
			return r1 -= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator-(const Tp& r1, const rational& r2)
		{
			rational temp;
			_mp_ops<Tp>::sub_b(temp._x, r1, r2._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator-(const Tp& r1, rational&& r2)
		{
			_mp_ops<Tp>::sub_b(r2._x, r1, r2._x);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator*(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::mul(temp._x, r1._x, r2);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator*(rational&& r1, const Tp& r2)
		{
			return r1 *= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator*(const Tp& r1, const rational& r2) noexcept
		{
			return r2 * r1;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator*(const Tp& r1, rational&& r2)
		{
			return r2 *= r1;
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator/(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::div_a(temp._x, r1._x, r2);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator/(rational&& r1, const Tp& r2)
		{
			return r1 /= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator/(const Tp& r1, const rational& r2)
		{
			rational temp;
			_mp_ops<Tp>::div_b(temp._x, r1, r2._x);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend inline rational operator/(const Tp& r1, rational&& r2)
		{
			_mp_ops<Tp>::div_b(r2._x, r1, r2._x);
			mpq_canonicalize(r2._x);
			return std::move(r2);
		}

		inline rational operator+() const& { return *this; }
		inline rational operator+() && { return std::move(*this); }
		inline rational operator-() const&
		{
			rational temp;
			mpq_neg(temp._x, _x);
			return temp;
		}
		inline rational operator-() &&
		{
			mpq_neg(_x, _x);
			return std::move(*this);
		}

		friend integer floor(const rational& q) { return _mp_ops<integer>::get(q._x); }
		friend inline rational mpfr::_mp_ops<rational>::get(mpfr_srcptr rop, mpfr_rnd_t rnd);

		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const rational& q)
		{
			return s << q._x;
		}
		template <class Char, class Traits>
		friend std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& s, rational& q)
		{
			s >> q._x;
			mpq_canonicalize(q._x);
			return s;
		}
	};
	inline integer mpfr::_mp_ops<integer>::get(mpq_srcptr rop)
	{
		integer z;
		mpz_fdiv_q(z._x, mpq_numref(rop), mpq_denref(rop));
		return z;
	}
	inline rational mpfr::_mp_ops<rational>::get(mpfr_srcptr rop, [[maybe_unused]] mpfr_rnd_t rnd)
	{
		rational q;
		mpfr_get_q(q._x, rop);
		return q;
	}
}  // namespace mpfr

#endif  // QBOOT_RATIONAL_HPP_
