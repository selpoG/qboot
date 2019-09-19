#ifndef QBOOT_RATIONAL_HPP_
#define QBOOT_RATIONAL_HPP_

#include <array>        // for array
#include <cassert>      // for assert
#include <istream>      // for basic_istream
#include <optional>     // for optional
#include <ostream>      // for basic_ostream
#include <stdexcept>    // for runtime_error
#include <string>       // for to_string, string_literals
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, is_integral_v, is_signed_v, is_unsigned_v
#include <utility>      // for move
#include <vector>       // for vector

#include "gmpxx.h"
#include "mpfr.h"

#include "integer.hpp"

namespace mp
{
	template <class Tp>
	inline constexpr bool _mpq_is_other_operands = _mpz_is_other_operands<Tp> || std::is_same_v<Tp, integer>;

	inline integer ceil(const rational& z);
	inline integer floor(const rational& z);
	inline integer truncate(const rational& z);

	class rational
	{
		mpq_t _x;

	public:
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

		explicit rational(mpq_srcptr o) : rational() { mpq_set(_x, o); }

		[[nodiscard]] rational clone() const { return *this; }
		[[nodiscard]] bool iszero() const { return mpq_sgn(_x) == 0; }
		void negate() & { mpq_neg(_x, _x); }
		void invert() & { mpq_inv(_x, _x); }
		// numerator
		integer num() const { return integer(mpq_numref(_x)); }
		// denominator
		integer den() const { return integer(mpq_denref(_x)); }
		// return {numerator, denominator}
		std::array<integer, 2> numden() const { return std::array{integer(mpq_numref(_x)), integer(mpq_denref(_x))}; }
		bool isinteger() const { return mpz_cmp_ui(mpq_denref(_x), 1) == 0; }

		template <class T>
		[[nodiscard]] rational eval([[maybe_unused]] const T& x) const
		{
			return *this;
		}

		[[nodiscard]] std::string str() const { return std::string(mpq_get_str(nullptr, 10, _x)); }
		// replace '/' by slash, '-' by minus
		// for example, rational(-17, 31u).str('#', '%') = "%17#31"
		[[nodiscard]] std::string str(char slash, char minus = '-') const
		{
			std::string s;
			if (*this < 0)
			{
				s += minus;
				s += (-num()).str();
			}
			else
				s += num().str();
			if (!isinteger())
			{
				s += slash;
				s += den().str();
			}
			return s;
		}
		static std::optional<rational> _parse(std::string_view str)
		{
			std::string s(str);
			rational x;
			if (mpq_set_str(x._x, s.data(), 10) == -1) return {};
			mpq_canonicalize(x._x);
			return std::optional{x};
		}

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

		rational(const integer& num, const integer& den) : rational()
		{
			_mp_ops<integer>::set(_x, num, den);
			mpq_canonicalize(_x);
		}

		explicit rational(std::string_view o) : rational()
		{
			std::string s(o);
			using namespace std::string_literals;
			if (mpq_set_str(_x, s.data(), 10) == -1)
			{
				mpq_clear(_x);
				throw std::runtime_error("in mp::rational(string_view):\n  invalid input format "s += o);
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
			std::string s(o);
			using namespace std::string_literals;
			if (mpq_set_str(_x, s.data(), 10) == -1)
				throw std::runtime_error("in mp::rational(string_view):\n  invalid input format "s += o);
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
		rational& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_add(_x, _x, o._x);
			else
				_mp_ops<Tp>::add(_x, _x, o);
			return *this;
		}

		template <class Tp>
		rational& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_sub(_x, _x, o._x);
			else
				_mp_ops<Tp>::sub_a(_x, _x, o);
			return *this;
		}

		template <class Tp>
		rational& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_mul(_x, _x, o._x);
			else
			{
				_mp_ops<Tp>::mul(_x, _x, o);
				mpq_canonicalize(_x);
			}
			return *this;
		}

		template <class Tp>
		rational& operator/=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, rational>)
				mpq_div(_x, _x, o._x);
			else
			{
				_mp_ops<Tp>::div_a(_x, _x, o);
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
		friend rational operator+(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::add(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator+(rational&& r1, const Tp& r2)
		{
			return r1 += r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator+(const Tp& r1, const rational& r2)
		{
			return r2 + r1;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator+(const Tp& r1, rational&& r2)
		{
			return r2 += r1;
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator-(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::sub_a(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator-(rational&& r1, const Tp& r2)
		{
			return r1 -= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator-(const Tp& r1, const rational& r2)
		{
			rational temp;
			_mp_ops<Tp>::sub_b(temp._x, r1, r2._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator-(const Tp& r1, rational&& r2)
		{
			_mp_ops<Tp>::sub_b(r2._x, r1, r2._x);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator*(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::mul(temp._x, r1._x, r2);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator*(rational&& r1, const Tp& r2)
		{
			return r1 *= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator*(const Tp& r1, const rational& r2) noexcept
		{
			return r2 * r1;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator*(const Tp& r1, rational&& r2)
		{
			return r2 *= r1;
		}

		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator/(const rational& r1, const Tp& r2)
		{
			rational temp;
			_mp_ops<Tp>::div_a(temp._x, r1._x, r2);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator/(rational&& r1, const Tp& r2)
		{
			return r1 /= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator/(const Tp& r1, const rational& r2)
		{
			rational temp;
			_mp_ops<Tp>::div_b(temp._x, r1, r2._x);
			mpq_canonicalize(temp._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpq_is_other_operands<Tp>>>
		friend rational operator/(const Tp& r1, rational&& r2)
		{
			_mp_ops<Tp>::div_b(r2._x, r1, r2._x);
			mpq_canonicalize(r2._x);
			return std::move(r2);
		}

		rational operator+() const& { return *this; }
		rational operator+() && { return std::move(*this); }
		rational operator-() const&
		{
			rational temp;
			mpq_neg(temp._x, _x);
			return temp;
		}
		rational operator-() &&
		{
			mpq_neg(_x, _x);
			return std::move(*this);
		}

		friend std::optional<rational> parse(std::string_view str);
		friend integer ceil(const rational& q) { return _mp_ops<integer>::ceil(q._x); }
		friend integer floor(const rational& q) { return _mp_ops<integer>::get(q._x); }
		friend integer truncate(const rational& q) { return _mp_ops<integer>::truncate(q._x); }
		friend mpq_srcptr mp::_mp_ops<rational>::data(const rational& rop);
		friend rational mp::_mp_ops<rational>::get(mpfr_srcptr rop, mpfr_rnd_t rnd);

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
		friend int sgn(const rational& x);
		friend rational abs(const rational& r);
		friend rational abs(rational&& r);
		friend rational pochhammer(const rational& x, _ulong n);
	};
	inline bool isinteger(const rational& r) { return r.isinteger(); }
	inline bool iszero(const rational& r) { return r.iszero(); }
	inline int sgn(const rational& x) { return mpq_sgn(x._x); }
	inline mpq_srcptr _mp_ops<rational>::data(const rational& op) { return op._x; }
	inline rational _mp_ops<rational>::get(mpfr_srcptr rop, [[maybe_unused]] mpfr_rnd_t rnd)
	{
		rational q;
		mpfr_get_q(q._x, rop);
		return q;
	}
	inline rational abs(const rational& r)
	{
		rational temp;
		mpq_abs(temp._x, r._x);
		return temp;
	}
	inline rational abs(rational&& r)
	{
		mpq_abs(r._x, r._x);
		return std::move(r);
	}
	inline rational pochhammer(const rational& x, _ulong n)
	{
		if (n == 0) return rational(1);
		if (n == 1) return x;
		integer xnum = x.num(), xden = x.den();
		integer num(1), den = pow(xden, n);
		for (_ulong i = 0;;)
		{
			num *= xnum;
			if (++i < n)
				xnum += xden;
			else
				break;
		}
		return rational(num, den);
	}
	inline std::vector<integer> continued_fraction(const rational& r)
	{
		if (r < 0)
		{
			auto v = continued_fraction(-r);
			for (auto& x : v) x.negate();
			return v;
		}
		if (r.isinteger()) return {r.num()};
		auto n = mp::floor(r);
		std::vector v{n};
		auto q = 1 / (r - n);
		while (true)
		{
			if (q.isinteger())
			{
				v.push_back(q.num());
				break;
			}
			n = mp::floor(q);
			v.push_back(n);
			q = 1 / (q - n);
		}
		return v;
	}
	// calculate simple rational y which approimates x
	// y satisfies |x - y| < error if not use_relerr
	// if use_relerr, |1 - y / x| < error
	inline rational approx(const rational& x, const rational& error, bool use_relerr = false)
	{
		if (x.isinteger() || error.iszero()) return x;
		if (x < 0) return -approx(-x, error, use_relerr);
		assert(error > 0);
		auto cf = continued_fraction(x);
		auto [a, b] = x.numden();
		// let x = a / b
		// if use_relerr, |1 - p / (q x)| <= error means |a q - p b| <= q a error
		// else, |x - p / q| <= error means |a q - p b| <= q b error
		auto e = (use_relerr ? a : b) * error;
		integer p = cf[0], q(1), p1(1), q1(0);
		for (uint32_t i = 1; i < cf.size() && mp::abs(a * q - p * b) >= q * e; ++i)
		{
			auto np = cf[i] * p + p1;
			auto nq = cf[i] * q + q1;
			p1 = p;
			q1 = q;
			p = np;
			q = nq;
		}
		return {p, q};
	}
	// read str as a rational, even it contains a floating point '.' or an exponential mark 'e' or 'E'.
	// if parse fails, return std::nullopt.
	// str should not contain spaces.
	//   read("123") == rational(123)
	//   read("123/456") == rational(123, 456u) == rational(41, 152u)
	//   read("1.23") == rational(123, 100u)
	//   read("1.25") == rational(125, 100u) == rational(5, 4u)
	//   read("3.2e-8") == rational(1, 31250000u)
	//   read("-5e3") == rational(-5000)
	//   read("12 / 59") -> fail!
	//   read("12/59 ") -> fail!
	inline std::optional<rational> parse(std::string_view str)
	{
		constexpr auto npos = std::string::npos;
		auto t = str.find_first_not_of("+-0123456789/");
		if (t == npos) return rational::_parse(str);
		auto i = str.find_first_of("eE", t);
		if (i == std::string::npos) return _parse_mantisa(str);
		if (str.find_first_not_of("+-0123456789", i + 1) != npos) return {};
		auto mant = _parse_mantisa(str.substr(0, i));
		if (!mant) return {};
		auto exp = integer::_parse(str.substr(i + 1));
		if (!exp) return {};
		if (exp.value() >= 0) return mant.value() * pow(10u, _ulong(exp.value()));
		return mant.value() / pow(10u, _ulong(-exp.value()));
	}
	inline std::optional<rational> _parse_mantisa(std::string_view str)
	{
		constexpr auto npos = std::string::npos;
		auto sg = str.find_first_of("+-");
		if (sg != npos)
		{
			if (sg != 0) return {};
			auto n = _parse_mantisa(str.substr(1));
			if (!n) return {};
			return (str[sg] == '+' ? 1 : -1) * n.value();
		}
		auto i = str.find('.');
		if (i == npos)
		{
			auto x = integer::_parse(str);
			if (!x.has_value()) return {};
			return rational(x.value());
		}
		auto a = i == 0 ? integer() : integer::_parse(str.substr(0, i));
		if (!a) return {};
		if (i + 1 == str.size()) return rational(a.value());
		auto b = integer::_parse(str.substr(i + 1));
		if (!b) return {};
		return a.value() + rational(b.value(), pow(10u, str.size() - i - 1));
	}
}  // namespace mp

#endif  // QBOOT_RATIONAL_HPP_
