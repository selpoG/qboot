//
// SHORT DESCRIPTION
// =================
//
// mpfr::real is a C++ interface to the GNU MPFR library
// version 3.0.0 or later.
//
// COPYRIGHT/LICENSE
// =================
//
// Copyright 2010,2011,2012 Christian Schneider <software(at)chschneider(dot)eu>
//
// Version: 0.0.9-alpha
//
// This file is part of mpfr::real.
//
// mpfr::real is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3 of the License, NOT any later
// version.
//
// mpfr::real is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mpfr::real.  If not, see <http://www.gnu.org/licenses/>.
//

// http://chschneider.eu/programming/mpfr_real/

#ifndef QBOOT_REAL_HPP_
#define QBOOT_REAL_HPP_

#include <cstddef>      // for size_t
#include <cstdint>      // for intmax_t
#include <cstring>      // for strlen
#include <iomanip>      // for setprecision
#include <ios>          // for ios_base, streamsize
#include <iostream>     // for basic_ostram, basic_istream
#include <limits>       // for numeric_limits
#include <locale>       // for use_facet, ctype
#include <optional>     // for optional
#include <sstream>      // for ostringstream
#include <stdexcept>    // for runtime_error
#include <string>       // for string, to_string, basic_string, string_literals
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, is_integral_v, is_signed_v
#include <utility>      // for move

#include "gmpxx.h"
#include "mpfr.h"

#include "integer.hpp"
#include "rational.hpp"

namespace mp
{
	// global variables
	// internal (binary) precision of mp::real
	extern mpfr_prec_t global_prec;
	// rounding mode of mp::real
	extern mpfr_rnd_t global_rnd;

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& _helper_ostream(std::basic_ostream<Char, Traits>& s, mpfr_ptr x,
	                                                         mpfr_rnd_t rnd);

	inline bool _is_nan(mpfr_srcptr x) { return mpfr_nan_p(x) != 0; }    // NOLINT
	inline bool _is_inf(mpfr_srcptr x) { return mpfr_inf_p(x) != 0; }    // NOLINT
	inline bool _is_zero(mpfr_srcptr x) { return mpfr_zero_p(x) != 0; }  // NOLINT
	inline int _sgn(mpfr_srcptr x) { return mpfr_sgn(x); }               // NOLINT
	inline int _signbit(mpfr_srcptr x) { return mpfr_signbit(x); }       // NOLINT

	template <class Tp>
	inline constexpr bool _mpfr_is_other_operands = _mpz_is_other_operands<Tp> || std::is_same_v<Tp, integer> ||
	                                                std::is_same_v<Tp, rational> || std::is_same_v<Tp, double>;

	class real
	{
		mpfr_t _x;

	public:
		/////////////////////////////////////////////////////////////////
		// default and copy constructors, default assignment operator, destructor
		/////////////////////////////////////////////////////////////////

		// default and copy constructor

		real()
		{
			mpfr_init2(_x, global_prec);
			mpfr_set_zero(_x, +1);
		}

		real(const real& o)
		{
			mpfr_init2(_x, global_prec);
			mpfr_set(_x, o._x, global_rnd);
		}

		real(real&& o) noexcept
		{
			mpfr_init2(_x, global_prec);
			mpfr_swap(_x, o._x);
		}

		// default assignment operator

		real& operator=(const real& o) &
		{
			if (this == &o) return *this;
			mpfr_set(_x, o._x, global_rnd);
			return *this;
		}

		real& operator=(real&& o) & noexcept
		{
			if (this == &o) return *this;
			mpfr_swap(_x, o._x);
			return *this;
		}

		// destructor

		~real() { mpfr_clear(_x); }

		void swap(real& o) & { mpfr_swap(_x, o._x); }

		[[nodiscard]] real clone() const { return *this; }

		[[nodiscard]] bool iszero() const { return _is_zero(_x); }
		[[nodiscard]] bool isinf() const { return _is_inf(_x); }
		[[nodiscard]] bool isnan() const { return _is_nan(_x); }

		void negate() & { *this *= -1; }

		[[nodiscard]] real norm() const& { return *this * *this; }
		[[nodiscard]] real norm() && { return std::move(*this) * *this; }

		template <class T>
		[[nodiscard]] real eval([[maybe_unused]] const T& x) const
		{
			return *this;
		}

		[[nodiscard]] std::string str(int32_t precision) const
		{
			std::ostringstream os;
			os << std::setprecision(precision) << *this;
			return os.str();
		}
		static std::optional<real> _parse(std::string_view str)
		{
			std::string s(str);
			mpfr_t x;
			mpfr_init2(x, global_prec);
			if (mpfr_set_str(x, s.data(), 0, global_rnd) == -1) return {};
			return real(x);
		}

		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const real& r)
		{
			return _helper_ostream(s, r.clone()._x, global_rnd);
		}
		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, real&& r)
		{
			return _helper_ostream(s, r._x, global_rnd);
		}

		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << *this;
			return os.str();
		}

		/////////////////////////////////////////////////////////////////
		// converting constructors and converting assignment operators
		/////////////////////////////////////////////////////////////////

		explicit real(mpfr_srcptr o)
		{
			mpfr_init2(_x, global_prec);
			mpfr_set(_x, o, global_rnd);
		}

		explicit real(std::string_view op)
		{
			std::string s(op);
			using namespace std::string_literals;
			mpfr_init2(_x, global_prec);
			if (mpfr_set_str(_x, s.data(), 0, global_rnd) == -1)
			{
				mpfr_clear(_x);
				throw std::runtime_error("in mp::real(string_view):\n  invalid input format "s += op);
			}
		}

		template <class T, class = std::enable_if_t<std::is_integral_v<T> || _mpfr_is_other_operands<T>>>
		explicit real(const T& o)
		{
			mpfr_init2(_x, global_prec);
			if constexpr (_mpfr_is_other_operands<T>)
				_mp_ops<T>::set(_x, o, global_rnd);
			else if constexpr (std::is_signed_v<T>)
				mpfr_set_sj(_x, o, global_rnd);
			else
				mpfr_set_uj(_x, o, global_rnd);
		}

		template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
		real(T op, mpfr_exp_t e)
		{
			mpfr_init2(_x, global_prec);
			if constexpr (std::is_signed_v<T>)
			{
				if constexpr (is_included_v<T, _long>)
					mpfr_set_si_2exp(_x, op, e, global_rnd);
				else
					mpfr_set_sj_2exp(_x, op, e, global_rnd);
			}
			else
			{
				if constexpr (is_included_v<T, _ulong>)
					mpfr_set_ui_2exp(_x, op, e, global_rnd);
				else
					mpfr_set_uj_2exp(_x, op, e, global_rnd);
			}
		}

		// converting assignment operators

		template <class T, class = std::enable_if_t<std::is_integral_v<T> || _mpfr_is_other_operands<T>>>
		real& operator=(const T& o) &
		{
			if constexpr (_mpfr_is_other_operands<T>)
				_mp_ops<T>::set(_x, o, global_rnd);
			else if constexpr (std::is_signed_v<T>)
				mpfr_set_sj(_x, o, global_rnd);
			else
				mpfr_set_uj(_x, o, global_rnd);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
			return *this;
#pragma GCC diagnostic pop
		}

		/////////////////////////////////////////////////////////////////
		// generic operators
		/////////////////////////////////////////////////////////////////

		// _cmp(a, b) returns the sign of a - b

		friend int _cmp(const real& r1, const real& r2) { return mpfr_cmp(r1._x, r2._x); }
		template <class T, class = std::enable_if_t<_mpfr_is_other_operands<T>>>
		friend int _cmp(const real& r1, T r2)
		{
			return _mp_ops<T>::cmp(r1._x, r2);
		}
		template <class T, class = std::enable_if_t<_mpfr_is_other_operands<T>>>
		friend int _cmp(T r1, const real& r2)
		{
			return -_cmp(r2, r1);
		}

		template <class Tp>
		real& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_add(_x, _x, o._x, global_rnd);
			else
				_mp_ops<Tp>::add(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		real& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_sub(_x, _x, o._x, global_rnd);
			else
				_mp_ops<Tp>::sub_a(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		real& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_mul(_x, _x, o._x, global_rnd);
			else
				_mp_ops<Tp>::mul(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		real& operator/=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_div(_x, _x, o._x, global_rnd);
			else
				_mp_ops<Tp>::div_a(_x, _x, o, global_rnd);
			return *this;
		}

		/////////////////////////////////////////////////////////////////
		// optimized operators
		/////////////////////////////////////////////////////////////////

		friend real mul(const real& r1, const real& r2) { return r1 * r2; }
		friend real mul(real&& r1, const real& r2) { return std::move(r1) * r2; }
		friend real mul(const real& r1, real&& r2) { return r1 * std::move(r2); }
		friend real mul(real&& r1, real&& r2) { return std::move(r1) * r2; }
		friend real mul_scalar(const real& r1, const real& r2) { return r1 * r2; }
		friend real mul_scalar(real&& r1, const real& r2) { return std::move(r1) * r2; }
		friend real mul_scalar(const real& r1, real&& r2) { return r1 * std::move(r2); }
		friend real mul_scalar(real&& r1, real&& r2) { return std::move(r1) * r2; }

		friend real operator+(const real& r1, const real& r2)
		{
			real temp;
			mpfr_add(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator+(real&& r1, const real& r2) { return r1 += r2; }
		friend real operator+(const real& r1, real&& r2) { return r2 += r1; }
		friend real operator+(real&& r1, real&& r2) { return r1 += r2; }

		friend real operator-(const real& r1, const real& r2)
		{
			real temp;
			mpfr_sub(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator-(real&& r1, const real& r2) { return r1 -= r2; }
		friend real operator-(const real& r1, real&& r2)
		{
			mpfr_sub(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator-(real&& r1, real&& r2) { return r1 -= r2; }

		friend real operator*(const real& r1, const real& r2)
		{
			real temp;
			mpfr_mul(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator*(real&& r1, const real& r2) { return r1 *= r2; }
		friend real operator*(const real& r1, real&& r2) { return r2 *= r1; }
		friend real operator*(real&& r1, real&& r2) { return r1 *= r2; }

		friend real operator/(const real& r1, const real& r2)
		{
			real temp;
			mpfr_div(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator/(real&& r1, const real& r2) { return r1 /= r2; }
		friend real operator/(const real& r1, real&& r2)
		{
			mpfr_div(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator/(real&& r1, real&& r2) { return r1 /= r2; }

		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator+(const real& r1, const Tp& r2)
		{
			real temp;
			_mp_ops<Tp>::add(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator+(real&& r1, const Tp& r2)
		{
			return r1 += r2;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator+(const Tp& r1, const real& r2)
		{
			return r2 + r1;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator+(const Tp& r1, real&& r2)
		{
			return r2 += r1;
		}

		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator-(const real& r1, const Tp& r2)
		{
			real temp;
			_mp_ops<Tp>::sub_a(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator-(real&& r1, const Tp& r2)
		{
			return r1 -= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator-(const Tp& r1, const real& r2)
		{
			real temp;
			_mp_ops<Tp>::sub_b(temp._x, r1, r2._x, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator-(const Tp& r1, real&& r2)
		{
			_mp_ops<Tp>::sub_b(r2._x, r1, r2._x, global_rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator*(const real& r1, const Tp& r2)
		{
			real temp;
			_mp_ops<Tp>::mul(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator*(real&& r1, const Tp& r2)
		{
			return r1 *= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator*(const Tp& r1, const real& r2) noexcept
		{
			return r2 * r1;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator*(const Tp& r1, real&& r2)
		{
			return r2 *= r1;
		}

		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator/(const real& r1, const Tp& r2)
		{
			real temp;
			_mp_ops<Tp>::div_a(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator/(real&& r1, const Tp& r2)
		{
			return r1 /= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator/(const Tp& r1, const real& r2)
		{
			real temp;
			_mp_ops<Tp>::div_b(temp._x, r1, r2._x, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpfr_is_other_operands<Tp>>>
		friend real operator/(const Tp& r1, real&& r2)
		{
			_mp_ops<Tp>::div_b(r2._x, r1, r2._x, global_rnd);
			return std::move(r2);
		}

		/////////////////////////////////////////////////////////////////
		// conversion operators
		/////////////////////////////////////////////////////////////////

		template <class T, class = std::enable_if_t<std::is_integral_v<T> || _mpfr_is_other_operands<T>>>
		explicit operator T() const
		{
			if constexpr (_mpfr_is_other_operands<T>)
				return _mp_ops<T>::get(_x, global_rnd);
			else if constexpr (std::is_signed_v<T>)
				return mpfr_get_sj(_x, global_rnd);
			else
				return mpfr_get_uj(_x, global_rnd);
		}

		/////////////////////////////////////////////////////////////////
		// unary operators
		/////////////////////////////////////////////////////////////////

		real operator+() const& { return *this; }
		real operator+() && { return std::move(*this); }
		real operator-() const&
		{
			real temp;
			mpfr_neg(temp._x, _x, global_rnd);
			return temp;
		}
		real operator-() &&
		{
			mpfr_neg(_x, _x, global_rnd);
			return std::move(*this);
		}
		friend real zero(const int n);
		friend real inf(const int n);
		friend real nan();
		friend real const_log2();
		friend real const_pi() noexcept;
		friend real sqrt(_ulong r);
		friend bool iszero(const real& r);
		friend int sgn(const real& r);
		friend real ceil(const real& r);
		friend real ceil(real&& r);
		friend real round(const real& r);
		friend real round(real&& r);
		friend real floor(const real& r);
		friend real floor(real&& r);
		friend bool isinteger(const real& r);
		friend real exp(const real& r);
		friend real exp(real&& r);
		friend real abs(const real& r);
		friend real abs(real&& r);
		friend real sin(const real& r);
		friend real sin(real&& r);
		friend real cos(const real& r);
		friend real cos(real&& r);
		friend real tan(const real& r);
		friend real tan(real&& r);
		friend real sec(const real& r);
		friend real sec(real&& r);
		friend real csc(const real& r);
		friend real csc(real&& r);
		friend real cot(const real& r);
		friend real cot(real&& r);
		friend real log(const real& r);
		friend real log(real&& r);
		friend real sqrt(const real& r);
		friend real sqrt(real&& r);
		friend real gamma_inc(const real& r1, const real& r2);
		friend real gamma_inc(real&& r1, const real& r2);
		friend real gamma_inc(const real& r1, real&& r2);
		friend real gamma_inc(real&& r1, real&& r2);
		friend real fdim(const real& r1, const real& r2);
		friend real fmax(const real& r1, const real& r2);
		friend real fmin(const real& r1, const real& r2);
		friend real pow(const real& r1, const real& r2);
		friend real pow(real&& r1, const real& r2);
		friend real pow(const real& r1, real&& r2);
		friend real pow(real&& r1, real&& r2);
		friend real pow(const real& op1, _ulong op2);
		friend real pow(real&& op1, _ulong op2);
		friend real pow(const real& op1, unsigned int op2);
		friend real pow(real&& op1, unsigned int op2);
		friend real pow(const real& op1, _long op2);
		friend real pow(real&& op1, _long op2);
		friend real pow(const real& op1, int op2);
		friend real pow(real&& op1, int op2);
		friend real pow(_ulong op1, const real& op2);
		friend real pow(_ulong op1, real&& op2);
		friend real pochhammer(const real& x, _ulong n);
		friend int cmpabs(const real& r1, const real& r2);
		template <class Char, class Traits>
		friend std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& in, real& r);
	};

	// if sign is nonnegative, return +0
	// if sign is negative, return -0
	inline real zero(const int n)
	{
		real temp;
		mpfr_set_zero(temp._x, n);
		return temp;
	}

	// if sign is nonnegative, return +inf
	// if sign is negative, return -inf
	inline real inf(const int n)
	{
		real temp;
		mpfr_set_inf(temp._x, n);
		return temp;
	}

	inline real nan()
	{
		real temp;
		mpfr_set_nan(temp._x);
		return temp;
	}

	inline real const_log2()
	{
		real temp;
		mpfr_const_log2(temp._x, global_rnd);
		return temp;
	}

	inline real const_pi() noexcept
	{
		real temp;
		mpfr_const_pi(temp._x, global_rnd);
		return temp;
	}

	inline real sqrt(_ulong r)
	{
		real temp;
		mpfr_sqrt_ui(temp._x, r, global_rnd);
		return temp;
	}

	inline bool iszero(const real& r) { return _is_zero(r._x); }

	inline int sgn(const real& r) { return _sgn(r._x); }

	inline real ceil(const real& r)
	{
		real temp;
		mpfr_ceil(temp._x, r._x);
		return temp;
	}
	inline real ceil(real&& r)
	{
		mpfr_ceil(r._x, r._x);
		return std::move(r);
	}

	inline real round(const real& r)
	{
		real temp;
		mpfr_round(temp._x, r._x);
		return temp;
	}
	inline real round(real&& r)
	{
		mpfr_round(r._x, r._x);
		return std::move(r);
	}

	inline real floor(const real& r)
	{
		real temp;
		mpfr_floor(temp._x, r._x);
		return temp;
	}
	inline real floor(real&& r)
	{
		mpfr_floor(r._x, r._x);
		return std::move(r);
	}

	inline bool isinteger(const real& r) { return mpfr_integer_p(r._x) != 0; }

	inline real exp(const real& r)
	{
		real temp;
		mpfr_exp(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real exp(real&& r)
	{
		mpfr_exp(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline real abs(const real& r)
	{
		real temp;
		mpfr_abs(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real abs(real&& r)
	{
		mpfr_abs(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline real sin(const real& r)
	{
		real temp;
		mpfr_sin(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real sin(real&& r)
	{
		mpfr_sin(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline real cos(const real& r)
	{
		real temp;
		mpfr_cos(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real cos(real&& r)
	{
		mpfr_cos(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline real tan(const real& r)
	{
		real temp;
		mpfr_tan(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real tan(real&& r)
	{
		mpfr_tan(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline real sec(const real& r)
	{
		real temp;
		mpfr_sec(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real sec(real&& r)
	{
		mpfr_sec(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline real csc(const real& r)
	{
		real temp;
		mpfr_csc(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real csc(real&& r)
	{
		mpfr_csc(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline real cot(const real& r)
	{
		real temp;
		mpfr_cot(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real cot(real&& r)
	{
		mpfr_cot(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline real log(const real& r)
	{
		real temp;
		mpfr_log(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real log(real&& r)
	{
		mpfr_log(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline real sqrt(const real& r)
	{
		real temp;
		mpfr_sqrt(temp._x, r._x, global_rnd);
		return temp;
	}
	inline real sqrt(real&& r)
	{
		mpfr_sqrt(r._x, r._x, global_rnd);
		return std::move(r);
	}

	/////////////////////////////////////////////////////////////////
	// mathematical functions (definitions for multiple "real" arguments)
	/////////////////////////////////////////////////////////////////

	inline real gamma_inc(const real& r1, const real& r2)
	{
		real temp;
		mpfr_gamma_inc(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}
	inline real gamma_inc(real&& r1, const real& r2)
	{
		mpfr_gamma_inc(r1._x, r1._x, r2._x, global_rnd);
		return std::move(r1);
	}
	inline real gamma_inc(const real& r1, real&& r2)
	{
		mpfr_gamma_inc(r2._x, r1._x, r2._x, global_rnd);
		return std::move(r2);
	}
	inline real gamma_inc(real&& r1, real&& r2) { return gamma_inc(std::move(r1), r2); }

	inline real fdim(const real& r1, const real& r2)
	{
		real temp;
		mpfr_dim(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline real fmax(const real& r1, const real& r2)
	{
		real temp;
		mpfr_max(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline real fmin(const real& r1, const real& r2)
	{
		real temp;
		mpfr_min(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline real pow(const real& r1, const real& r2)
	{
		real temp;
		mpfr_pow(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}
	inline real pow(real&& r1, const real& r2)
	{
		mpfr_pow(r1._x, r1._x, r2._x, global_rnd);
		return std::move(r1);
	}
	inline real pow(const real& r1, real&& r2)
	{
		mpfr_pow(r2._x, r1._x, r2._x, global_rnd);
		return std::move(r2);
	}
	inline real pow(real&& r1, real&& r2) { return pow(std::move(r1), r2); }

	inline real pow(const real& op1, _ulong op2)
	{
		real temp;
		mpfr_pow_ui(temp._x, op1._x, op2, global_rnd);
		return temp;
	}
	inline real pow(real&& op1, _ulong op2)
	{
		mpfr_pow_ui(op1._x, op1._x, op2, global_rnd);
		return std::move(op1);
	}
	inline real pow(const real& op1, unsigned int op2) { return pow(op1, _ulong(op2)); }
	inline real pow(real&& op1, unsigned int op2) { return pow(std::move(op1), _ulong(op2)); }

	inline real pow(const real& op1, _long op2)
	{
		real temp;
		mpfr_pow_si(temp._x, op1._x, op2, global_rnd);
		return temp;
	}
	inline real pow(real&& op1, _long op2)
	{
		mpfr_pow_si(op1._x, op1._x, op2, global_rnd);
		return std::move(op1);
	}
	inline real pow(const real& op1, int op2) { return pow(op1, _long(op2)); }
	inline real pow(real&& op1, int op2) { return pow(std::move(op1), _long(op2)); }

	inline real pow(_ulong op1, const real& op2)
	{
		real temp;
		mpfr_ui_pow(temp._x, op1, op2._x, global_rnd);
		return temp;
	}
	inline real pow(_ulong op1, real&& op2)
	{
		mpfr_ui_pow(op2._x, op1, op2._x, global_rnd);
		return std::move(op2);
	}

	// returns (x)_(n) = x (x + 1) ... (x + n - 1)
	inline real pochhammer(const real& x, _ulong n)
	{
		real temp(1);
		for (_ulong i = 0; i < n; ++i) temp *= x + i;
		return temp;
	}
	// do not need to define rvalue version because we cannot reuse x in the above algorithm

	inline int cmpabs(const real& r1, const real& r2) { return mpfr_cmpabs(r1._x, r2._x); }
}  // namespace mp

// io functions
namespace mp
{
	/////////////////////////////////////////////////////////////////
	// helper functions
	/////////////////////////////////////////////////////////////////

	enum class _pad_direction : int
	{
		RIGHT = 0,
		INTERNAL = 1,
		LEFT = 2
	};
	inline _pad_direction _get_direction(std::ios_base::fmtflags flags)
	{
		if ((flags & std::ios_base::left) != 0) return _pad_direction::LEFT;
		if ((flags & std::ios_base::internal) != 0) return _pad_direction::INTERNAL;
		return _pad_direction::RIGHT;
	}
	enum class _exp_style : int
	{
		DEFAULT_FLOAT = 0,
		SCIENTIFIC = 1,
		FIXED = 2
	};
	inline _exp_style _get_style(std::ios_base::fmtflags flags)
	{
		if ((flags & std::ios_base::scientific) != 0) return _exp_style::SCIENTIFIC;
		if ((flags & std::ios_base::fixed) != 0) return _exp_style::FIXED;
		return _exp_style::DEFAULT_FLOAT;
	}
	inline size_t _needed_precision(_exp_style st, std::streamsize prec)
	{
		if (st != _exp_style::SCIENTIFIC) return prec == 0 ? 1u : static_cast<size_t>(prec);
		return static_cast<size_t>(prec) + 1;
	}
	inline size_t _strlen(const char* s) { return std::strlen(s); }
	inline size_t _strlen(const std::string& s) { return s.length(); }

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& _helper_ostream_const(std::basic_ostream<Char, Traits>& s,
	                                                               std::string_view abs, bool is_negative)
	{
		auto flags = s.flags();
		auto showpos = (flags & std::ios_base::showpos) != 0;
		auto dir = _get_direction(flags);
		auto len = int(abs.length());
		if (showpos || is_negative) ++len;
		auto cnt = int(s.width());
		s.width(0);
		if (dir == _pad_direction::RIGHT)
			for (int i = len; i < cnt; ++i) s << s.fill();
		if (is_negative)
			s << '-';
		else if (showpos)
			s << '+';
		if (dir == _pad_direction::INTERNAL)
			for (int i = len; i < cnt; ++i) s << s.fill();
		s << abs;
		if (dir == _pad_direction::LEFT)
			for (int i = len; i < cnt; ++i) s << s.fill();
		return s;
	}

	inline std::string _as_scientific(std::string&& s, intmax_t exp, std::ios_base::fmtflags flags)
	{
		if (s.length() > 1) s.insert(1, 1, '.');
		--exp;

		if ((flags & std::ios_base::uppercase) != 0)
			s += 'E';
		else
			s += 'e';
		if (exp >= 0)
			s += '+';
		else
		{
			s += '-';
			exp = -exp;
		}
		if (-9 <= exp && exp <= 9) s += '0';
		return s += std::to_string(exp);
	}

	inline std::string _as_default(std::string&& s, intmax_t exp, std::ios_base::fmtflags flags, size_t prec)
	{
		// as fixed
		if (-3 <= exp && exp <= intmax_t(prec))
		{
			if (exp <= 0)
			{
				s.insert(0, size_t(2 - exp), '0');
				s[1] = '.';
			}
			else
				s.insert(size_t(exp), 1, '.');
			while (s.back() == '0') s.pop_back();
			if (s.back() == '.') s.pop_back();
			return std::move(s);
		}
		while (s.back() == '0') s.pop_back();
		return _as_scientific(std::move(s), exp, flags);
	}

	inline void _shrink_to_fit(std::string* s, size_t len)
	{
		if (s->length() > len) s->erase(len, s->length() - len);
	}

	inline std::string _as_fixed(std::string&& s, intmax_t exp, size_t prec)
	{
		if (prec == 0)
		{
			if (exp <= 0) return "0";
			return std::move(s);
		}
		if (exp <= 0)
		{
			s.insert(0, size_t(2 - exp), '0');
			s[1] = '.';
			_shrink_to_fit(&s, 2 + prec);
		}
		else
		{
			s.insert(size_t(exp), 1, '.');
			_shrink_to_fit(&s, size_t(exp) + 1 + prec);
		}
		return std::move(s);
	}

	// TODO(selpo): handle ios_base::showpoint
	// TODO(selpo): handle ios_base::hexfloat

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& _helper_ostream(std::basic_ostream<Char, Traits>& s, mpfr_ptr x,
	                                                         mpfr_rnd_t rnd)
	{
		if (_is_nan(x)) return _helper_ostream_const(s, "@NaN@", false);
		if (_is_inf(x)) return _helper_ostream_const(s, "@Inf@", _sgn(x) < 0);
		if (_is_zero(x)) return _helper_ostream_const(s, "0", _signbit(x) != 0);

		auto style = _get_style(s.flags());
		auto prec = _needed_precision(style, s.precision());

		bool is_negative = _sgn(x) < 0;
		if (is_negative) mpfr_neg(x, x, rnd);
		mpfr_exp_t exp, exp0 = 0;
		if (style == _exp_style::FIXED)
		{
			auto ch = mpfr_get_str(nullptr, &exp0, 10, 1, x, rnd);
			if (ch == nullptr)
			{
				if (is_negative) mpfr_neg(x, x, rnd);
				throw std::runtime_error(
				    "in std::ostream& operator<<(std::ostream& s, const real& r):\n  conversion failed");
			}
			mpfr_free_str(ch);
			if (exp0 < 0) exp0 = 0;
			exp0 += 3;
		}
		auto ch = mpfr_get_str(nullptr, &exp, 10, prec + size_t(exp0), x, rnd);
		if (is_negative) mpfr_neg(x, x, rnd);
		if (ch == nullptr)
			throw std::runtime_error(
			    "in std::ostream& operator<<(std::ostream& s, const real& r):\n  conversion failed");
		std::string t(ch);
		mpfr_free_str(ch);

		if (style == _exp_style::SCIENTIFIC)
			_helper_ostream_const(s, _as_scientific(std::move(t), exp, s.flags()), is_negative);
		else if (style == _exp_style::DEFAULT_FLOAT)
			_helper_ostream_const(s, _as_default(std::move(t), exp, s.flags(), prec), is_negative);
		else
			_helper_ostream_const(s, _as_fixed(std::move(t), exp, prec), is_negative);

		return s;
	}

	enum class _state : int
	{
		ERROR = -2,
		CANCEL = -1,
		BEGIN = 0,
		SIGN,
		POINT,
		EXP,
		END
	};

	template <class Char, class Traits>
	inline bool _helper_extract_float(std::basic_istream<Char, Traits>& s, std::string* num)
	{
		const auto& fac = std::use_facet<std::ctype<Char>>(s.getloc());
		const char pos = fac.widen('+'), neg = fac.widen('-'), point = fac.widen('.'), exp_e = fac.widen('e'),
		           exp_E = fac.widen('E');

		bool has_digits = false, has_edigits = false;
		auto state = _state::BEGIN;

		auto next = [pos, neg, point, exp_e, exp_E](_state st, Char c) -> _state {
			if (c == point) return int(st) < int(_state::POINT) ? _state::POINT : _state::CANCEL;
			if (c == exp_e || c == exp_E) return int(st) < int(_state::EXP) ? _state::EXP : _state::ERROR;
			if ((c == pos || c == neg) && st != _state::BEGIN && st != _state::EXP) return _state::CANCEL;
			switch (st)
			{
			case _state::BEGIN:
			case _state::SIGN: return _state::SIGN;
			case _state::POINT: return _state::POINT;
			case _state::EXP:
			case _state::END: return _state::END;
			case _state::ERROR:
			case _state::CANCEL:
			default: return _state::ERROR;
			}
		};

		// [+-]? [0-9]* [.]? [0-9]* ( [eE] [+-]? [0-9]* )?
		Char c;
		while (s.get(c))
		{
			auto isdigit = fac.is(fac.digit, c);
			if (c != pos && c != neg && !isdigit && c != exp_e && c != exp_E && c != point)
				state = _state::CANCEL;
			else
			{
				if (int(state) < int(_state::EXP))
					has_digits |= isdigit;
				else
					has_edigits |= isdigit;
				state = next(state, c);
			}
			if (int(state) < 0)
			{
				s.putback(c);
				break;
			}
			num->push_back(fac.narrow(c, '\0'));
		}
		return state != _state::ERROR && has_digits && (int(state) < int(_state::EXP) || has_edigits);
	}

	/////////////////////////////////////////////////////////////////
	// std::istream and std::ostream operators
	/////////////////////////////////////////////////////////////////

	template <class Char, class Traits>
	std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& in, real& r)
	{
		try
		{
			typename std::basic_istream<Char, Traits>::sentry s(in, false);
			if (s && !in.eof())
			{
				std::string num;
				if (auto ok = _helper_extract_float(in, &num); ok && !num.empty())
				{
					if (mpfr_set_str(r._x, num.c_str(), 0, global_rnd) == -1)
					{
						in.setstate(std::ios_base::failbit);
						mpfr_set_zero(r._x, +1);
						throw std::runtime_error(
						    "in std::istream& operator>>(std::istream& s, real& r):\n  invalid input "
						    "format " +
						    num);
					}
				}
				else
				{
					in.setstate(std::ios_base::failbit);
					mpfr_set_zero(r._x, +1);
				}
			}
		}
		catch (std::ios_base::failure&)
		{
			mpfr_set_zero(r._x, +1);
			throw;
		}
		return in;
	}
}  // namespace mp

#endif  // QBOOT_REAL_HPP_
