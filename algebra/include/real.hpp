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

/////////////////////////////////////////////////////////////////
// inclusion of headers
/////////////////////////////////////////////////////////////////

#include <cstddef>      // for size_t
#include <cstdint>      // for intmax_t
#include <cstring>      // for strlen
#include <iomanip>      // for setprecision
#include <ios>          // for ios_base, streamsize
#include <iostream>     // for basic_ostram, basic_istream
#include <limits>       // for numeric_limits
#include <locale>       // for use_facet, ctype
#include <sstream>      // for ostringstream
#include <stdexcept>    // for runtime_error
#include <string>       // for string, to_string, basic_string
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, true_type, false_type, is_integral_v, is_signed_v
#include <utility>      // for move

#include "mpfr.h"

namespace mpfr
{
	extern mpfr_prec_t global_prec;
	extern mpfr_rnd_t global_rnd;

	/////////////////////////////////////////////////////////////////
	// type definitions
	/////////////////////////////////////////////////////////////////

	// clang-tidy warns to use int32_t or int64_t instead of long int
	using _long = long int;            // NOLINT
	using _ulong = unsigned long int;  // NOLINT

	class real;

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& helper_ostream(std::basic_ostream<Char, Traits>& s, mpfr_t x,
	                                                        mpfr_rnd_t rnd);

	inline bool _is_nan(mpfr_srcptr x) { return mpfr_nan_p(x) != 0; }    // NOLINT
	inline bool _is_inf(mpfr_srcptr x) { return mpfr_inf_p(x) != 0; }    // NOLINT
	inline bool _is_zero(mpfr_srcptr x) { return mpfr_zero_p(x) != 0; }  // NOLINT
	inline int _sgn(mpfr_srcptr x) { return mpfr_sgn(x); }               // NOLINT
	inline int _signbit(mpfr_srcptr x) { return mpfr_signbit(x); }       // NOLINT

	/////////////////////////////////////////////////////////////////
	// type traits
	/////////////////////////////////////////////////////////////////

	template <class Tp>
	inline constexpr bool is_other_operands =
	    std::is_same_v<Tp, int> || std::is_same_v<Tp, unsigned int> || std::is_same_v<Tp, _long> ||
	    std::is_same_v<Tp, _ulong> || std::is_same_v<Tp, double>;

	template <class Tp>
	struct type_traits;

	template <>
	struct type_traits<_ulong>
	{
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const _ulong op2) { return mpfr_cmp_ui(op1, op2); }
	};

	template <>
	struct type_traits<_long>
	{
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const _long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const _long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const _long op2) { return mpfr_cmp_si(op1, op2); }
	};

	template <>
	struct type_traits<unsigned int> : type_traits<_ulong>
	{
	};

	template <>
	struct type_traits<int> : type_traits<_long>
	{
	};

	template <>
	struct type_traits<double>
	{
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_d(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_d(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_d_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_d(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_d(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_d_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const double op2) { return mpfr_cmp_d(op1, op2); }
	};

	// check all values in integral class I1 are included in integral class I2
	// and I1, I2 has the same signed property
	template <class I1, class I2,
	          class = std::enable_if_t<std::is_integral_v<I1> && std::is_integral_v<I2> &&
	                                   std::is_signed_v<I1> == std::is_signed_v<I2>>>
	inline constexpr bool is_included_v = std::numeric_limits<I2>::min() <= std::numeric_limits<I1>::min() &&
	                                      std::numeric_limits<I1>::max() <= std::numeric_limits<I2>::max();

	/////////////////////////////////////////////////////////////////
	// definitions of relational operators
	/////////////////////////////////////////////////////////////////

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator==(const Tp1& r1, const Tp2& r2)
	{
		if constexpr (std::is_same_v<Tp1, real> && std::is_same_v<Tp2, real>)
			return mpfr_equal_p(r1._x, r2._x) != 0;
		else if constexpr (std::is_same_v<Tp1, real>)
			return type_traits<Tp2>::cmp(r1._x, r2) == 0;
		else
			return r2 == r1;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator!=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 == r2);
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator<(const Tp1& r1, const Tp2& r2)
	{
		if constexpr (std::is_same_v<Tp1, real> && std::is_same_v<Tp2, real>)
			return mpfr_less_p(r1._x, r2._x) != 0;  // NOLINT
		else if constexpr (std::is_same_v<Tp1, real>)
			return type_traits<Tp2>::cmp(r1._x, r2) < 0;
		else
			return type_traits<Tp1>::cmp(r2._x, r1) > 0;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator>(const Tp1& r1, const Tp2& r2)
	{
		return r2 < r1;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator<=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 > r2);
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<std::is_same_v<Tp1, real> || std::is_same_v<Tp2, real>>>
	inline bool operator>=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 < r2);
	}

	/////////////////////////////////////////////////////////////////
	// class definition
	/////////////////////////////////////////////////////////////////

	class real
	{
	public:
		mpfr_t _x;  // NOLINT

		/////////////////////////////////////////////////////////////////
		// default and copy constructors, default assignment operator, destructor
		/////////////////////////////////////////////////////////////////

		// default and copy constructor

		inline real()
		{
			mpfr_init2(_x, global_prec);
			mpfr_set_zero(_x, +1);
		}

		inline real(const real& o)
		{
			mpfr_init2(_x, global_prec);
			mpfr_set(_x, o._x, global_rnd);
		}

		inline real(real&& o) noexcept
		{
			mpfr_init2(_x, global_prec);
			mpfr_swap(_x, o._x);
		}

		// default assignment operator

		inline real& operator=(const real& o) &
		{
			if (this == &o) return *this;
			mpfr_set(_x, o._x, global_rnd);
			return *this;
		}

		inline real& operator=(real&& o) & noexcept
		{
			if (this == &o) return *this;
			mpfr_swap(_x, o._x);
			return *this;
		}

		// destructor

		inline ~real() { mpfr_clear(_x); }

		void swap(real& o) & { mpfr_swap(_x, o._x); }

		[[nodiscard]] real clone() const { return *this; }

		[[nodiscard]] bool iszero() const { return _is_zero(_x); }
		[[nodiscard]] bool isinf() const { return _is_inf(_x); }
		[[nodiscard]] bool isnan() const { return _is_nan(_x); }

		void negate() & { mpfr_mul_si(_x, _x, -1, global_rnd); }

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

		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const real& r)
		{
			return helper_ostream(s, const_cast<mpfr_t&>(r._x), global_rnd);  // NOLINT
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

		explicit inline real(const mpfr_t& o)
		{
			mpfr_init2(_x, global_prec);
			mpfr_set(_x, o, global_rnd);
		}

		explicit inline real(std::string_view op)
		{
			using namespace std::string_literals;
			mpfr_init2(_x, global_prec);
			if (auto err = mpfr_set_str(_x, op.data(), 0, global_rnd); err == -1)
				throw std::runtime_error("in mpfr::real(const char*):\n  invalid input format "s += op);
		}

		explicit inline real(double o)
		{
			mpfr_init2(_x, global_prec);
			mpfr_set_d(_x, o, global_rnd);
		}
		template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
		explicit inline real(T o)
		{
			mpfr_init2(_x, global_prec);
			if constexpr (std::is_signed_v<T>)
			{
				if constexpr (is_included_v<T, _long>)
					mpfr_set_si(_x, o, global_rnd);
				else
					mpfr_set_sj(_x, o, global_rnd);
			}
			else
			{
				if constexpr (is_included_v<T, _ulong>)
					mpfr_set_ui(_x, o, global_rnd);
				else
					mpfr_set_uj(_x, o, global_rnd);
			}
		}

		template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
		inline real(T op, mpfr_exp_t e)
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

		inline real& operator=(double o) &
		{
			mpfr_set_d(_x, o, global_rnd);
			return *this;
		}
		template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
		inline real& operator=(T o) &
		{
			if constexpr (std::is_signed_v<T>)
			{
				if constexpr (is_included_v<T, _long>)
					mpfr_set_si(_x, o, global_rnd);
				else
					mpfr_set_sj(_x, o, global_rnd);
			}
			else
			{
				if constexpr (is_included_v<T, _ulong>)
					mpfr_set_ui(_x, o, global_rnd);
				else
					mpfr_set_uj(_x, o, global_rnd);
			}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
			return *this;
#pragma GCC diagnostic pop
		}

		/////////////////////////////////////////////////////////////////
		// generic operators
		/////////////////////////////////////////////////////////////////

		template <class Tp>
		inline real& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_add(_x, _x, o._x, global_rnd);
			else
				type_traits<Tp>::add(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_sub(_x, _x, o._x, global_rnd);
			else
				type_traits<Tp>::sub_a(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_mul(_x, _x, o._x, global_rnd);
			else
				type_traits<Tp>::mul(_x, _x, o, global_rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator/=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, real>)
				mpfr_div(_x, _x, o._x, global_rnd);
			else
				type_traits<Tp>::div_a(_x, _x, o, global_rnd);
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
		friend real operator+(real&& r1, const real& r2)
		{
			mpfr_add(r1._x, r1._x, r2._x, global_rnd);
			return std::move(r1);
		}
		friend real operator+(const real& r1, real&& r2)
		{
			mpfr_add(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator+(real&& r1, real&& r2) { return std::move(r1) + r2; }

		friend real operator-(const real& r1, const real& r2)
		{
			real temp;
			mpfr_sub(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator-(real&& r1, const real& r2)
		{
			mpfr_sub(r1._x, r1._x, r2._x, global_rnd);
			return std::move(r1);
		}
		friend real operator-(const real& r1, real&& r2)
		{
			mpfr_sub(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator-(real&& r1, real&& r2) { return std::move(r1) - r2; }

		friend real operator*(const real& r1, const real& r2)
		{
			real temp;
			mpfr_mul(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator*(real&& r1, const real& r2)
		{
			mpfr_mul(r1._x, r1._x, r2._x, global_rnd);
			return std::move(r1);
		}
		friend real operator*(const real& r1, real&& r2)
		{
			mpfr_mul(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator*(real&& r1, real&& r2) { return std::move(r1) * r2; }

		friend real operator/(const real& r1, const real& r2)
		{
			real temp;
			mpfr_div(temp._x, r1._x, r2._x, global_rnd);
			return temp;
		}
		friend real operator/(real&& r1, const real& r2)
		{
			mpfr_div(r1._x, r1._x, r2._x, global_rnd);
			return std::move(r1);
		}
		friend real operator/(const real& r1, real&& r2)
		{
			mpfr_div(r2._x, r1._x, r2._x, global_rnd);
			return std::move(r2);
		}
		friend real operator/(real&& r1, real&& r2) { return std::move(r1) / r2; }

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::add(r1._x, r1._x, r2, global_rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r2._x, r1, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::add(r2._x, r2._x, r1, global_rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::sub_a(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::sub_a(r1._x, r1._x, r2, global_rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::sub_b(temp._x, r1, r2._x, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::sub_b(r2._x, r1, r2._x, global_rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::mul(r1._x, r1._x, r2, global_rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const Tp& r1, const real& r2) noexcept
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r2._x, r1, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::mul(r2._x, r2._x, r1, global_rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::div_a(temp._x, r1._x, r2, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::div_a(r1._x, r1._x, r2, global_rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::div_b(temp._x, r1, r2._x, global_rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::div_b(r2._x, r1, r2._x, global_rnd);
			return std::move(r2);
		}

		/////////////////////////////////////////////////////////////////
		// conversion operators
		/////////////////////////////////////////////////////////////////

		explicit inline operator double() const { return mpfr_get_d(_x, global_rnd); }
		template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
		explicit inline operator T() const
		{
			if constexpr (std::is_signed_v<T>)
			{
				if constexpr (is_included_v<T, _long>) return mpfr_get_si(_x, global_rnd);
				return mpfr_get_sj(_x, global_rnd);
			}
			if constexpr (is_included_v<T, _ulong>) return mpfr_get_ui(_x, global_rnd);
			return mpfr_get_uj(_x, global_rnd);
		}

		/////////////////////////////////////////////////////////////////
		// unary operators
		/////////////////////////////////////////////////////////////////

		inline real operator+() const& { return *this; }
		inline real operator+() && { return std::move(*this); }
		inline real operator-() const&
		{
			real temp;
			mpfr_neg(temp._x, _x, global_rnd);
			return temp;
		}
		inline real operator-() &&
		{
			mpfr_neg(_x, _x, global_rnd);
			return std::move(*this);
		}
	};  // class real

	// if sign is nonnegative, return +0
	// if sign is negative, return -0
	inline auto zero(const int n)
	{
		real temp;
		mpfr_set_zero(temp._x, n);
		return temp;
	}

	// if sign is nonnegative, return +inf
	// if sign is negative, return -inf
	inline auto inf(const int n)
	{
		real temp;
		mpfr_set_inf(temp._x, n);
		return temp;
	}

	inline auto nan()
	{
		real temp;
		mpfr_set_nan(temp._x);
		return temp;
	}

	inline auto const_log2()
	{
		real temp;
		mpfr_const_log2(temp._x, global_rnd);
		return temp;
	}

	inline auto const_pi() noexcept
	{
		real temp;
		mpfr_const_pi(temp._x, global_rnd);
		return temp;
	}

	// return n! = Gamma(n + 1)
	inline auto factorial(const _ulong n)
	{
		real temp;
		mpfr_fac_ui(temp._x, n, global_rnd);
		return temp;
	}

	inline auto sqrt(_ulong r)
	{
		real temp;
		mpfr_sqrt_ui(temp._x, r, global_rnd);
		return temp;
	}

	inline auto pow(_ulong op1, _ulong op2)
	{
		real temp;
		mpfr_ui_pow_ui(temp._x, op1, op2, global_rnd);
		return temp;
	}

	inline bool iszero(const real& r) { return _is_zero(r._x); }

	inline int sgn(const real& r) { return _sgn(r._x); }

	inline auto ceil(const real& r)
	{
		real temp;
		mpfr_ceil(temp._x, r._x);
		return temp;
	}
	inline auto ceil(real&& r)
	{
		mpfr_ceil(r._x, r._x);
		return std::move(r);
	}

	inline auto round(const real& r)
	{
		real temp;
		mpfr_round(temp._x, r._x);
		return temp;
	}
	inline auto round(real&& r)
	{
		mpfr_round(r._x, r._x);
		return std::move(r);
	}

	inline auto floor(const real& r)
	{
		real temp;
		mpfr_floor(temp._x, r._x);
		return temp;
	}
	inline auto floor(real&& r)
	{
		mpfr_floor(r._x, r._x);
		return std::move(r);
	}

	inline bool isinteger(const real& r) { return mpfr_integer_p(r._x) != 0; }

	inline auto exp(const real& r)
	{
		real temp;
		mpfr_exp(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto exp(real&& r)
	{
		mpfr_exp(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline auto abs(const real& r)
	{
		real temp;
		mpfr_abs(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto abs(real&& r)
	{
		mpfr_abs(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline auto sin(const real& r)
	{
		real temp;
		mpfr_sin(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto sin(real&& r)
	{
		mpfr_sin(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline auto cos(const real& r)
	{
		real temp;
		mpfr_cos(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto cos(real&& r)
	{
		mpfr_cos(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline auto tan(const real& r)
	{
		real temp;
		mpfr_tan(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto tan(real&& r)
	{
		mpfr_tan(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline auto sec(const real& r)
	{
		real temp;
		mpfr_sec(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto sec(real&& r)
	{
		mpfr_sec(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline auto csc(const real& r)
	{
		real temp;
		mpfr_csc(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto csc(real&& r)
	{
		mpfr_csc(r._x, r._x, global_rnd);
		return std::move(r);
	}
	inline auto cot(const real& r)
	{
		real temp;
		mpfr_cot(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto cot(real&& r)
	{
		mpfr_cot(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline auto log(const real& r)
	{
		real temp;
		mpfr_log(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto log(real&& r)
	{
		mpfr_log(r._x, r._x, global_rnd);
		return std::move(r);
	}

	inline auto sqrt(const real& r)
	{
		real temp;
		mpfr_sqrt(temp._x, r._x, global_rnd);
		return temp;
	}
	inline auto sqrt(real&& r)
	{
		mpfr_sqrt(r._x, r._x, global_rnd);
		return std::move(r);
	}

	/////////////////////////////////////////////////////////////////
	// mathematical functions (definitions for multiple "real" arguments)
	/////////////////////////////////////////////////////////////////

	// operands must have same pricision!
	// if you need to compute f(x, y) where x and y have different precisions, you must cast x or y

	inline auto gamma_inc(const real& r1, const real& r2)
	{
		real temp;
		mpfr_gamma_inc(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}
	inline auto gamma_inc(real&& r1, const real& r2)
	{
		mpfr_gamma_inc(r1._x, r1._x, r2._x, global_rnd);
		return std::move(r1);
	}
	inline auto gamma_inc(const real& r1, real&& r2)
	{
		mpfr_gamma_inc(r2._x, r1._x, r2._x, global_rnd);
		return std::move(r2);
	}
	inline auto gamma_inc(real&& r1, real&& r2) { return gamma_inc(std::move(r1), r2); }

	inline auto fdim(const real& r1, const real& r2)
	{
		real temp;
		mpfr_dim(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline auto fmax(const real& r1, const real& r2)
	{
		real temp;
		mpfr_max(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline auto fmin(const real& r1, const real& r2)
	{
		real temp;
		mpfr_min(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}

	inline auto pow(const real& r1, const real& r2)
	{
		real temp;
		mpfr_pow(temp._x, r1._x, r2._x, global_rnd);
		return temp;
	}
	inline auto pow(real&& r1, const real& r2)
	{
		mpfr_pow(r1._x, r1._x, r2._x, global_rnd);
		return std::move(r1);
	}
	inline auto pow(const real& r1, real&& r2)
	{
		mpfr_pow(r2._x, r1._x, r2._x, global_rnd);
		return std::move(r2);
	}
	inline auto pow(real&& r1, real&& r2) { return pow(std::move(r1), r2); }

	inline auto pow(const real& op1, _ulong op2)
	{
		real temp;
		mpfr_pow_ui(temp._x, op1._x, op2, global_rnd);
		return temp;
	}
	inline auto pow(real&& op1, _ulong op2)
	{
		mpfr_pow_ui(op1._x, op1._x, op2, global_rnd);
		return std::move(op1);
	}
	inline auto pow(const real& op1, unsigned int op2) { return pow(op1, _ulong(op2)); }
	inline auto pow(real&& op1, unsigned int op2) { return pow(std::move(op1), _ulong(op2)); }

	inline auto pow(const real& op1, _long op2)
	{
		real temp;
		mpfr_pow_si(temp._x, op1._x, op2, global_rnd);
		return temp;
	}
	inline auto pow(real&& op1, _long op2)
	{
		mpfr_pow_si(op1._x, op1._x, op2, global_rnd);
		return std::move(op1);
	}
	inline auto pow(const real& op1, int op2) { return pow(op1, _long(op2)); }
	inline auto pow(real&& op1, int op2) { return pow(std::move(op1), _long(op2)); }

	inline auto pow(_ulong op1, const real& op2)
	{
		real temp;
		mpfr_ui_pow(temp._x, op1, op2._x, global_rnd);
		return temp;
	}
	inline auto pow(_ulong op1, real&& op2)
	{
		mpfr_ui_pow(op2._x, op1, op2._x, global_rnd);
		return std::move(op2);
	}

	// returns (x)_(n) = x (x + 1) ... (x + n - 1)
	inline auto pochhammer(const real& x, _ulong n)
	{
		real temp(1);
		for (_ulong i = 0; i < n; ++i) temp *= x + i;
		return temp;
	}
	// do not need to define rvalue version because we cannot reuse x in the above algorithm

	inline int cmpabs(const real& r1, const real& r2) { return mpfr_cmpabs(r1._x, r2._x); }
}  // namespace mpfr

// io functions
namespace mpfr
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
	inline std::basic_ostream<Char, Traits>& helper_ostream_const(std::basic_ostream<Char, Traits>& s,
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
	inline std::basic_ostream<Char, Traits>& helper_ostream(std::basic_ostream<Char, Traits>& s, mpfr_t x,
	                                                        mpfr_rnd_t rnd)
	{
		if (_is_nan(x)) return helper_ostream_const(s, "@NaN@", false);
		if (_is_inf(x)) return helper_ostream_const(s, "@Inf@", _sgn(x) < 0);
		if (_is_zero(x)) return helper_ostream_const(s, "0", _signbit(x) != 0);

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
			helper_ostream_const(s, _as_scientific(std::move(t), exp, s.flags()), is_negative);
		else if (style == _exp_style::DEFAULT_FLOAT)
			helper_ostream_const(s, _as_default(std::move(t), exp, s.flags(), prec), is_negative);
		else
			helper_ostream_const(s, _as_fixed(std::move(t), exp, prec), is_negative);

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
	inline bool helper_extract_float(std::basic_istream<Char, Traits>& s, std::string* num)
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
				if (auto ok = helper_extract_float(in, &num); ok && !num.empty())
				{
					if (auto err = mpfr_set_str(r._x, num.c_str(), 0, global_rnd); err == -1)
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
}  // namespace mpfr

#endif  // QBOOT_REAL_HPP_
