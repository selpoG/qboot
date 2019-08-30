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

#include <algorithm>    // for max
#include <cstddef>      // for size_t
#include <cstdint>      // for intmax_t
#include <cstring>      // for strlen
#include <ios>          // for ios_base, streamsize
#include <iostream>     // for basic_ostram, basic_istream
#include <locale>       // for use_facet, ctype
#include <sstream>      // for ostringstream
#include <stdexcept>    // for runtime_error
#include <string>       // for string, to_string, basic_string
#include <type_traits>  // for void_t, enable_if, enable_if_t, is_same_v, true_type, false_type
#include <utility>      // for move

#include "mpfr.h"

namespace mpfr
{
	/////////////////////////////////////////////////////////////////
	// type definitions
	/////////////////////////////////////////////////////////////////

	// clang-tidy warns to use int32_t or int64_t instead of long int
	using mpfr_old_long = long int;            // NOLINT
	using mpfr_old_ulong = unsigned long int;  // NOLINT

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	class real;

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
	    std::is_same_v<Tp, int> || std::is_same_v<Tp, unsigned int> || std::is_same_v<Tp, mpfr_old_long> ||
	    std::is_same_v<Tp, mpfr_old_ulong> || std::is_same_v<Tp, double>;

	template <class Tp>
	struct type_traits;

	template <>
	struct type_traits<mpfr_old_ulong>
	{
		inline static int set(mpfr_ptr rop, const mpfr_old_ulong op, mpfr_rnd_t rnd)
		{
			return mpfr_set_ui(rop, op, rnd);
		}
		inline static mpfr_old_ulong get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_ui(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const mpfr_old_ulong op2) { return mpfr_cmp_ui(op1, op2); }
	};

	template <>
	struct type_traits<mpfr_old_long>
	{
		inline static int set(mpfr_ptr rop, const mpfr_old_long op, mpfr_rnd_t rnd)
		{
			return mpfr_set_si(rop, op, rnd);
		}
		inline static mpfr_old_long get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_si(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const mpfr_old_long op2) { return mpfr_cmp_si(op1, op2); }
	};

	template <>
	struct type_traits<unsigned int> : type_traits<mpfr_old_ulong>
	{
	};

	template <>
	struct type_traits<int> : type_traits<mpfr_old_long>
	{
	};

	template <>
	struct type_traits<double>
	{
		inline static int set(mpfr_ptr rop, const double op, mpfr_rnd_t rnd) { return mpfr_set_d(rop, op, rnd); }
		inline static double get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_d(op, rnd); }
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
}  // namespace mpfr

namespace mpfr
{
	/////////////////////////////////////////////////////////////////
	// basic meta-programming
	/////////////////////////////////////////////////////////////////
	template <class T>
	struct is_mpfr_real : std::false_type
	{
	};
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct is_mpfr_real<real<prec, rnd>> : std::true_type
	{
	};
	// T is mpfr::real or not
	template <class T>
	inline constexpr bool is_mpfr_real_v = is_mpfr_real<T>::value;
}  // namespace mpfr

namespace mpfr
{
	/////////////////////////////////////////////////////////////////
	// definitions of arithmetic operators
	/////////////////////////////////////////////////////////////////

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		return r1 * r2;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		return std::move(r1) * r2;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		return r1 * std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) * r2;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul_scalar(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		return r1 * r2;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul_scalar(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		return std::move(r1) * r2;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul_scalar(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		return r1 * std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto mul_scalar(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) * r2;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator+(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_add(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator+(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_add(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator+(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_add(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator+(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) + r2;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator-(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_sub(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator-(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_sub(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator-(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_sub(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator-(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) - r2;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator*(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_mul(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator*(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_mul(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator*(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_mul(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator*(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) * r2;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator/(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_div(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator/(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_div(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator/(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_div(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto operator/(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return std::move(r1) / r2;
	}

	/////////////////////////////////////////////////////////////////
	// definitions of relational operators
	/////////////////////////////////////////////////////////////////

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator==(const Tp1& r1, const Tp2& r2)
	{
		if constexpr (is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>)
			return mpfr_equal_p(r1._x, r2._x) != 0;
		else if constexpr (is_mpfr_real_v<Tp1>)
			return type_traits<Tp2>::cmp(r1._x, r2) == 0;
		else
			return r2 == r1;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator!=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 == r2);
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator<(const Tp1& r1, const Tp2& r2)
	{
		if constexpr (is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>)
			return mpfr_less_p(r1._x, r2._x) != 0;  // NOLINT
		else if constexpr (is_mpfr_real_v<Tp1>)
			return type_traits<Tp2>::cmp(r1._x, r2) < 0;
		else
			return type_traits<Tp1>::cmp(r2._x, r1) > 0;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator>(const Tp1& r1, const Tp2& r2)
	{
		return r2 < r1;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator<=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 > r2);
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> || is_mpfr_real_v<Tp2>>>
	inline bool operator>=(const Tp1& r1, const Tp2& r2)
	{
		return !(r1 < r2);
	}

	/////////////////////////////////////////////////////////////////
	// mathematical constants
	/////////////////////////////////////////////////////////////////

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto zero(const int n)
	{
		real<_prec, _rnd> temp;
		mpfr_set_zero(temp._x, n);
		return temp;
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto inf(const int n)
	{
		real<_prec, _rnd> temp;
		mpfr_set_inf(temp._x, n);
		return temp;
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto nan()
	{
		real<_prec, _rnd> temp;
		mpfr_set_nan(temp._x);
		return temp;
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto const_log2()
	{
		real<_prec, _rnd> temp;
		mpfr_const_log2(temp._x, _rnd);
		return temp;
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto const_pi() noexcept
	{
		real<_prec, _rnd> temp;
		mpfr_const_pi(temp._x, _rnd);
		return temp;
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto factorial(const mpfr_old_ulong n)
	{
		real<_prec, _rnd> temp;
		mpfr_fac_ui(temp._x, n, _rnd);
		return temp;
	}

	/////////////////////////////////////////////////////////////////
	// class definition
	/////////////////////////////////////////////////////////////////

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	class real
	{
	public:
		mpfr_t _x;  // NOLINT

		static constexpr mpfr_prec_t prec = _prec;
		static constexpr mpfr_rnd_t rnd = _rnd;
		using type = real;
		// return
		static real log2() { return const_log2<_prec, _rnd>(); }
		static real pi() { return const_pi<_prec, _rnd>(); }
		static real nan() { return mpfr::nan<_prec, _rnd>(); }
		// if sign is nonnegative, return +inf
		// if sign is negative, return -inf
		static real inf(int sign = 1) { return mpfr::inf<_prec, _rnd>(sign); }
		// if sign is nonnegative, return +0
		// if sign is negative, return -0
		static real zero(int sign = 1) { return mpfr::zero<_prec, _rnd>(sign); }
		// return n! = Gamma(n + 1)
		static real factorial(mpfr_old_ulong n) { return mpfr::factorial<_prec, _rnd>(n); }

		/////////////////////////////////////////////////////////////////
		// default and copy constructors, default assignment operator, destructor
		/////////////////////////////////////////////////////////////////

		// default and copy constructor

		inline real()
		{
			mpfr_init2(_x, _prec);
			mpfr_set_zero(_x, +1);
		}

		inline real(const real& o)
		{
			mpfr_init2(_x, _prec);
			mpfr_set(_x, o._x, _rnd);
		}

		inline real(real&& o) noexcept
		{
			mpfr_init2(_x, _prec);
			mpfr_swap(_x, o._x);
		}

		// default assignment operator

		inline real& operator=(const real& o) &
		{
			if (this == &o) return *this;
			mpfr_set(_x, o._x, _rnd);
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

		void negate() & { mpfr_mul_si(_x, _x, -1, _rnd); }

		[[nodiscard]] real norm() const& { return *this * *this; }
		[[nodiscard]] real norm() && { return std::move(*this) * *this; }

		template <class T>
		real eval([[maybe_unused]] const T& x) const
		{
			return *this;
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
			mpfr_init2(_x, _prec);
			mpfr_set(_x, o, _rnd);
		}

		explicit inline real(const char* op)
		{
			mpfr_init2(_x, _prec);
			if (auto err = mpfr_set_str(_x, op, 0, _rnd); err == -1)
				throw std::runtime_error(std::string("in mpfr::real(const char*):\n  invalid input format ") + op);
		}

		explicit inline real(const std::string& op)
		{
			mpfr_init2(_x, _prec);
			if (auto err = mpfr_set_str(_x, op.c_str(), 0, _rnd); err == -1)
				throw std::runtime_error("in mpfr::real(const std::string&):\n  invalid input format " + op);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		explicit inline real(const Tp& o)
		{
			mpfr_init2(_x, _prec);
			type_traits<Tp>::set(_x, o, _rnd);
		}

		template <mpfr_prec_t _prec2, mpfr_rnd_t _rnd2, class = std::enable_if_t<_prec != _prec2 || _rnd != _rnd2>>
		explicit inline real(const real<_prec2, _rnd2>& o)
		{
			mpfr_init2(_x, _prec);
			mpfr_set(_x, o._x, _rnd);
		}
		template <mpfr_prec_t _prec2, mpfr_rnd_t _rnd2, class = std::enable_if_t<_prec != _prec2 || _rnd != _rnd2>>
		explicit inline real(real<_prec2, _rnd2>&& o)
		{
			mpfr_init2(_x, _prec);
			mpfr_swap(_x, o._x);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_ulong, uintmax_t>>>>
		explicit inline real(uintmax_t o)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_uj(_x, o, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_long, intmax_t>>>>
		explicit inline real(intmax_t o)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_sj(_x, o, _rnd);
		}

		inline real(mpfr_old_ulong op, mpfr_exp_t e)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_ui_2exp(_x, op, e, _rnd);
		}

		inline real(mpfr_old_long op, mpfr_exp_t e)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_si_2exp(_x, op, e, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<uintmax_t, mpfr_old_ulong>>>>
		inline real(uintmax_t op, mpfr_exp_t e)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_uj_2exp(_x, op, e, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<intmax_t, mpfr_old_long>>>>
		inline real(intmax_t op, mpfr_exp_t e)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_sj_2exp(_x, op, e, _rnd);
		}

		// converting assignment operators

		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		inline real& operator=(const Tp& o) &
		{
			type_traits<Tp>::set(_x, o, _rnd);
			return *this;
		}

		template <mpfr_prec_t _prec1, mpfr_rnd_t _rnd1, class = std::enable_if_t<_prec != _prec1 || _rnd != _rnd1>>
		inline real& operator=(const real<_prec1, _rnd1>& o) &
		{
			mpfr_set(_x, o._x, _rnd);
			return *this;
		}
		template <mpfr_prec_t _prec1, mpfr_rnd_t _rnd1, class = std::enable_if_t<_prec != _prec1 || _rnd != _rnd1>>
		inline real& operator=(real<_prec1, _rnd1>&& o) &
		{
			mpfr_swap(_x, o._x);
			return *this;
		}

		/////////////////////////////////////////////////////////////////
		// generic operators
		/////////////////////////////////////////////////////////////////

		template <class Tp>
		inline real& operator+=(const Tp& o) &
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_add(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::add(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator-=(const Tp& o) &
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_sub(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::sub_a(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator*=(const Tp& o) &
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_mul(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::mul(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator/=(const Tp& o) &
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_div(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::div_a(_x, _x, o, _rnd);
			return *this;
		}

		/////////////////////////////////////////////////////////////////
		// optimized operators
		/////////////////////////////////////////////////////////////////

		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator+(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r1._x, r2, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator+(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::add(r1._x, r1._x, r2, _rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r2._x, r1, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::add(r2._x, r2._x, r1, _rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::sub_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator-(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::sub_a(r1._x, r1._x, r2, _rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::sub_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::sub_b(r2._x, r1, r2._x, _rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r1._x, r2, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator*(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::mul(r1._x, r1._x, r2, _rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const Tp& r1, const real& r2) noexcept
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r2._x, r1, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::mul(r2._x, r2._x, r1, _rnd);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::div_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator/(real&& r1, const Tp& r2)
		{
			type_traits<Tp>::div_a(r1._x, r1._x, r2, _rnd);
			return std::move(r1);
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::div_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const Tp& r1, real&& r2)
		{
			type_traits<Tp>::div_b(r2._x, r1, r2._x, _rnd);
			return std::move(r2);
		}

		/////////////////////////////////////////////////////////////////
		// conversion operators
		/////////////////////////////////////////////////////////////////

		template <class Tp>
		explicit inline operator Tp() const
		{
			static_assert(is_mpfr_real_v<Tp> || is_other_operands<Tp> || std::is_same_v<Tp, uintmax_t> ||
			              std::is_same_v<Tp, intmax_t>);
			if constexpr (is_other_operands<Tp>)
				return type_traits<Tp>::get(_x, _rnd);
			else if constexpr (is_mpfr_real_v<Tp>)
			{
				Tp temp;
				mpfr_set(temp._x, _x, Tp::rnd);
				return temp;
			}
			else if constexpr (std::is_same_v<Tp, uintmax_t>)
				return mpfr_get_uj(_x, _rnd);
			else
				return mpfr_get_sj(_x, _rnd);
		}

		/////////////////////////////////////////////////////////////////
		// unary operators
		/////////////////////////////////////////////////////////////////

		inline real operator+() const& { return *this; }
		inline real operator+() && { return std::move(*this); }
		inline real operator-() const&
		{
			real<_prec, _rnd> temp;
			mpfr_neg(temp._x, _x, _rnd);
			return temp;
		}
		inline real operator-() &&
		{
			mpfr_neg(_x, _x, _rnd);
			return std::move(*this);
		}
	};  // class real

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline bool iszero(const real<prec, rnd>& r)
	{
		return _is_zero(r._x);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline int sgn(const real<prec, rnd>& r)
	{
		return _sgn(r._x);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto ceil(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_ceil(temp._x, r._x);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto ceil(real<prec, rnd>&& r)
	{
		mpfr_ceil(r._x, r._x);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto round(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_round(temp._x, r._x);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto round(real<prec, rnd>&& r)
	{
		mpfr_round(r._x, r._x);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto floor(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_floor(temp._x, r._x);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto floor(real<prec, rnd>&& r)
	{
		mpfr_floor(r._x, r._x);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline bool isinteger(const real<prec, rnd>& r)
	{
		return mpfr_integer_p(r._x) != 0;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto exp(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_exp(temp._x, r._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto exp(real<prec, rnd>&& r)
	{
		mpfr_exp(r._x, r._x, rnd);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto abs(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_abs(temp._x, r._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto abs(real<prec, rnd>&& r)
	{
		mpfr_abs(r._x, r._x, rnd);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto log(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_log(temp._x, r._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto log(real<prec, rnd>&& r)
	{
		mpfr_log(r._x, r._x, rnd);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto log2(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_log2(temp._x, r._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto log2(real<prec, rnd>&& r)
	{
		mpfr_log2(r._x, r._x, rnd);
		return std::move(r);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto sqrt(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_sqrt(temp._x, r._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto sqrt(real<prec, rnd>&& r)
	{
		mpfr_sqrt(r._x, r._x, rnd);
		return std::move(r);
	}

	/////////////////////////////////////////////////////////////////
	// mathematical functions (definitions for multiple "real" arguments)
	/////////////////////////////////////////////////////////////////

	// operands must have same pricision!
	// if you need to compute f(x, y) where x and y have different precisions, you must cast x or y

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto gamma_inc(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_gamma_inc(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto gamma_inc(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_gamma_inc(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto gamma_inc(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_gamma_inc(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto gamma_inc(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return gamma_inc(std::move(r1), r2);
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto fdim(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_dim(temp._x, r1._x, r2._x, rnd);
		return temp;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto fmax(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_max(temp._x, r1._x, r2._x, rnd);
		return temp;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto fmin(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_min(temp._x, r1._x, r2._x, rnd);
		return temp;
	}

	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto pow(const real<prec, rnd>& r1, const real<prec, rnd>& r2)
	{
		real<prec, rnd> temp;
		mpfr_pow(temp._x, r1._x, r2._x, rnd);
		return temp;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto pow(real<prec, rnd>&& r1, const real<prec, rnd>& r2)
	{
		mpfr_pow(r1._x, r1._x, r2._x, rnd);
		return std::move(r1);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto pow(const real<prec, rnd>& r1, real<prec, rnd>&& r2)
	{
		mpfr_pow(r2._x, r1._x, r2._x, rnd);
		return std::move(r2);
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	inline auto pow(real<prec, rnd>&& r1, real<prec, rnd>&& r2)
	{
		return pow(std::move(r1), r2);
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(const real<_prec, _rnd>& op1, mpfr_old_ulong op2)
	{
		real<_prec, _rnd> temp;
		mpfr_pow_ui(temp._x, op1._x, op2, _rnd);
		return temp;
	}
	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(real<_prec, _rnd>&& op1, mpfr_old_ulong op2)
	{
		mpfr_pow_ui(op1._x, op1._x, op2, _rnd);
		return std::move(op1);
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(const real<_prec, _rnd>& op1, mpfr_old_long op2)
	{
		real<_prec, _rnd> temp;
		mpfr_pow_si(temp._x, op1._x, op2, _rnd);
		return temp;
	}
	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(real<_prec, _rnd>&& op1, mpfr_old_long op2)
	{
		mpfr_pow_si(op1._x, op1._x, op2, _rnd);
		return std::move(op1);
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(mpfr_old_ulong op1, const real<_prec, _rnd>& op2)
	{
		real<_prec, _rnd> temp;
		mpfr_ui_pow(temp._x, op1, op2._x, _rnd);
		return temp;
	}
	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(mpfr_old_ulong op1, real<_prec, _rnd>&& op2)
	{
		mpfr_ui_pow(op2._x, op1, op2._x, _rnd);
		return std::move(op2);
	}

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pow(mpfr_old_ulong op1, mpfr_old_ulong op2)
	{
		real<_prec, _rnd> temp;
		mpfr_ui_pow_ui(temp._x, op1, op2, _rnd);
		return temp;
	}

	// returns (x)_(n) = x (x + 1) ... (x + n - 1)
	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd>
	inline auto pochhammer(const real<_prec, _rnd>& x, mpfr_old_ulong n)
	{
		real<_prec, _rnd> temp(1);
		for (mpfr_old_ulong i = 0; i < n; ++i) temp *= x + i;
		return temp;
	}
	// do not need to define rvalue version because we cannot reuse x in the above algorithm

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline int cmpabs(const Tp1& r1, const Tp2& r2)
	{
		return mpfr_cmpabs(r1._x, r2._x);
	}
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

	template <class Char, class Traits, class String>
	inline std::basic_ostream<Char, Traits>& helper_ostream_const(std::basic_ostream<Char, Traits>& s,
	                                                              const String& abs, bool is_negative)
	{
		auto flags = s.flags();
		auto showpos = (flags & std::ios_base::showpos) != 0;
		auto dir = _get_direction(flags);
		auto len = int(_strlen(abs));
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
				    "in std::ostream& operator<<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
			}
			mpfr_free_str(ch);
			if (exp0 < 0) exp0 = 0;
			exp0 += 3;
		}
		auto ch = mpfr_get_str(nullptr, &exp, 10, prec + size_t(exp0), x, rnd);
		if (is_negative) mpfr_neg(x, x, rnd);
		if (ch == nullptr)
			throw std::runtime_error(
			    "in std::ostream& operator<<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
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

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd, class Char, class Traits>
	std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& in, real<_prec, _rnd>& r)
	{
		try
		{
			typename std::basic_istream<Char, Traits>::sentry s(in, false);
			if (s && !in.eof())
			{
				std::string num;
				if (auto ok = helper_extract_float(in, &num); ok && !num.empty())
				{
					if (auto err = mpfr_set_str(r._x, num.c_str(), 0, _rnd); err == -1)
					{
						in.setstate(std::ios_base::failbit);
						mpfr_set_zero(r._x, +1);
						throw std::runtime_error(
						    "in std::istream& operator>>(std::istream& s, real<_prec, _rnd>& r):\n  invalid input "
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

	template <mpfr_prec_t _prec, mpfr_rnd_t _rnd, class Char, class Traits>
	std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const real<_prec, _rnd>& r)
	{
		return helper_ostream(s, const_cast<mpfr_t&>(r._x), _rnd);  // NOLINT
	}
}  // namespace mpfr

#endif  // QBOOT_REAL_HPP_
