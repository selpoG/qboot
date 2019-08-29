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

// //////////////////////////////////////////////////////////////////
// inclusion of headers
// //////////////////////////////////////////////////////////////////

#include <sstream>      // for ostringstream
#include <string>       // for string
#include <type_traits>  // for void_t, enable_if, enable_if_t, is_same_v, true_type, false_type
#include <utility>      // for move

#include "mpfr.h"

// TODO: implement rvalue reference parameters
// example:
//   real operator+(real&& x, const real& y) {
//     return std::move(x += y);
//   }
//   real exp(real&& x) {
//     mpfr_exp(x._x, x._x, rnd);
//     return std::move(x);
//   }

namespace mpfr
{
	// //////////////////////////////////////////////////////////////////
	// type definitions
	// //////////////////////////////////////////////////////////////////

	using real_prec_t = mpfr_prec_t;
	using real_exp_t = mp_exp_t;
	using real_rnd_t = mpfr_rnd_t;

	using mpfr_old_long = long int;
	using mpfr_old_ulong = unsigned long int;

	template <real_prec_t _prec, real_rnd_t _rnd>
	class real;
}  // namespace mpfr

namespace mpfr
{
	// //////////////////////////////////////////////////////////////////
	// type traits
	// //////////////////////////////////////////////////////////////////

	template <class Tp>
	constexpr bool is_other_operands =
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
	// //////////////////////////////////////////////////////////////////
	// basic meta-programming
	// //////////////////////////////////////////////////////////////////

	// //////////////////////////////////////////////////////////////////
	// calculated result type for two arguments
	// //////////////////////////////////////////////////////////////////

	// At least one argument must be of type "real" and all "real"s must have
	// the same rounding ("_rnd").  The result type is the "real" with highest
	// precision.

	template <class Tp1, class Tp2>
	struct result_type
	{
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class Tp>
	struct result_type<real<_prec, _rnd>, Tp>
	{
		typedef real<_prec, _rnd> type;
		static constexpr real_rnd_t rnd = _rnd;
		static constexpr real_prec_t prec = _prec;
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class Tp>
	struct result_type<Tp, real<_prec, _rnd>>
	{
		typedef real<_prec, _rnd> type;
		static constexpr real_rnd_t rnd = _rnd;
		static constexpr real_prec_t prec = _prec;
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd>
	struct result_type<real<_prec1, _rnd>, real<_prec2, _rnd>>
	{
		typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
		static constexpr real_rnd_t rnd = _rnd;
		static constexpr real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
	};

	template <class Tp1, class Tp2>
	using result_type_t = typename result_type<Tp1, Tp2>::type;

	// T is mpfr::real or not
	template <class T>
	struct is_mpfr_real;
	template <class T>
	inline constexpr bool is_mpfr_real_v = is_mpfr_real<T>::value;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct is_mpfr_real<mpfr::real<prec, rnd>> : std::true_type
	{
	};
	template <class T>
	struct is_mpfr_real : std::false_type
	{
	};
}  // namespace mpfr

namespace mpfr
{
	// //////////////////////////////////////////////////////////////
	// generic operators (definitions of binary operators)
	// //////////////////////////////////////////////////////////////

	template <class Tp1, class Tp2>
	inline result_type_t<Tp1, Tp2> mul(const Tp1& r1, const Tp2& r2)
	{
		return r1 * r2;
	}
	template <class Tp1, class Tp2>
	inline result_type_t<Tp1, Tp2> mul_scalar(const Tp1& r1, const Tp2& r2)
	{
		return r1 * r2;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto operator+(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_add(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto operator-(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_sub(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto operator*(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_mul(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline result_type_t<Tp1, Tp2> operator/(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_div(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

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
			return mpfr_less_p(r1._x, r2._x) != 0;
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

	// //////////////////////////////////////////////////////////////
	// mathematical constants
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto zero(const int n)
	{
		real<_prec, _rnd> temp;
		mpfr_set_zero(temp._x, n);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto inf(const int n)
	{
		real<_prec, _rnd> temp;
		mpfr_set_inf(temp._x, n);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto nan()
	{
		real<_prec, _rnd> temp;
		mpfr_set_nan(temp._x);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto const_log2()
	{
		real<_prec, _rnd> temp;
		mpfr_const_log2(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto const_pi() noexcept
	{
		real<_prec, _rnd> temp;
		mpfr_const_pi(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto factorial(const mpfr_old_ulong n)
	{
		real<_prec, _rnd> temp;
		mpfr_fac_ui(temp._x, n, _rnd);
		return temp;
	}

	// //////////////////////////////////////////////////////////////////
	// class definition
	// //////////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	class real
	{
	public:
		mpfr_t _x;

		static constexpr real_prec_t prec = _prec;
		static constexpr real_rnd_t rnd = _rnd;
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

		// //////////////////////////////////////////////////////////////
		// default and copy constructors, default assignment operator, destructor
		// //////////////////////////////////////////////////////////////

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

		inline real(real&& o)
		{
			mpfr_init2(_x, _prec);
			mpfr_swap(_x, o._x);
		}

		// default assignment operator

		inline real& operator=(const real& o)
		{
			if (&o != this) mpfr_set(_x, o._x, _rnd);
			return *this;
		}

		inline real& operator=(real&& o)
		{
			if (&o != this) mpfr_swap(_x, o._x);
			return *this;
		}

		// destructor

		inline ~real() { mpfr_clear(_x); }

		void swap(real& o) { mpfr_swap(_x, o._x); }

		real clone() const { return *this; }

		bool iszero() const { return mpfr_zero_p(_x); }

		void negate() { mpfr_mul_si(_x, _x, -1, _rnd); }

		real norm() const { return *this * *this; }

		template <class T>
		real eval(const T&) const
		{
			return *this;
		}

		std::string str() const
		{
			std::ostringstream os;
			os << *this;
			return os.str();
		}

		// //////////////////////////////////////////////////////////////
		// converting constructors and converting assignment operators
		// //////////////////////////////////////////////////////////////

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

		template <real_prec_t _prec2, real_rnd_t _rnd2, class = std::enable_if_t<_prec != _prec2>>
		explicit inline real(const real<_prec2, _rnd2>& o)
		{
			mpfr_init2(_x, _prec);
			mpfr_set(_x, o._x, _rnd);
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

		inline real(unsigned long op, mpfr_exp_t e)
		{
			mpfr_init2(_x, _prec);
			mpfr_set_ui_2exp(_x, op, e, _rnd);
		}

		inline real(long op, mpfr_exp_t e)
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
		inline real& operator=(const Tp& o)
		{
			type_traits<Tp>::set(_x, o, _rnd);
			return *this;
		}

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		inline real& operator=(const real<_prec1, _rnd1>& o)
		{
			mpfr_set(_x, o._x, _rnd);
			return *this;
		}

		// //////////////////////////////////////////////////////////////
		// generic operators
		// //////////////////////////////////////////////////////////////

		template <class Tp>
		inline real& operator+=(const Tp& o)
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_add(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::add(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator-=(const Tp& o)
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_sub(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::sub_a(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator*=(const Tp& o)
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_mul(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::mul(_x, _x, o, _rnd);
			return *this;
		}

		template <class Tp>
		inline real& operator/=(const Tp& o)
		{
			if constexpr (is_mpfr_real_v<Tp>)
				mpfr_div(_x, _x, o._x, _rnd);
			else
				type_traits<Tp>::div_a(_x, _x, o, _rnd);
			return *this;
		}

		// //////////////////////////////////////////////////////////////
		// optimized operators
		// //////////////////////////////////////////////////////////////

		template <class Tp, class = std::enable_if_t<!is_mpfr_real_v<Tp>>>
		friend inline real operator+(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator+(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::add(temp._x, r2._x, r1, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::sub_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator-(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::sub_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator*(const Tp& r1, const real& r2) noexcept
		{
			real temp;
			type_traits<Tp>::mul(temp._x, r2._x, r1, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const real& r1, const Tp& r2)
		{
			real temp;
			type_traits<Tp>::div_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class Tp, class = std::enable_if_t<is_other_operands<Tp>>>
		friend inline real operator/(const Tp& r1, const real& r2)
		{
			real temp;
			type_traits<Tp>::div_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}

		// //////////////////////////////////////////////////////////////
		// conversion operators
		// //////////////////////////////////////////////////////////////

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
				mpfr_set(temp._x, _x, _rnd);
				return temp;
			}
			else if constexpr (std::is_same_v<Tp, uintmax_t>)
				return mpfr_get_uj(_x, _rnd);
			else
				return mpfr_get_sj(_x, _rnd);
		}

		// //////////////////////////////////////////////////////////////
		// unary operators
		// //////////////////////////////////////////////////////////////

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

	template <real_prec_t prec, real_rnd_t rnd>
	inline bool iszero(const real<prec, rnd>& r)
	{
		return mpfr_zero_p(r._x);
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline int sgn(const real<prec, rnd>& r)
	{
		return mpfr_sgn(r._x);
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto ceil(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_ceil(temp._x, r._x);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto round(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_round(temp._x, r._x);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto floor(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_floor(temp._x, r._x);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline bool isinteger(const real<prec, rnd>& r)
	{
		return mpfr_integer_p(r._x) != 0;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto exp(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_exp(temp._x, r._x, rnd);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto abs(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_abs(temp._x, r._x, rnd);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto log(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_log(temp._x, r._x, rnd);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto log2(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_log2(temp._x, r._x, rnd);
		return temp;
	}

	template <real_prec_t prec, real_rnd_t rnd>
	inline auto sqrt(const real<prec, rnd>& r)
	{
		real<prec, rnd> temp;
		mpfr_sqrt(temp._x, r._x, rnd);
		return temp;
	}

	// //////////////////////////////////////////////////////////////
	// mathematical functions (definitions for multiple "real" arguments)
	// //////////////////////////////////////////////////////////////

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto gamma_inc(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_gamma_inc(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto fdim(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_dim(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto fmax(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_max(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto fmin(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_min(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline auto pow(const Tp1& r1, const Tp2& r2)
	{
		result_type_t<Tp1, Tp2> temp;
		mpfr_pow(temp._x, r1._x, r2._x, Tp1::rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto pow(const real<_prec, _rnd>& op1, unsigned long op2)
	{
		real<_prec, _rnd> temp;
		mpfr_pow_ui(temp._x, op1._x, op2, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto pow(const real<_prec, _rnd>& op1, long op2)
	{
		real<_prec, _rnd> temp;
		mpfr_pow_si(temp._x, op1._x, op2, _rnd);
		return temp;
	}
	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto pow(unsigned long op1, const real<_prec, _rnd>& op2)
	{
		real<_prec, _rnd> temp;
		mpfr_ui_pow(temp._x, op1, op2._x, _rnd);
		return temp;
	}
	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto pow(unsigned long op1, unsigned long op2)
	{
		real<_prec, _rnd> temp;
		mpfr_ui_pow_ui(temp._x, op1, op2, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline auto pochhammer(const real<_prec, _rnd>& op1, unsigned long op2)
	{
		real<_prec, _rnd> temp(1);
		for (unsigned long i = 0; i < op2; ++i) temp *= op1 + i;
		return temp;
	}

	template <class Tp1, class Tp2, class = std::enable_if_t<is_mpfr_real_v<Tp1> && is_mpfr_real_v<Tp2>>>
	inline int cmpabs(const Tp1& r1, const Tp2& r2)
	{
		return mpfr_cmpabs(r1._x, r2._x);
	}
}  // namespace mpfr

#endif  // QBOOT_REAL_HPP_
