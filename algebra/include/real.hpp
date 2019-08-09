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

#ifndef REAL_HPP_
#define REAL_HPP_

// //////////////////////////////////////////////////////////////////
// inclusion of headers
// //////////////////////////////////////////////////////////////////

#include <limits>
#include <sstream>
#include "mpfr.h"

// //////////////////////////////////////////////////////////////////
// default template arguments
// //////////////////////////////////////////////////////////////////

#ifndef MPFR_REAL_CLASS_PREC_DFLT
#define MPFR_REAL_CLASS_PREC_DFLT 53
#endif  // MPFR_REAL_CLASS_PREC_DFLT

#ifndef MPFR_REAL_CLASS_RND_DFLT
#define MPFR_REAL_CLASS_RND_DFLT MPFR_RNDN
#endif  // MPFR_REAL_CLASS_RND_DFLT

// //////////////////////////////////////////////////////////////////
// namespace of MPFR functions
// //////////////////////////////////////////////////////////////////

#ifndef MPFR_NS
#define MPFR_NS
#endif  // MPFR_NS

namespace mpfr
{
	// //////////////////////////////////////////////////////////////////
	// type definitions
	// //////////////////////////////////////////////////////////////////

	typedef mpfr_prec_t real_prec_t;
	typedef mp_exp_t real_exp_t;
	typedef mpfr_rnd_t real_rnd_t;

	using mpfr_old_short = short int;
	using mpfr_old_long = long int;
	using mpfr_old_ushort = unsigned short int;
	using mpfr_old_ulong = unsigned long int;
}  // namespace mpfr

namespace mpfr
{
	// //////////////////////////////////////////////////////////////
	// declaration of real
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	class real;

	// //////////////////////////////////////////////////////////////////
	// type traits
	// //////////////////////////////////////////////////////////////////

	template <class _Tp1, class _Tp2, bool _overwrite>
	struct type_traits
	{
		typedef _Tp1 real_type;
		typedef _Tp2 other_type;

		// support level in class real
		static const bool enable_impl_ctor = false;   // implicit ctor (else expl.)
		static const bool enable_assign_op = false;   // assignment operator
		static const bool enable_conv_func = false;   // conversion member function
		static const bool enable_arithm_ops = false;  // arithmetic operators
		static const bool enable_compar_ops = false;  // comparison operators
		static const bool enable_math_funcs = false;  // mathematical functions

		// support in MPFR library
		// (Note: "has_get_a" beats "has_get_b", if both are "true".)
		static const bool has_set = false;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_old_ulong, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_old_ulong other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const mpfr_old_ulong op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_ui(rop, op, rnd);
		}
		inline static mpfr_old_ulong get_a(mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_ui(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ulong op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_old_long, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_old_long other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const mpfr_old_long op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_si(rop, op, rnd);
		}
		inline static mpfr_old_long get_a(mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_si(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_long op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, unsigned int, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef unsigned int other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const unsigned int op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_ui(rop, op, rnd);
		}
		inline static unsigned int get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<unsigned int>(MPFR_NS mpfr_get_ui(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const unsigned int op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const unsigned int op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, int, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef int other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const int op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_si(rop, op, rnd); }
		inline static int get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<int>(MPFR_NS mpfr_get_si(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const int op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const int op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const int op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_old_ushort, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_old_ushort other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const mpfr_old_ushort op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_ui(rop, op, rnd);
		}
		inline static mpfr_old_ushort get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<mpfr_old_ushort>(MPFR_NS mpfr_get_ui(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ushort op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ushort op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_ushort op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ushort op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_ushort op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_ushort op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_old_short, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_old_short other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const mpfr_old_short op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_si(rop, op, rnd);
		}
		inline static mpfr_old_short get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<mpfr_old_short>(MPFR_NS mpfr_get_si(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_short op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_short op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const mpfr_old_short op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_short op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const mpfr_old_short op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const mpfr_old_short op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, unsigned char, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef unsigned char other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const unsigned char op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_ui(rop, op, rnd);
		}
		inline static unsigned char get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<unsigned char>(MPFR_NS mpfr_get_ui(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const unsigned char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const unsigned char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const unsigned char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_ui_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, signed char, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef signed char other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const signed char op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_si(rop, op, rnd);
		}
		inline static signed char get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<signed char>(MPFR_NS mpfr_get_si(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const signed char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const signed char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const signed char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, char, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef char other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const char op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_si(rop, op, rnd); }
		inline static char get_a(mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return static_cast<char>(MPFR_NS mpfr_get_si(op, rnd));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const char op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const char op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_si_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, double, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef double other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const double op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_d(rop, op, rnd);
		}
		inline static double get_a(mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_d(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_d(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_d(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_d_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_d(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const double op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_d(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, const double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_d_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, float, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef float other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const float op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_flt(rop, op, rnd);
		}
		inline static float get_a(mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_flt(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_d(rop, op1, static_cast<double>(op2), rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_d(rop, op1, static_cast<double>(op2), rnd);
		}
		inline static int sub_b(mpfr_ptr rop, const float op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_d_sub(rop, static_cast<double>(op1), op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_d(rop, op1, static_cast<double>(op2), rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, const float op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_d(rop, op1, static_cast<double>(op2), rnd);
		}
		inline static int div_b(mpfr_ptr rop, const float op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_d_div(rop, static_cast<double>(op1), op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, long double, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef long double other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = true;
		static const bool has_get_b = false;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, const long double op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_set_ld(rop, op, rnd);
		}
		inline static long double get_a(mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_ld(op, rnd); }
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpz_t, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpz_t other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_z(rop, op, rnd); }
		inline static int get_b(mpz_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_z(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_z(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_z(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_z(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_z(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpq_t, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpq_t other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_q(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_q(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_q(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_q(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_q(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpf_t, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpf_t other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_f(rop, op, rnd); }
		inline static int get_b(mpf_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_get_f(rop, op, rnd); }
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_t, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_t other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set(rop, op, rnd); }
		inline static int get_b(mpfr_t rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpz_ptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpz_ptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_z(rop, op, rnd); }
		inline static int get_b(mpz_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_get_z(rop, op, rnd);
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_z(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_z(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_z(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_z(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpq_ptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpq_ptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_q(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_q(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_q(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_q(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_q(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpf_ptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpf_ptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_f(rop, op, rnd); }
		inline static int get_b(mpf_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_get_f(rop, op, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_ptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_ptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = true;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set(rop, op, rnd); }
		inline static int get_b(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpz_srcptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpz_srcptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpz_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_z(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_z(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_z(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_z(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpz_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_z(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpq_srcptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpq_srcptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = false;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpq_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_q(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add_q(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub_q(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul_q(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpq_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div_q(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpf_srcptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpf_srcptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpf_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set_f(rop, op, rnd); }
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, mpfr_srcptr, _overwrite>
	{
		typedef real<_prec, _rnd> real_type;
		typedef mpfr_srcptr other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		static const bool has_set = true;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = true;
		static const bool has_sub_a = true;
		static const bool has_sub_b = true;
		static const bool has_mul = true;
		static const bool has_div_a = true;
		static const bool has_div_b = true;

		// functions in MPFR library
		// (must be defined if corresponding "has_..." boolean is set to "true")
		inline static int set(mpfr_ptr rop, mpfr_srcptr op, mpfr_rnd_t rnd) { return MPFR_NS mpfr_set(rop, op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_add(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_mul(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, mpfr_srcptr op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return MPFR_NS mpfr_div(rop, op1, op2, rnd);
		}
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec1, _rnd>, real<_prec2, _rnd>, _overwrite>
	{
		typedef real<_prec1, _rnd> real_type;
		typedef real<_prec2, _rnd> other_type;

		// support level in class real
		static const bool enable_impl_ctor = true;
		static const bool enable_assign_op = true;
		static const bool enable_conv_func = true;
		static const bool enable_arithm_ops = true;
		static const bool enable_compar_ops = true;
		static const bool enable_math_funcs = true;

		// support in MPFR library
		// (all "has_..." booleans *MUST* be set to "false" for "real<_prec, _rnd>")
		static const bool has_set = false;
		static const bool has_get_a = false;
		static const bool has_get_b = false;
		static const bool has_add = false;
		static const bool has_sub_a = false;
		static const bool has_sub_b = false;
		static const bool has_mul = false;
		static const bool has_div_a = false;
		static const bool has_div_b = false;

		// functions in MPFR library
		// (there should be no function definitions for "real<_prec, _rnd>")
	};

}  // namespace mpfr

namespace mpfr
{
	// //////////////////////////////////////////////////////////////
	// declaration of real
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	class real;

	// //////////////////////////////////////////////////////////////////
	// basic meta-programming
	// //////////////////////////////////////////////////////////////////

	// enable_if

	template <bool, class _Tp = void>
	struct enable_if
	{
	};

	template <class _Tp>
	struct enable_if<true, _Tp>
	{
		typedef _Tp type;
	};

	// //////////////////////////////////////////////////////////////////
	// calculated result type for two arguments
	// //////////////////////////////////////////////////////////////////

	// At least one argument must be of type "real" and all "real"s must have
	// the same rounding ("_rnd").  The result type is the "real" with highest
	// precision.

	template <class _Tp1, class _Tp2, bool _overwrite>
	struct result_type2
	{
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp, bool _overwrite>
	struct result_type2<real<_prec, _rnd>, _Tp, _overwrite>
	{
		typedef real<_prec, _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = _prec;
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp, bool _overwrite>
	struct result_type2<_Tp, real<_prec, _rnd>, _overwrite>
	{
		typedef real<_prec, _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = _prec;
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, bool _overwrite>
	struct result_type2<real<_prec1, _rnd>, real<_prec2, _rnd>, _overwrite>
	{
		typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, bool _overwrite>
	struct result_type2<real<_prec1, _rnd1>, real<_prec2, _rnd2>, _overwrite>
	{
	};

	// //////////////////////////////////////////////////////////////////
	// calculated result type for three arguments
	// //////////////////////////////////////////////////////////////////

	// At least one argument must be of type "real" and all "real"s must have
	// the same rounding ("_rnd").  The result type is the "real" with highest
	// precision.

	template <class _Tp1, class _Tp2, class _Tp3, bool _overwrite>
	struct result_type3
	{
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2, bool _overwrite>
	struct result_type3<real<_prec, _rnd>, _Tp1, _Tp2, _overwrite>
	{
		typedef real<_prec, _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = _prec;
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2, bool _overwrite>
	struct result_type3<_Tp2, real<_prec, _rnd>, _Tp1, _overwrite>
	{
		typedef real<_prec, _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = _prec;
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp1, class _Tp2, bool _overwrite>
	struct result_type3<_Tp1, _Tp2, real<_prec, _rnd>, _overwrite>
	{
		typedef real<_prec, _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = _prec;
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp, bool _overwrite>
	struct result_type3<real<_prec1, _rnd>, real<_prec2, _rnd>, _Tp, _overwrite>
	{
		typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp, bool _overwrite>
	struct result_type3<real<_prec1, _rnd1>, real<_prec2, _rnd2>, _Tp, _overwrite>
	{
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp, bool _overwrite>
	struct result_type3<_Tp, real<_prec1, _rnd>, real<_prec2, _rnd>, _overwrite>
	{
		typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp, bool _overwrite>
	struct result_type3<_Tp, real<_prec1, _rnd1>, real<_prec2, _rnd2>, _overwrite>
	{
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd, class _Tp, bool _overwrite>
	struct result_type3<real<_prec2, _rnd>, _Tp, real<_prec1, _rnd>, _overwrite>
	{
		typedef real<((_prec1 < _prec2) ? _prec2 : _prec1), _rnd> type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec = ((_prec1 < _prec2) ? _prec2 : _prec1);
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2, class _Tp, bool _overwrite>
	struct result_type3<real<_prec2, _rnd2>, _Tp, real<_prec1, _rnd1>, _overwrite>
	{
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_prec_t _prec3, real_rnd_t _rnd, bool _overwrite>
	struct result_type3<real<_prec1, _rnd>, real<_prec2, _rnd>, real<_prec3, _rnd>, _overwrite>
	{
		typedef real<
		    ((((_prec3 < _prec2) ? _prec2 : _prec3) < _prec1) ? _prec1 : ((_prec3 < _prec2) ? _prec2 : _prec3)), _rnd>
		    type;
		static const real_rnd_t rnd = _rnd;
		static const real_prec_t prec =
		    ((((_prec3 < _prec2) ? _prec2 : _prec3) < _prec1) ? _prec1 : ((_prec3 < _prec2) ? _prec2 : _prec3));
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_prec_t _prec3, real_rnd_t _rnd1, real_rnd_t _rnd2,
	          real_rnd_t _rnd3, bool _overwrite>
	struct result_type3<real<_prec1, _rnd1>, real<_prec2, _rnd2>, real<_prec3, _rnd3>, _overwrite>
	{
	};

	// //////////////////////////////////////////////////////////////////
	// promotion to real
	// //////////////////////////////////////////////////////////////////

	template <class _Tp1, class _Tp2>
	struct promote
	{
	};

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
	struct promote<real<_prec, _rnd>, _Tp>
	{
		typedef const real<_prec, _rnd> type;
	};

	template <real_prec_t _prec1, real_prec_t _prec2, real_rnd_t _rnd1, real_rnd_t _rnd2>
	struct promote<real<_prec1, _rnd1>, real<_prec2, _rnd2>>
	{
		typedef const real<_prec2, _rnd2>& type;
	};

	// //////////////////////////////////////////////////////////////////
	// check for equal types
	// //////////////////////////////////////////////////////////////////

	template <class _Tp1, class _Tp2>
	struct equal_types2
	{
		static const bool val = false;
	};

	template <class _Tp>
	struct equal_types2<_Tp, _Tp>
	{
		static const bool val = true;
	};

}  // namespace mpfr

namespace mpfr
{
	// //////////////////////////////////////////////////////////////
	// class declaration
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	class real;

	// //////////////////////////////////////////////////////////////
	// generic operators (definitions of binary operators)
	// //////////////////////////////////////////////////////////////

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
	    const typename result_type2<_Tp1, _Tp2, true>::type>::type
	operator+(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_add(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
	    const typename result_type2<_Tp1, _Tp2, true>::type>::type
	operator-(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_sub(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
	    const typename result_type2<_Tp1, _Tp2, true>::type>::type
	operator*(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_mul(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
	    const typename result_type2<_Tp1, _Tp2, true>::type>::type
	operator/(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_div(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator==(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_equal_p(temp1._x, temp2._x) != 0;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator!=(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_lessgreater_p(temp1._x, temp2._x) != 0;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator<(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_less_p(temp1._x, temp2._x) != 0;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator<=(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_lessequal_p(temp1._x, temp2._x) != 0;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator>(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_greater_p(temp1._x, temp2._x) != 0;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
	    bool>::type
	operator>=(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_greaterequal_p(temp1._x, temp2._x) != 0;
	}

	// //////////////////////////////////////////////////////////////
	// mathematical functions (definitions for single "real" argument)
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	isfinite(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_number_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	isinf(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_inf_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	isnan(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_nan_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	isnormal(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_regular_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	signbit(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_signbit(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	acos(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_acos(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	acosh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_acosh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	asin(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_asin(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	asinh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_asinh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	atan(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_atan(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	atanh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_atanh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cbrt(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_cbrt(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	ceil(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_ceil(temp._x, r._x);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cos(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_cos(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cosh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_cosh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	erf(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_erf(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	erfc(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_erfc(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	exp(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_exp(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	exp2(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_exp2(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	expm1(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_expm1(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	fabs(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_abs(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	abs(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_abs(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	floor(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_floor(temp._x, r._x);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline bool is_integer(const real<_prec, _rnd>& r)
	{
		return r == floor(r);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	log(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_log(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	log10(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_log10(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	log1p(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_log1p(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	log2(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_log2(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	nearbyint(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_rint(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	rint(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_rint(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	round(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_round(temp._x, r._x);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	sin(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_sin(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	sinh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_sinh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	sqrt(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_sqrt(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	tan(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_tan(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	tanh(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_tanh(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	tgamma(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_gamma(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	trunc(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_trunc(temp._x, r._x);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	j0(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_j0(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	j1(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_j1(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	y0(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_y0(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	y1(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_y1(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	ai(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_ai(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cot(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_cot(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	coth(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_coth(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	csc(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_csc(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	csch(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_csch(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	digamma(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_digamma(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	exp10(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_exp10(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	expint(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_eint(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	frac(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_frac(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	isinteger(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_integer_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	iszero(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_zero_p(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	li2(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_li2(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	rec_sqrt(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_rec_sqrt(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	sec(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_sec(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	sech(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_sech(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs, int>::type
	sgn(const real<_prec, _rnd>& r)
	{
		return MPFR_NS mpfr_sgn(r._x);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	zeta(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_zeta(temp._x, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	frexp(const real<_prec, _rnd>& r, real_exp_t* exp)
	{
		if (MPFR_NS mpfr_zero_p(r._x))
		{
			*exp = 0;
			return r;
		}
		else if (MPFR_NS mpfr_inf_p(r._x) || MPFR_NS mpfr_nan_p(r._x))
		{
			// *exp = 0;
			return r;
		}
		else
		{
			real<_prec, _rnd> temp = r;
			*exp = MPFR_NS mpfr_get_exp(r._x);
			MPFR_NS mpfr_set_exp(temp._x, 0);
			return temp;
		}
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real_exp_t>::type
	ilogb(const real<_prec, _rnd>& r)
	{
		if (MPFR_NS mpfr_zero_p(r._x) || MPFR_NS mpfr_nan_p(r._x))
			return std::numeric_limits<real_exp_t>::min();
		else if (MPFR_NS mpfr_inf_p(r._x))
			return std::numeric_limits<real_exp_t>::max();
		else
		{
			real_exp_t temp = MPFR_NS mpfr_get_exp(r._x);
			if (temp != std::numeric_limits<real_exp_t>::min()) temp--;
			return temp;
		}
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	ldexp(const real<_prec, _rnd>& r, const mpfr_old_long exp)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_mul_2si(temp._x, r._x, exp, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	lgamma(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		int signp;
		MPFR_NS mpfr_lgamma(temp._x, &signp, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	logb(const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		if (MPFR_NS mpfr_zero_p(r._x))
			MPFR_NS mpfr_set_inf(temp._x, -1);
		else if (MPFR_NS mpfr_nan_p(r._x))
			MPFR_NS mpfr_set_nan(temp._x);
		else if (MPFR_NS mpfr_inf_p(r._x))
			MPFR_NS mpfr_set_inf(temp._x, 1);
		else
		{
			temp = MPFR_NS mpfr_get_exp(r._x);
			--temp;
		}
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	scalbln(const real<_prec, _rnd>& r, const mpfr_old_long exp)
	{
		return ldexp(r, exp);  // FLT_RADIX == 2???
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	scalbn(const real<_prec, _rnd>& r, const int exp)
	{
		return ldexp(r, static_cast<mpfr_old_long>(exp));  // FLT_RADIX == 2???
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	jn(const mpfr_old_long n, const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_jn(temp._x, n, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	yn(const mpfr_old_long n, const real<_prec, _rnd>& r)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_yn(temp._x, n, r._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cyl_bessel_j(const mpfr_old_long n, const real<_prec, _rnd>& r)
	{
		return jn(n, r);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	cyl_neumann(const mpfr_old_long n, const real<_prec, _rnd>& r)
	{
		return yn(n, r);
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	root(const real<_prec, _rnd>& r, const mpfr_old_ulong n)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_root(temp._x, r._x, n, _rnd);
		return temp;
	}

	// //////////////////////////////////////////////////////////////
	// mathematical functions (definitions for multiple "real" arguments)
	// //////////////////////////////////////////////////////////////

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	isgreater(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_greater_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	isgreaterequal(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_greaterequal_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	isless(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_less_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	islessequal(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_lessequal_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	islessgreater(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_lessgreater_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	isunordered(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_unordered_p(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	atan2(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_atan2(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	copysign(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_copysign(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	fdim(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_dim(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2, class _Tp3>
	inline typename enable_if<
	    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp2, true>::enable_math_funcs &&
	        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp3, true>::enable_math_funcs,
	    typename result_type3<_Tp1, _Tp2, _Tp3, true>::type>::type
	fma(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3)
	{
		typedef typename result_type3<_Tp1, _Tp2, _Tp3, true>::type temp_type;
		const real_rnd_t rnd = result_type3<_Tp1, _Tp2, _Tp3, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		typename promote<temp_type, _Tp3>::type temp3(r3);
		temp_type temp;
		MPFR_NS mpfr_fma(temp._x, temp1._x, temp2._x, temp3._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	fmax(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_max(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	fmin(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_min(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	fmod(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_fmod(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	hypot(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_hypot(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	pow(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_pow(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pow(const real<_prec, _rnd>& op1, unsigned long op2)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_pow_ui(temp._x, op1._x, op2, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pow(const real<_prec, _rnd>& op1, long op2)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_pow_si(temp._x, op1._x, op2, _rnd);
		return temp;
	}
	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pow(const real<_prec, _rnd>& op1, mpz_t op2)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_pow_z(temp._x, op1._x, op2, _rnd);
		return temp;
	}
	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pow(unsigned long op1, const real<_prec, _rnd>& op2)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_ui_pow(temp._x, op1, op2._x, _rnd);
		return temp;
	}
	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pow(unsigned long op1, unsigned long op2)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_ui_pow_ui(temp._x, op1, op2, _rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd>
	inline real<_prec, _rnd> pochhammer(const real<_prec, _rnd>& op1, unsigned long op2)
	{
		real<_prec, _rnd> temp(1);
		for (unsigned long i = 0; i < op2; i++) temp *= op1 + i;
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	remainder(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_remainder(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	agm(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_agm(temp._x, temp1._x, temp2._x, rnd);
		return temp;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    int>::type
	cmpabs(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		return MPFR_NS mpfr_cmpabs(temp1._x, temp2._x);
	}

	template <class _Tp1, class _Tp2, class _Tp3>
	inline typename enable_if<
	    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp2, true>::enable_math_funcs &&
	        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp3, true>::enable_math_funcs,
	    typename result_type3<_Tp1, _Tp2, _Tp3, true>::type>::type
	fms(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3)
	{
		typedef typename result_type3<_Tp1, _Tp2, _Tp3, true>::type temp_type;
		const real_rnd_t rnd = result_type3<_Tp1, _Tp2, _Tp3, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		typename promote<temp_type, _Tp3>::type temp3(r3);
		temp_type temp;
		MPFR_NS mpfr_fms(temp._x, temp1._x, temp2._x, temp3._x, rnd);
		return temp;
	}

	template <real_prec_t _prec, real_rnd_t _rnd, class _Tp>
	inline typename enable_if<
	    type_traits<typename result_type2<real<_prec, _rnd>, _Tp, true>::type, _Tp, true>::enable_math_funcs &&
	        type_traits<typename result_type2<real<_prec, _rnd>, _Tp, true>::type, real<_prec, _rnd>,
	                    true>::enable_math_funcs,
	    typename result_type2<real<_prec, _rnd>, _Tp, true>::type>::type
	modf(const _Tp& r, real<_prec, _rnd>* iptr)
	{
		typedef typename result_type2<real<_prec, _rnd>, _Tp, true>::type temp_type;
		const real_rnd_t rnd = result_type2<real<_prec, _rnd>, _Tp, true>::rnd;
		typename promote<temp_type, _Tp>::type temp(r);
		temp_type temp1;
		MPFR_NS mpfr_modf(iptr->_x, temp1._x, temp._x, rnd);
		return temp1;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	nextafter(const _Tp1& r1, const _Tp2& r2)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		temp_type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		MPFR_NS mpfr_nexttoward(temp1._x, temp2._x);
		return temp1;
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	nexttoward(const _Tp1& r1, const _Tp2& r2)
	{
		return nextafter(r1, r2);
	}

	template <class _Tp1, class _Tp2>
	inline typename enable_if<
	    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
	        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
	    typename result_type2<_Tp1, _Tp2, true>::type>::type
	remquo(const _Tp1& r1, const _Tp2& r2, mpfr_old_long* quo)
	{
		typedef typename result_type2<_Tp1, _Tp2, true>::type temp_type;
		const real_rnd_t rnd = result_type2<_Tp1, _Tp2, true>::rnd;
		typename promote<temp_type, _Tp1>::type temp1(r1);
		typename promote<temp_type, _Tp2>::type temp2(r2);
		temp_type temp;
		MPFR_NS mpfr_remquo(temp._x, quo, temp1._x, temp2._x, rnd);
		return temp;
	}

	// //////////////////////////////////////////////////////////////
	// mathematical constants
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	zero(const int n)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_set_zero(temp._x, n);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	inf(const int n)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_set_inf(temp._x, n);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	nan([[maybe_unused]] const char* tagp)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_set_nan(temp._x);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	const_log2()
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_const_log2(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	const_pi() noexcept
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_const_pi(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	const_euler()
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_const_euler(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	const_catalan()
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_const_catalan(temp._x, _rnd);
		return temp;
	}

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	inline typename enable_if<type_traits<real<_prec, _rnd>, real<_prec, _rnd>, true>::enable_math_funcs,
	                          const real<_prec, _rnd>>::type
	factorial(const mpfr_old_ulong n)
	{
		real<_prec, _rnd> temp;
		MPFR_NS mpfr_fac_ui(temp._x, n, _rnd);
		return temp;
	}

	// //////////////////////////////////////////////////////////////////
	// class definition
	// //////////////////////////////////////////////////////////////////

	template <real_prec_t _prec = MPFR_REAL_CLASS_PREC_DFLT, real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT>
	class real
	{
	public:
		mpfr_t _x;

		static constexpr real_prec_t prec = _prec;
		static constexpr real_rnd_t rnd = _rnd;
		using type = real;

		// //////////////////////////////////////////////////////////////
		// default and copy constructors, default assignment operator, destructor
		// //////////////////////////////////////////////////////////////

		// default and copy constructor

		inline real()
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_zero(_x, +1);
		}

		inline real(const real& o)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set(_x, o._x, _rnd);
		}

		inline real(real&& o)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_swap(_x, o._x);
		}

		// default assignment operator

		inline real& operator=(const real& o)
		{
			if (&o != this) MPFR_NS mpfr_set(_x, o._x, _rnd);
			return *this;
		}

		inline real& operator=(real&& o)
		{
			if (&o != this) MPFR_NS mpfr_swap(_x, o._x);
			return *this;
		}

		// destructor

		inline ~real() { MPFR_NS mpfr_clear(_x); }

		void swap(real& o) { MPFR_NS mpfr_swap(_x, o._x); }

		real clone() const { return *this; }

		bool iszero() const { return mpfr::iszero(*this); }

		real sqrt() const { return mpfr::sqrt(*this); }

		void negate() { MPFR_NS mpfr_mul_si(_x, _x, -1, _rnd); }

		std::string str() const
		{
			std::ostringstream os;
			os << *this;
			return os.str();
		}

		// //////////////////////////////////////////////////////////////
		// converting constructors and converting assignment operators
		// //////////////////////////////////////////////////////////////

		// friend of other reals

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend class real;

		// implicit conversion constructors

		template <class _Tp>
		explicit inline real(const _Tp& o,
		                     [[maybe_unused]]
		                     typename enable_if<type_traits<real, _Tp, true>::has_set &&
		                                        type_traits<real, _Tp, true>::enable_impl_ctor>::type* dummy = nullptr)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			type_traits<real, _Tp, true>::set(_x, o, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_ulong, uintmax_t>>>>
		explicit inline real(uintmax_t o)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_uj(_x, o, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_long, intmax_t>>>>
		explicit inline real(intmax_t o)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_sj(_x, o, _rnd);
		}

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		inline real(const real<_prec1, _rnd1>& o,
		            [[maybe_unused]]
		            typename enable_if<type_traits<real, real<_prec1, _rnd1>, true>::enable_impl_ctor>::type* dummy =
		                nullptr)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set(_x, o._x, _rnd);
		}

		// explicit conversion constructors

		template <class _Tp>
		inline explicit real(
		    const _Tp& o, [[maybe_unused]]
		                  typename enable_if<type_traits<real, _Tp, true>::has_set &&
		                                     (!type_traits<real, _Tp, true>::enable_impl_ctor)>::type* dummy = nullptr)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			type_traits<real, _Tp, true>::set(_x, o, _rnd);
		}

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		inline explicit real(
		    const real<_prec1, _rnd1>& o,
		    [[maybe_unused]]
		    typename enable_if<(!type_traits<real, real<_prec1, _rnd1>, true>::enable_impl_ctor)>::type* dummy =
		        nullptr)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set(_x, o._x, _rnd);
		}

		inline real(unsigned long op, mpfr_exp_t e)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_ui_2exp(_x, op, e, _rnd);
		}
		inline real(long op, mpfr_exp_t e)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_si_2exp(_x, op, e, _rnd);
		}
		template <class = std::void_t<std::enable_if<!std::is_same_v<uintmax_t, mpfr_old_ulong>>>>
		inline real(uintmax_t op, mpfr_exp_t e)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_uj_2exp(_x, op, e, _rnd);
		}
		template <class = std::void_t<std::enable_if<!std::is_same_v<intmax_t, mpfr_old_long>>>>
		inline real(intmax_t op, mpfr_exp_t e)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_sj_2exp(_x, op, e, _rnd);
		}
		inline real(mpz_t op, mpfr_exp_t e)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_z_2exp(_x, op, e, _rnd);
		}

	private:
		// conversion constructor for int
		// ensures a working constructor independent of the configuration in the
		// type_traits of int

		inline real(const int& o, const bool)
		{
			MPFR_NS mpfr_init2(_x, _prec);
			MPFR_NS mpfr_set_si(_x, o, _rnd);
		}

	public:
		// converting assignment operators

		template <class _Tp>
		inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_assign_op && type_traits<real, _Tp, true>::has_set,
		                       real&>::type
		    operator=(const _Tp& o)
		{
			type_traits<real, _Tp, true>::set(_x, o, _rnd);
			return *this;
		}

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		inline typename enable_if<type_traits<real, real<_prec1, _rnd1>, true>::enable_assign_op, real&>::type
		operator=(const real<_prec1, _rnd1>& o)
		{
			MPFR_NS mpfr_set(_x, o._x, _rnd);
			return *this;
		}

		// //////////////////////////////////////////////////////////////
		// generic operators
		// //////////////////////////////////////////////////////////////

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && (!type_traits<real, _Tp, true>::has_add), real&>::type
		operator+=(const _Tp& o)
		{
			typename promote<real, _Tp>::type temp(o);
			MPFR_NS mpfr_add(_x, _x, temp._x, _rnd);
			return *this;
		}

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && (!type_traits<real, _Tp, true>::has_sub_a), real&>::type
		operator-=(const _Tp& o)
		{
			typename promote<real, _Tp>::type temp(o);
			MPFR_NS mpfr_sub(_x, _x, temp._x, _rnd);
			return *this;
		}

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && (!type_traits<real, _Tp, true>::has_mul), real&>::type
		operator*=(const _Tp& o)
		{
			typename promote<real, _Tp>::type temp(o);
			MPFR_NS mpfr_mul(_x, _x, temp._x, _rnd);
			return *this;
		}

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && (!type_traits<real, _Tp, true>::has_div_a), real&>::type
		operator/=(const _Tp& o)
		{
			typename promote<real, _Tp>::type temp(o);
			MPFR_NS mpfr_div(_x, _x, temp._x, _rnd);
			return *this;
		}

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
		    const typename result_type2<_Tp1, _Tp2, true>::type>::type
		operator+(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
		    const typename result_type2<_Tp1, _Tp2, true>::type>::type
		operator-(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
		    const typename result_type2<_Tp1, _Tp2, true>::type>::type
		operator*(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_arithm_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_arithm_ops,
		    const typename result_type2<_Tp1, _Tp2, true>::type>::type
		operator/(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator==(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator!=(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator<(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator<=(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator>(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_compar_ops &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_compar_ops,
		    bool>::type
		operator>=(const _Tp1& r1, const _Tp2& r2);

		// //////////////////////////////////////////////////////////////
		// optimized operators
		// //////////////////////////////////////////////////////////////

		template <class _Tp>
		inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_add,
		                       real&>::type
		    operator+=(const _Tp& o)
		{
			type_traits<real, _Tp, true>::add(_x, _x, o, _rnd);
			return *this;
		}

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_sub_a, real&>::type
		operator-=(const _Tp& o)
		{
			type_traits<real, _Tp, true>::sub_a(_x, _x, o, _rnd);
			return *this;
		}

		template <class _Tp>
		inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_mul,
		                       real&>::type
		    operator*=(const _Tp& o)
		{
			type_traits<real, _Tp, true>::mul(_x, _x, o, _rnd);
			return *this;
		}

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_div_a, real&>::type
		operator/=(const _Tp& o)
		{
			type_traits<real, _Tp, true>::div_a(_x, _x, o, _rnd);
			return *this;
		}

		template <class _Tp>
		friend inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_add,
		                       const real>::type
		    operator+(const real& r1, const _Tp& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::add(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_add,
		                       const real>::type
		    operator+(const _Tp& r1, const real& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::add(temp._x, r2._x, r1, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops &&
		                                     type_traits<real, _Tp, true>::has_sub_a,
		                                 const real>::type
		operator-(const real& r1, const _Tp& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::sub_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops &&
		                                     type_traits<real, _Tp, true>::has_sub_b,
		                                 const real>::type
		operator-(const _Tp& r1, const real& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::sub_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_mul,
		                       const real>::type
		    operator*(const real& r1, const _Tp& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::mul(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline
		    typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops && type_traits<real, _Tp, true>::has_mul,
		                       const real>::type
		    operator*(const _Tp& r1, const real& r2) noexcept
		{
			real temp;
			type_traits<real, _Tp, true>::mul(temp._x, r2._x, r1, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops &&
		                                     type_traits<real, _Tp, true>::has_div_a,
		                                 const real>::type
		operator/(const real& r1, const _Tp& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::div_a(temp._x, r1._x, r2, _rnd);
			return temp;
		}

		template <class _Tp>
		friend inline typename enable_if<type_traits<real, _Tp, true>::enable_arithm_ops &&
		                                     type_traits<real, _Tp, true>::has_div_b,
		                                 const real>::type
		operator/(const _Tp& r1, const real& r2)
		{
			real temp;
			type_traits<real, _Tp, true>::div_b(temp._x, r1, r2._x, _rnd);
			return temp;
		}

		// //////////////////////////////////////////////////////////////
		// conversion operators and functions
		// //////////////////////////////////////////////////////////////

		// conversion operators

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_ulong, uintmax_t>>>>
		explicit inline operator uintmax_t() const
		{
			return MPFR_NS mpfr_get_uj(_x, _rnd);
		}

		template <class = std::void_t<std::enable_if<!std::is_same_v<mpfr_old_long, intmax_t>>>>
		explicit inline operator intmax_t() const
		{
			return MPFR_NS mpfr_get_sj(_x, _rnd);
		}

		explicit inline operator mpfr_old_ulong() const { return MPFR_NS mpfr_get_ui(_x, _rnd); }

		explicit inline operator mpfr_old_long() const { return MPFR_NS mpfr_get_si(_x, _rnd); }

		explicit inline operator unsigned int() const
		{
			return static_cast<unsigned int>(MPFR_NS mpfr_get_ui(_x, _rnd));
		}

		explicit inline operator int() const { return static_cast<int>(MPFR_NS mpfr_get_si(_x, _rnd)); }

		explicit inline operator mpfr_old_ushort() const
		{
			return static_cast<mpfr_old_ushort>(MPFR_NS mpfr_get_ui(_x, _rnd));
		}

		explicit inline operator mpfr_old_short() const
		{
			return static_cast<mpfr_old_short>(MPFR_NS mpfr_get_si(_x, _rnd));
		}

		explicit inline operator unsigned char() const
		{
			return static_cast<unsigned char>(MPFR_NS mpfr_get_ui(_x, _rnd));
		}

		explicit inline operator signed char() const { return static_cast<signed char>(MPFR_NS mpfr_get_si(_x, _rnd)); }

		explicit inline operator char() const { return static_cast<char>(MPFR_NS mpfr_get_si(_x, _rnd)); }

		explicit inline operator double() const { return MPFR_NS mpfr_get_d(_x, _rnd); }

		explicit inline operator float() const { return MPFR_NS mpfr_get_flt(_x, _rnd); }

		explicit inline operator long double() const { return MPFR_NS mpfr_get_ld(_x, _rnd); }

		// conversion functions

		template <class _Tp>
		inline typename enable_if<
		    type_traits<real, _Tp, true>::enable_conv_func && type_traits<real, _Tp, true>::has_get_a, void>::type
		conv(_Tp* o) const
		{
			*o = type_traits<real, _Tp, true>::get_a(_x, _rnd);
		}

		template <class _Tp>
		inline typename enable_if<type_traits<real, _Tp, true>::enable_conv_func &&
		                              (!type_traits<real, _Tp, true>::has_get_a) &&
		                              type_traits<real, _Tp, true>::has_get_b,
		                          void>::type
		conv(_Tp* o) const
		{
			type_traits<real, _Tp, true>::get_b(*o, _x, _rnd);
		}

		// //////////////////////////////////////////////////////////////
		// increment, decrement, and negation operators
		// //////////////////////////////////////////////////////////////

		// increment operators

		inline real& operator++()
		{
			static const real<_prec, _rnd> _one(1, true);
			*this += _one;
			return *this;
		}

		inline const real operator++(int)
		{
			real<_prec, _rnd> temp = *this;
			++(*this);
			return temp;
		}

		// decrement operators

		inline real& operator--()
		{
			static const real<_prec, _rnd> _one(1, true);
			*this -= _one;
			return *this;
		}

		inline const real operator--(int)
		{
			real<_prec, _rnd> temp = *this;
			--(*this);
			return temp;
		}

		// NOTE: The unary member operator- is declared after any template
		// binary friend operator-, because the latter may be unqualified
		// in the code above.  This way we make sure that binary - operations
		// do not match the unary member operator- (in any case).

		inline const real operator-() const
		{
			real<_prec, _rnd> temp;
			MPFR_NS mpfr_neg(temp._x, _x, _rnd);
			return temp;
		}

		// //////////////////////////////////////////////////////////////
		// mathematical functions with single "real" argument
		// //////////////////////////////////////////////////////////////

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		isfinite(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		isinf(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		isnan(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		isnormal(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		signbit(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		acos(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		acosh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		asin(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		asinh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		atan(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		atanh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cbrt(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		ceil(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cos(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cosh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		erf(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		erfc(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		exp(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		exp2(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		expm1(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		fabs(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		abs(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		floor(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		log(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		log10(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		log1p(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		log2(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		nearbyint(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		rint(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		round(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		sin(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		sinh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		sqrt(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		tan(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		tanh(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		tgamma(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		trunc(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		j0(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		j1(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		y0(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		y1(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		ai(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cot(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		coth(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		csc(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		csch(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		digamma(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		exp10(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		expint(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		frac(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		isinteger(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		iszero(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		li2(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		rec_sqrt(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		sec(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		sech(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          int>::type
		sgn(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		zeta(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		frexp(const real<_prec1, _rnd1>& r, real_exp_t* exp);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real_exp_t>::type
		ilogb(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		ldexp(const real<_prec1, _rnd1>& r, const mpfr_old_long exp);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		lgamma(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		logb(const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		scalbln(const real<_prec1, _rnd1>& r, const mpfr_old_long exp);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		scalbn(const real<_prec1, _rnd1>& r, const int exp);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		jn(const mpfr_old_long n, const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		yn(const mpfr_old_long n, const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cyl_bessel_j(const mpfr_old_long n, const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		cyl_neumann(const mpfr_old_long n, const real<_prec1, _rnd1>& r);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		root(const real<_prec1, _rnd1>& r, const mpfr_old_ulong n);

		// //////////////////////////////////////////////////////////////
		// mathematical functions with multiple "real" arguments
		// //////////////////////////////////////////////////////////////

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		isgreater(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		isgreaterequal(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		isless(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		islessequal(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		islessgreater(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		isunordered(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		atan2(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		copysign(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		fdim(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2, class _Tp3>
		friend typename enable_if<
		    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp2, true>::enable_math_funcs &&
		        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp3, true>::enable_math_funcs,
		    typename result_type3<_Tp1, _Tp2, _Tp3, true>::type>::type
		fma(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		fmax(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		fmin(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		fmod(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		hypot(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		pow(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		remainder(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		agm(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    int>::type
		cmpabs(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2, class _Tp3>
		friend typename enable_if<
		    type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp2, true>::enable_math_funcs &&
		        type_traits<typename result_type3<_Tp1, _Tp2, _Tp3, true>::type, _Tp3, true>::enable_math_funcs,
		    typename result_type3<_Tp1, _Tp2, _Tp3, true>::type>::type
		fms(const _Tp1& r1, const _Tp2& r2, const _Tp3& r3);

		template <real_prec_t _prec1, real_rnd_t _rnd1, class _Tp>
		friend typename enable_if<
		    type_traits<typename result_type2<real<_prec1, _rnd1>, _Tp, true>::type, _Tp, true>::enable_math_funcs &&
		        type_traits<typename result_type2<real<_prec1, _rnd1>, _Tp, true>::type, real<_prec1, _rnd1>,
		                    true>::enable_math_funcs,
		    typename result_type2<real<_prec1, _rnd1>, _Tp, true>::type>::type
		modf(const _Tp& r1, real<_prec1, _rnd1>* iptr);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		nextafter(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		nexttoward(const _Tp1& r1, const _Tp2& r2);

		template <class _Tp1, class _Tp2>
		friend typename enable_if<
		    type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp1, true>::enable_math_funcs &&
		        type_traits<typename result_type2<_Tp1, _Tp2, true>::type, _Tp2, true>::enable_math_funcs,
		    typename result_type2<_Tp1, _Tp2, true>::type>::type
		remquo(const _Tp1& r1, const _Tp2& r2, mpfr_old_long* quo);

		// //////////////////////////////////////////////////////////////
		// mathematical constants
		// //////////////////////////////////////////////////////////////

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		zero(const int n);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		inf(const int n);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		nan([[maybe_unused]] const char* tagp);

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		const_log2();

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		const_pi() noexcept;

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		const_euler();

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		const_catalan();

		template <real_prec_t _prec1, real_rnd_t _rnd1>
		friend typename enable_if<type_traits<real<_prec1, _rnd1>, real<_prec1, _rnd1>, true>::enable_math_funcs,
		                          const real<_prec1, _rnd1>>::type
		factorial(const mpfr_old_long n);
	};  // class real
}  // namespace mpfr

#endif  // REAL_HPP_
