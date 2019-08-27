#ifndef REAL_IO_HPP_
#define REAL_IO_HPP_

#include "real.hpp"

#include <algorithm>  // for max
#include <cstddef>    // for size_t
#include <cstdint>    // for intmax_t
#include <cstring>    // for strlen
#include <ios>        // for ios_base, streamsize
#include <iostream>   // for basic_ostram, basic_istream
#include <locale>     // for use_facet, ctype
#include <stdexcept>  // for runtime_error
#include <string>     // for string, to_string, basic_string
#include <utility>    // for move

namespace mpfr
{
	// //////////////////////////////////////////////////////////////////
	// helper functions
	// //////////////////////////////////////////////////////////////////

	inline int helper_set_stdstr(mpfr_ptr rop, const std::string& op, mpfr_rnd_t rnd) noexcept
	{
		return MPFR_NS mpfr_set_str(rop, op.c_str(), 0, rnd);
	}

	inline int helper_set_charptr(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) noexcept(false)
	{
		const int err = MPFR_NS mpfr_set_str(rop, op, 0, rnd);
		if (err == -1)
			throw std::runtime_error(
			    std::string(
			        "in mpfr::helper_set_charptr(mpfr_ptr, const char*, mpfr_rnd_t):\n  invalid input format ") +
			    op);
		return err;
	}

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

	inline void _shrink_to_fit(std::string& s, size_t len)
	{
		if (s.length() > len) s.erase(len, s.length() - len);
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
			_shrink_to_fit(s, 2 + prec);
		}
		else
		{
			s.insert(size_t(exp), 1, '.');
			_shrink_to_fit(s, size_t(exp) + 1 + prec);
		}
		return std::move(s);
	}

	// missing: handling ios_base::showpoint

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& helper_ostream(std::basic_ostream<Char, Traits>& s, mpfr_t x,
	                                                        mpfr_rnd_t rnd)
	{
		if (MPFR_NS mpfr_nan_p(x)) return helper_ostream_const(s, "@NaN@", false);
		if (MPFR_NS mpfr_inf_p(x)) return helper_ostream_const(s, "@Inf@", MPFR_NS mpfr_sgn(x) < 0);
		if (MPFR_NS mpfr_zero_p(x)) return helper_ostream_const(s, "0", MPFR_NS mpfr_signbit(x) != 0);

		auto style = _get_style(s.flags());
		auto prec = _needed_precision(style, s.precision());

		bool is_negative = MPFR_NS mpfr_sgn(x) < 0;
		if (is_negative) MPFR_NS mpfr_neg(x, x, rnd);
		real_exp_t exp, exp0 = 0;
		if (style == _exp_style::FIXED)
		{
			auto ch = MPFR_NS mpfr_get_str(nullptr, &exp0, 10, 1, x, rnd);
			if (ch == nullptr)
			{
				if (is_negative) MPFR_NS mpfr_neg(x, x, rnd);
				throw std::runtime_error(
				    "in std::ostream& operator<<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
			}
			MPFR_NS mpfr_free_str(ch);
			if (exp0 < 0) exp0 = 0;
			exp0 += 3;
		}
		auto ch = MPFR_NS mpfr_get_str(nullptr, &exp, 10, prec + size_t(exp0), x, rnd);
		if (is_negative) MPFR_NS mpfr_neg(x, x, rnd);
		if (ch == nullptr)
			throw std::runtime_error(
			    "in std::ostream& operator<<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
		std::string t(ch);
		MPFR_NS mpfr_free_str(ch);

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

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, std::string, _overwrite>
	{
		using real_type = real<_prec, _rnd>;
		using other_type = std::string;

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
		inline static int set(mpfr_ptr rop, const std::string& op, mpfr_rnd_t rnd)
		{
			auto res = helper_set_stdstr(rop, op, rnd);
			if (res == -1)
				throw std::runtime_error(
				    "in mpfr::helper_set_stdstr(mpfr_ptr, const std::string&, mpfr_rnd_t):\n  invalid input format " +
				    op);
			return res;
		}
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, char*, _overwrite>
	{
		using real_type = real<_prec, _rnd>;
		using other_type = char*;

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
		inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) { return helper_set_charptr(rop, op, rnd); }
	};

	template <real_prec_t _prec, real_rnd_t _rnd, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, const char*, _overwrite>
	{
		using real_type = real<_prec, _rnd>;
		using other_type = const char*;

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
		inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) { return helper_set_charptr(rop, op, rnd); }
	};

	template <real_prec_t _prec, real_rnd_t _rnd, size_t _i, bool _overwrite>
	struct type_traits<real<_prec, _rnd>, char[_i], _overwrite>
	{
		using real_type = real<_prec, _rnd>;
		using other_type = char[_i];

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
		inline static int set(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd) { return helper_set_charptr(rop, op, rnd); }
	};

	// //////////////////////////////////////////////////////////////
	// std::istream and std::ostream operators
	// //////////////////////////////////////////////////////////////

	template <real_prec_t _prec, real_rnd_t _rnd, class Char, class Traits>
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
					if (auto err = helper_set_stdstr(r._x, num, _rnd); err == -1)
					{
						in.setstate(std::ios_base::failbit);
						MPFR_NS mpfr_set_zero(r._x, +1);
						throw std::runtime_error(
						    "in std::istream& operator>>(std::istream& s, real<_prec, _rnd>& r):\n  invalid input "
						    "format " +
						    num);
					}
				}
				else
				{
					in.setstate(std::ios_base::failbit);
					MPFR_NS mpfr_set_zero(r._x, +1);
				}
			}
		}
		catch (std::ios_base::failure&)
		{
			MPFR_NS mpfr_set_zero(r._x, +1);
			throw;
		}
		return in;
	}

	template <real_prec_t _prec, real_rnd_t _rnd, class Char, class Traits>
	std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const real<_prec, _rnd>& r)
	{
		return helper_ostream(s, const_cast<mpfr_t&>(r._x), _rnd);
	}
}  // namespace mpfr

#endif  // REAL_IO_HPP_
