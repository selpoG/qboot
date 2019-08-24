#ifndef REAL_IO_HPP_
#define REAL_IO_HPP_

#include "real.hpp"

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

namespace mpfr
{
	// exception type

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wweak-vtables"
#endif
	class exception_real : public std::exception
	{
	public:
		explicit exception_real(const std::string& msg = "exception_real") noexcept : _msg(msg) {}
		~exception_real() noexcept override = default;
		exception_real& operator=(const exception_real&) = default;
		exception_real(const exception_real&) = default;
		exception_real& operator=(exception_real&&) noexcept = default;
		exception_real(exception_real&&) noexcept = default;
		// returns cause of error
		const char* what() const noexcept override { return _msg.c_str(); }

	private:
		std::string _msg;
	};
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	// //////////////////////////////////////////////////////////////////
	// helper functions
	// //////////////////////////////////////////////////////////////////

	inline int helper_set_stdstr(mpfr_ptr rop, const std::string& op, mpfr_rnd_t rnd)
	{
		const int err = MPFR_NS mpfr_set_str(rop, op.c_str(), 0, rnd);
		if (err == -1)
			throw exception_real(
			    std::string(
			        "in mpfr::helper_set_stdstr(mpfr_ptr, const std::string&, mpfr_rnd_t):\n  invalid input format ") +
			    op);
		return err;
	}

	inline int helper_set_charptr(mpfr_ptr rop, const char* op, mpfr_rnd_t rnd)
	{
		const int err = MPFR_NS mpfr_set_str(rop, op, 0, rnd);
		if (err == -1)
			throw exception_real(
			    std::string(
			        "in mpfr::helper_set_charptr(mpfr_ptr, const char*, mpfr_rnd_t):\n  invalid input format ") +
			    op);
		return err;
	}

	// there might be some room for improvements for the next
	// missing: handling of ios_base::fixed/ios_base::scientific and
	// ios_base::showpoint

	template <class Char, class Traits>
	inline std::basic_ostream<Char, Traits>& helper_ostream(std::basic_ostream<Char, Traits>& s, const mpfr_t x,
	                                                        mpfr_rnd_t rnd)
	{
		real_exp_t exp;
		auto ch = MPFR_NS mpfr_get_str(nullptr, &exp, 10, static_cast<size_t>(s.precision() + 1), x, rnd);
		if (ch == nullptr)
			throw exception_real(
			    "in std::ostream& operator <<(std::ostream& s, const real<_prec, _rnd>& r):\n  conversion failed");
		auto t = std::string(ch);
		MPFR_NS mpfr_free_str(ch);

		const std::ios_base::fmtflags flags = s.flags();
		auto t_iter = t.begin();

		if (*t_iter == '-') ++t_iter;

		// digit?
		if (*t_iter == '0' || *t_iter == '1' || *t_iter == '2' || *t_iter == '3' || *t_iter == '4' || *t_iter == '5' ||
		    *t_iter == '6' || *t_iter == '7' || *t_iter == '8' || *t_iter == '9')
		{
			// positive sign
			if ((t_iter == t.begin()) && (flags & std::ios_base::showpos) != 0)
			{
				t_iter = t.insert(t_iter, '+');
				++t_iter;
			}

			// decimal point
			++t_iter;
			t.insert(t_iter, '.');

			// fixing exponent after insertion of decimal point
			// why must life be so difficult? (any suggestions for improvements?)
			if (!MPFR_NS mpfr_zero_p(x))
			{
				const real_exp_t exp_prev = exp;
				volatile real_exp_t* exp_ptr = &exp;
				--exp;
				if (*exp_ptr > exp_prev)
					throw exception_real(
					    "in std::ostream& operator <<(std::ostream& s, const real<_prec, _rnd>& r):\n  exponent out of "
					    "range");
			}

			// composing of the exponent
			if ((flags & std::ios_base::uppercase) != 0)
				t += 'E';
			else
				t += 'e';
			if (exp >= 0)
				t += '+';
			else
			{
				t += '-';
				exp = -exp;
			}
			if (exp >= -9 && exp <= 9) t += '0';
			std::stringstream temp;
			temp << exp;
			t += temp.str();
		}

		// width and adjustment
		if (s.width() > 0 && static_cast<unsigned int>(s.width()) > t.size())
		{
			if ((flags & std::ios_base::left) != 0)
				t_iter = t.end();
			else if ((flags & std::ios_base::internal) != 0)
			{
				t_iter = t.begin();
				if (*t_iter == '+' || *t_iter == '-') ++t_iter;
			}
			else
				t_iter = t.begin();
			while (t.size() < static_cast<unsigned int>(s.width())) t_iter = t.insert(t_iter, s.fill());
		}

		s << t;

		return s;
	}
	enum
	{
		MANT_SIGN = 0x01u,    // leading sign
		MANT_DIGIT = 0x02u,   // digits before decimal point
		MANT_POINT = 0x04u,   // decimal point
		MANT_FDIGIT = 0x08u,  // digits after decimal point
		EXP_SYMBOL = 0x10u,   // symbol of exponent ('e' or 'E')
		EXP_SIGN = 0x20u,     // sign of exponent
		EXP_DIGIT = 0x40u,    // digits of exponent
		MASK_EXP = (EXP_SYMBOL | EXP_SIGN | EXP_DIGIT),
		MASK_NINT = (MANT_POINT | MANT_FDIGIT | MASK_EXP)  // non-integral
	};

	inline bool helper_extract_float(std::istream& s, std::string* num)
	{
		auto ok = true;
		unsigned parts = 0x00;

		char c;
		while (s.get(c) && ok)
		{
			if (c == '+' || c == '-')
			{
				// very beginning
				if (parts == 0x00)
				{
					*num += c;
					parts |= MANT_SIGN;
				}
				// has symbol of exponent, but not yet a sign or digit
				else if ((parts & MASK_EXP) == EXP_SYMBOL)
				{
					*num += c;
					parts |= EXP_SIGN;
				}
				// end of number
				else
				{
					s.putback(c);
					break;
				}
			}
			else if (c == '.')
			{
				// does not yet have a decimal point or anything after it
				if ((parts & MASK_NINT) == 0x00)
				{
					*num += c;
					parts |= MANT_POINT;
				}
				// end of number
				else
				{
					s.putback(c);
					break;
				}
			}
			else if (c == 'e' || c == 'E')
			{
				// must have a digit && must not yet have an expontential
				if ((parts & unsigned(MANT_DIGIT | MANT_FDIGIT)) != 0x00 && (parts & MASK_EXP) == 0x00)
				{
					*num += c;
					parts |= EXP_SYMBOL;
				}
				// bad syntax
				else
				{
					s.putback(c);
					ok = false;
				}
			}
			else if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' || c == '5' || c == '6' || c == '7' ||
			         c == '8' || c == '9')
			{
				// before decimal point
				if ((parts & MASK_NINT) == 0x00)
				{
					*num += c;
					parts |= MANT_DIGIT;
				}
				// after decimal point
				else if ((parts & MASK_EXP) == 0x00)
				{
					*num += c;
					parts |= MANT_FDIGIT;
				}
				// in exponent
				else if ((parts & EXP_SYMBOL) != 0x00)
				{
					*num += c;
					parts |= EXP_DIGIT;
				}
				// some strange error?
				else
				{
					s.putback(c);
					ok = false;
				}
			}
			// other character => end of parsing
			else
			{
				s.putback(c);
				break;
			}
		}  // while (s.good() && ok)

		// further syntax checks
		// must have a digit, if a character has been parsed
		if (parts != 0x00 && (parts & unsigned(MANT_DIGIT | MANT_FDIGIT)) == 0x00) ok = false;
		// must have a digit in exponent, if symbol of exponent is set
		else if ((parts & EXP_SYMBOL) != 0x00 && (parts & EXP_DIGIT) == 0x00)
			ok = false;

		return ok;
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
			return helper_set_stdstr(rop, op, rnd);
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
	std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& s, real<_prec, _rnd>& r)
	{
		std::istream::sentry cerberos(s, false);

		if (cerberos)
		{
			// extract number
			std::string num;
			bool ok = helper_extract_float(s, &num);

			// bad syntax
			if (!ok)
			{
				s.setstate(std::ios_base::failbit);
				MPFR_NS mpfr_set_zero(r._x, +1);
			}
			// not empty (which could be due to an EOF)
			else if (!num.empty())
			{
				// conversion to mpfr::real
				try
				{
					helper_set_stdstr(r._x, num, _rnd);
				}
				// should, in principle, never fail, but ...
				catch (...)
				{
					s.setstate(std::ios_base::failbit);
					MPFR_NS mpfr_set_zero(r._x, +1);
					throw exception_real(std::string("in std::istream& operator >>(std::istream& s, real<_prec, "
					                                 "_rnd>& r):\n  invalid input format ") +
					                     num);
				}
			}
		}

		return s;
	}

	template <real_prec_t _prec, real_rnd_t _rnd, class Char, class Traits>
	std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const real<_prec, _rnd>& r)
	{
		return helper_ostream(s, r._x, _rnd);
	}
}  // namespace mpfr

#endif  // REAL_IO_HPP_
