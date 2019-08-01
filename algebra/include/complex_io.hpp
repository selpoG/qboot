#ifndef COMPLEX_IO_HPP_
#define COMPLEX_IO_HPP_

#include <istream>  // for basic_istream
#include <ostream>  // for basic_ostream
#include <sstream>  // for basic_ostringstream

#include "complex.hpp"
#include "real_io.hpp"

namespace algebra
{
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd, class CharT, class Traits>
	inline std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os,
	                                                     const complex<_prec, _rnd>& z)
	{
		std::basic_ostringstream<CharT, Traits> s;
		s.flags(os.flags());
		s.imbue(os.getloc());
		s.precision(os.precision());
		s << '(' << z.real() << ',' << z.imag() << ')';
		return os << s.str();
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd, class CharT, class Traits>
	inline std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, complex<_prec, _rnd>& z)
	{
		CharT head = 0;
		mpfr::real<_prec, _rnd> real = 0, imag = 0;
		if (is >> head && head != '(')
		{
			is.putback(head);
			is >> real;
		}
		else if (is >> real >> head && head != ',')
		{
			if (head != ')')
			{
				is.putback(head);
				is.setstate(std::ios_base::failbit);
			}
		}
		else if (is >> imag >> head && head != ')')
		{
			is.putback(head);
			is.setstate(std::ios_base::failbit);
		}

		if (!is.fail()) z = complex(real, imag);
		return is;
	}
}  // namespace algebra

#endif  // COMPLEX_IO_HPP_
