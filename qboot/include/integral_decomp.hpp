#ifndef QBOOT_INTEGRAL_DECOMP_HPP_
#define QBOOT_INTEGRAL_DECOMP_HPP_

#include <cstdint>  // for uint32_t

#include "matrix.hpp"  // for Vector
#include "real.hpp"    // for real, gamma_inc, pow, log

namespace qboot
{
	// let pole_position = -p, base = exp(-k)
	// calculate \int_{0}^{\infty} e ^ {-k x} x ^ n / (x + p) dx, for n = 0, ..., pole_order_max
	// this integral equals to
	// n! p ^ n e ^ {p k} \Gamma(-n, p k)
	// = (-p) ^ n e ^ {p k} \Gamma(0, p k)
	//   + (1 / k ^ n) \sum_{i = 0}^{n - 1} (n - i - 1)! (-p k) ^ i
	// incomplete_gamma_factor = e ^ {p k} \Gamma(0, p k)
	template <class Real>
	algebra::Vector<Real> simple_pole_integral(uint32_t pole_order_max, const Real& base, const Real& pole_position)
	{
		Real incomplete_gamma = pole_position == 0 ? Real(Real::prec)
		                                           : mpfr::pow(base, pole_position) *
		                                                 mpfr::gamma_inc(Real(0), pole_position * mpfr::log(base));
		algebra::Vector<Real> result(pole_order_max + 1);
		result[0] = incomplete_gamma;
		Real tmp{}, pow = pole_position * incomplete_gamma;
		Real minus_pole_position = -pole_position;
		Real factorial(1);
		Real minus_log_base = -1 / mpfr::log(base);
		Real log_base_power = minus_log_base;
		for (uint32_t j = 1; j <= pole_order_max; ++j)
		{
			// factorial == (j - 1)!;
			// pow == incomplete_gamma_factor * pow(pole_position, j);
			// log_base_power == pow(minus_log_base, j);
			tmp = factorial * log_base_power + tmp * pole_position;
			// tmp == sum((j - k - 1)! * pow(pole_position, k) * pow(minus_log_base, j - k), 0 <= k < j);
			result[j] = tmp + pow;
			// result[j] == sum((j - k - 1)! * pow(pole_position, k) * pow(minus_log_base, j - k), 0 <= k < j)
			//              + incomplete_gamma_factor * pow(pole_position, j);

			if (j < pole_order_max)
			{
				pow *= pole_position;
				log_base_power *= minus_log_base;
				factorial *= j;
			}
		}
		return result;
	}
}  // namespace qboot

#endif  // QBOOT_INTEGRAL_DECOMP_HPP_
