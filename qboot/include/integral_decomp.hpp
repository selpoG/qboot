#ifndef INTEGRAL_DECOMP_HPP_
#define INTEGRAL_DECOMP_HPP_

#include <cstdint>  // for uint32_t

#include "matrix.hpp"  // for Vector
#include "real.hpp"    // for real

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> simple_pole_integral(uint32_t pole_order_max, const Real& base, const Real& pole_position,
	                                           const Real& incomplete_gamma_factor);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> double_pole_integral(uint32_t pole_order_max, const Real& base, const Real& pole_position,
	                                           const Real& incomplete_gamma_factor);
	// let pole_position = -p, base = exp(-k)
	// calculate \int_{0}^{\infty} e ^ {-k x} x ^ n / (x + p) dx, for n = 0, ..., pole_order_max
	// this integral equals to
	// n! p ^ n e ^ {p k} \Gamma(-n, p k)
	// = (-p) ^ n e ^ {p k} \Gamma(0, p k)
	//   + (1 / k ^ n) \sum_{i = 0}^{n - 1} (n - i - 1)! (-p k) ^ i
	// incomplete_gamma_factor = e ^ {p k} \Gamma(0, p k)
	template <class Real>
	algebra::Vector<Real> simple_pole_integral(uint32_t pole_order_max, const Real& base, const Real& pole_position,
	                                           const Real& incomplete_gamma_factor)
	{
		algebra::Vector<Real> result(pole_order_max + 1);
		result[0] = incomplete_gamma_factor;
		Real tmp{}, pow = pole_position * incomplete_gamma_factor;
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

	// let pole_position = -p, base = \exp(-k)
	// calculate \int_{0}^{\infty} e ^ {-k x} x ^ n / (x + p) ^ 2 dx, for n = 0, ..., pole_order_max
	template <class Real>
	algebra::Vector<Real> double_pole_integral(uint32_t pole_order_max, const Real& base, const Real& pole_position,
	                                           const Real& incomplete_gamma_factor)
	{
		algebra::Vector<Real> result(pole_order_max + 1);

		Real tmp = mpfr::log(base);
		Real minus_pole_position = -pole_position;
		Real minus_log_base = -1 / mpfr::log(base);
		Real log_base_power = minus_log_base;

		result[0] = incomplete_gamma_factor * tmp;
		tmp = 1 / minus_pole_position;
		result[0] += tmp;
		// result[0] == 1 / minus_pole_position + incomplete_gamma_factor * log(base);

		for (uint32_t i = 1; i <= pole_order_max; ++i) result[i] = result[i - 1] * pole_position;

		algebra::Vector<Real> factorial_times_power_lnb(0);
		algebra::Vector<Real> single_pole_coeffs(0);

		// x / (x + a) ^ 2 case

		if (pole_order_max >= 1)
		{
			single_pole_coeffs = algebra::Vector<Real>(pole_order_max);
			// x ^ (j + 1) = (
			// single_pole_coeffs[0] x ^ (j - 1) +
			// single_pole_coeffs[1] x ^ (j - 2) + ... +
			// single_pole_coeffs[j - 1] x^0
			// ) * (x - a) ^ 2 +
			//
			// single_pole_coeffs[j](x - a) * ((x - a) + a)
			//
			// + a ^ (j + 1)
			//
			// => single_pole_coeffs[j + 1]
			//
			// single_pole_coeffs[j + 1] = single_pole_coeffs[j] * a + a ^ j + 1
			// single_pole_coeffs[0] =
			if (pole_order_max >= 2)
			{
				factorial_times_power_lnb = algebra::Vector<Real>(pole_order_max - 1);
				factorial_times_power_lnb[0] = minus_log_base;
			}

			tmp = pole_position;
			// below tmp is used as pole_position ^ j

			single_pole_coeffs[0] = 1;
			result[1] += incomplete_gamma_factor;

			for (uint32_t j = 1; j + 1 <= pole_order_max; ++j)
			{
				single_pole_coeffs[j] = single_pole_coeffs[j - 1] * pole_position + tmp;

				result[j + 1] += single_pole_coeffs[j] * incomplete_gamma_factor;
				if (j + 2 <= pole_order_max)
				{
					tmp *= pole_position;
					factorial_times_power_lnb[j] = factorial_times_power_lnb[j - 1] * minus_log_base;
					factorial_times_power_lnb[j] *= j;
				}
			}
		}

		for (uint32_t j = 0; j + 2 <= pole_order_max; ++j)
			for (uint32_t k = 0; k <= j; ++k) result[j + 2] += factorial_times_power_lnb[k] * single_pole_coeffs[j - k];

		return result;
	}
}  // namespace qboot

#endif  // INTEGRAL_DECOMP_HPP_
