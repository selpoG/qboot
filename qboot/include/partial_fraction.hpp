#ifndef PARTIAL_FRACTION_H
#define PARTIAL_FRACTION_H

#include <vector>

#include "matrix.hpp"
#include "real.hpp"

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> fast_partial_fraction(const algebra::Vector<Real>& pole_locations,
	                                              const std::vector<bool>& is_double, size_t n_poles)
	{
		size_t expected_result_length = n_poles;
		for (size_t i = 0; i < n_poles; i++)
			if (is_double[i]) ++expected_result_length;
		algebra::Vector<Real> result(expected_result_length);
		Real tmp;
		for (size_t i = 0, pos = 0; i < n_poles; i++, pos++)
		{
			result[pos] = 1;
			// product (pole[i] - pole[j]) ^ (1 or 2), j != i
			for (size_t j = 0; j < n_poles; j++)
			{
				if (i == j) continue;
				tmp = pole_locations[i] - pole_locations[j];
				result[pos] *= tmp;
				if (is_double[j]) result[pos] *= tmp;
			}
			result[pos] = 1 / result[pos];
			if (is_double[i])
			{
				++pos;
				for (size_t j = 0; j < n_poles; j++)
				{
					if (i == j) continue;
					result[pos] += (is_double[j] ? 2 : 1) / (pole_locations[i] - pole_locations[j]);
				}
				result[pos] *= -result[pos - 1];
			}
		}
		return result;
	}
}  // namespace qboot

#endif  // PARTIAL_FRACTION_H
