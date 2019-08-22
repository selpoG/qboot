#ifndef PARTIAL_FRACTION_HPP_
#define PARTIAL_FRACTION_HPP_

#include <cstdint>  // for uint32_t

#include "matrix.hpp"  // for Vector
#include "real.hpp"    // for real

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> fast_partial_fraction(const algebra::Vector<Real>& poles)
	{
		algebra::Vector<Real> result(poles.size());
		for (uint32_t i = 0; i < poles.size(); ++i)
		{
			result[i] = 1;
			for (uint32_t j = 0; j < poles.size(); ++j)
				if (i != j) result[i] *= poles[i] - poles[j];
			result[i] = 1 / result[i];
		}
		return result;
	}
}  // namespace qboot

#endif  // PARTIAL_FRACTION_HPP_
