#ifndef CHOL_AND_INVERSE_H
#define CHOL_AND_INVERSE_H

#include <vector>

#include "matrix.hpp"
#include "real.hpp"

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Matrix<Real> anti_band_to_inverse(const std::vector<Real>& ab)
	{
		auto dim = (1 + ab.size()) / 2;
		algebra::Matrix<Real> A(dim, dim);
		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j < dim; j++) A.get(i, j) = ab[i + j];
		return A.cholesky_decomposition().lower_triangular_inverse();
	}
}  // namespace qboot

#endif  // CHOL_AND_INVERSE_H
