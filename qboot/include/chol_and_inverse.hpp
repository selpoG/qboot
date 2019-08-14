#ifndef CHOL_AND_INVERSE_HPP_
#define CHOL_AND_INVERSE_HPP_

#include <vector>  // for vector

#include "matrix.hpp"  // for Matrix
#include "real.hpp"    // for real

namespace qboot
{
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	algebra::Matrix<Ring> anti_band_to_inverse(const std::vector<Ring>& ab)
	{
		auto dim = (1 + ab.size()) / 2;
		algebra::Matrix<Ring> A(dim, dim);
		for (size_t i = 0; i < dim; ++i)
			for (size_t j = 0; j < dim; ++j) A.get(i, j) = ab[i + j];
		return A.cholesky_decomposition().lower_triangular_inverse();
	}
}  // namespace qboot

#endif  // CHOL_AND_INVERSE_HPP_
