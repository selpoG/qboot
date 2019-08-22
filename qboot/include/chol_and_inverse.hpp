#ifndef CHOL_AND_INVERSE_HPP_
#define CHOL_AND_INVERSE_HPP_

#include <cstdint>  // for uint32_t

#include "matrix.hpp"  // for Matrix, Vector
#include "real.hpp"    // for real

namespace qboot
{
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	algebra::Matrix<Ring> anti_band_to_inverse(const algebra::Vector<Ring>& ab)
	{
		auto dim = (1 + ab.size()) / 2;
		algebra::Matrix<Ring> A(dim, dim);
		for (uint32_t i = 0; i < dim; ++i)
			for (uint32_t j = 0; j < dim; ++j) A.at(i, j) = ab[i + j];
		return A.cholesky_decomposition().lower_triangular_inverse();
	}
}  // namespace qboot

#endif  // CHOL_AND_INVERSE_HPP_
