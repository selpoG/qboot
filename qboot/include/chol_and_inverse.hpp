#ifndef QBOOT_CHOL_AND_INVERSE_HPP_
#define QBOOT_CHOL_AND_INVERSE_HPP_

#include <cstdint>  // for uint32_t

#include "matrix.hpp"  // for Matrix, Vector
#include "real.hpp"    // for real

namespace qboot
{
	inline algebra::Matrix<mpfr::real> anti_band_to_inverse(const algebra::Vector<mpfr::real>& ab)
	{
		auto dim = (1 + ab.size()) / 2;
		algebra::Matrix<mpfr::real> A(dim, dim);
		for (uint32_t i = 0; i < dim; ++i)
			for (uint32_t j = 0; j < dim; ++j) A.at(i, j) = ab[i + j];
		return std::move(A).cholesky_decomposition().lower_triangular_inverse();
	}
}  // namespace qboot

#endif  // QBOOT_CHOL_AND_INVERSE_HPP_
