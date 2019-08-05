#include "matrix.hpp"

namespace algebra
{
	template class Vector<mpfr::real<1000, MPFR_RNDN>>;
	template class Matrix<mpfr::real<1000, MPFR_RNDN>>;
	template class Tensor<mpfr::real<1000, MPFR_RNDN>>;
}  // namespace algebra
