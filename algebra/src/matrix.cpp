#include "matrix.hpp"
#include "polynomial.hpp"

namespace algebra
{
	template class Vector<mpfr::real<1000, MPFR_RNDN>>;
	template class Matrix<mpfr::real<1000, MPFR_RNDN>>;
	template class Tensor<mpfr::real<1000, MPFR_RNDN>>;
	template class Polynomial<mpfr::real<1000, MPFR_RNDN>>;
	template class Vector<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Matrix<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Tensor<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Polynomial<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Vector<Polynomial<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Matrix<Polynomial<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Tensor<Polynomial<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Polynomial<Polynomial<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
}  // namespace algebra
