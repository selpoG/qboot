#include "matrix.hpp"
#include "complex_function.hpp"
#include "polynomial.hpp"
#include "real_function.hpp"

namespace algebra
{
	template class Vector<mpfr::real<1000, MPFR_RNDN>>;
	template class Matrix<mpfr::real<1000, MPFR_RNDN>>;
	template class Vector<Vector<mpfr::real<1000, MPFR_RNDN>>>;
	template class Vector<Matrix<mpfr::real<1000, MPFR_RNDN>>>;
	template class Polynomial<mpfr::real<1000, MPFR_RNDN>>;
	template class RealFunction<mpfr::real<1000, MPFR_RNDN>>;
	template class RealFunctionWithPower<mpfr::real<1000, MPFR_RNDN>>;
	template class RealConverter<mpfr::real<1000, MPFR_RNDN>>;
	template class ComplexFunction<mpfr::real<1000, MPFR_RNDN>>;
	template class Vector<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Matrix<Polynomial<mpfr::real<1000, MPFR_RNDN>>>;
	template class Vector<RealFunction<mpfr::real<1000, MPFR_RNDN>>>;
	template class Matrix<RealFunction<mpfr::real<1000, MPFR_RNDN>>>;
	template class Vector<ComplexFunction<mpfr::real<1000, MPFR_RNDN>>>;
	template class Matrix<ComplexFunction<mpfr::real<1000, MPFR_RNDN>>>;
	template class Vector<RealFunction<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Vector<ComplexFunction<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Matrix<RealFunction<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
	template class Matrix<ComplexFunction<Polynomial<mpfr::real<1000, MPFR_RNDN>>>>;
}  // namespace algebra
