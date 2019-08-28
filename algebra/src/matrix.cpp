#include "matrix.hpp"
#include "complex_function.hpp"
#include "polynomial.hpp"
#include "real_function.hpp"

using R = mpfr::real<1000, MPFR_RNDN>;

namespace algebra
{
	template class Vector<R>;
	template class Matrix<R>;
	template class Vector<Vector<R>>;
	template class Vector<Matrix<R>>;
	template class Polynomial<R>;
	template class RealFunction<R>;
	template class RealFunctionWithPower<R>;
	template class RealConverter<R>;
	template class ComplexFunction<R>;
	template class Vector<Polynomial<R>>;
	template class Matrix<Polynomial<R>>;
	template class Vector<RealFunction<R>>;
	template class Matrix<RealFunction<R>>;
	template class Vector<ComplexFunction<R>>;
	template class Matrix<ComplexFunction<R>>;
	template class Vector<RealFunction<Polynomial<R>>>;
	template class Vector<ComplexFunction<Polynomial<R>>>;
	template class Matrix<RealFunction<Polynomial<R>>>;
	template class Matrix<ComplexFunction<Polynomial<R>>>;
}  // namespace algebra
