#include "chol_and_inverse.hpp"
#include "complex_function.hpp"
#include "context_variables.hpp"
#include "damped_rational.hpp"
#include "hor_formula.hpp"
#include "hor_recursion.hpp"
#include "integral_decomp.hpp"
#include "matrix.hpp"
#include "partial_fraction.hpp"
#include "pole_data.hpp"
#include "polynomial.hpp"
#include "primary_op.hpp"
#include "real.hpp"
#include "real_function.hpp"
#include "real_io.hpp"

using R = mpfr::real<1000, MPFR_RNDN>;
using algebra::Polynomial, algebra::Vector, algebra::Matrix;
template <class R>
using P = Polynomial<R>;
template <class R>
using V = Vector<R>;

namespace qboot
{
	template class MobiusTransformation<mpfr::real<1000, MPFR_RNDN>>;
	template class ExpLike<mpfr::real<1000, MPFR_RNDN>>;
	template class DampedRational<mpfr::real<1000, MPFR_RNDN>>;
}  // namespace qboot

[[maybe_unused]] static void test()
{
	// (p * x + q) / (r * x + s)
	// q == r == 0 && p == s => Id
	// r == 0 && p == s => S(q / s)
	// p == s == 0 && q == r => F
	// q == r == 0 => M(p / s)
	// p == 0 && q == r => F.S(s / r)
	// s == 0 && q == r => S(p / r).F
	// r == 0 => S(q / s).M(p / s)
	// p == s == 0 => M(q / r).F
	// q == 0 && p == s => F.S(r / p).F
	// p == 0 => F.S(s / q).M(r / q)
	// (q - r) * r == p * s => S(p / r).F.S(s / r) (special version of order 4, M(...) = M(1) = Id)
	// s == 0 => S(p / r).M(q / r).F
	// otherwise => S(p / r).M(-(p * s - q * r) / (r * r)).F.S(s / r)
	// [order 0]
	// @ Id = x => x
	// [order 1]
	// @ S(a) = x => x + a
	// @ F = x => 1 / x
	// @ M(a) = x => a * x
	// [order 2]
	//   F.F = Id
	// @ F.S(a) = x => 1 / (x + a)
	//   F.M(a) = M(1 / a).f
	//   S(a).S(b) = S(a + b)
	// @ S(a).F = x => (a * x + 1) / x
	// @ S(a).M(b) = x => b * x + a
	//   M(a).M(b) = M(a * b)
	//   M(a).S(b) = S(a * b).M(a)
	// @ M(a).F = x => a / x
	// [order 3]
	//   F.F.S(a) = S(a)
	// @ F.S(a).F = x => x / (a * x + 1)
	// @ F.S(a).M(b) = x => 1 / (b * x + a)
	//   F.M(a).F = M(1 / a)
	// @ S(a).F.S(b) = x => (a * x + a * b + 1) / (x + b)
	//   S(a).S(b).F = S(a + b).F
	//   S(a).S(b).M(c) = S(a + b).M(c)
	// @ S(a).M(b).F = x => (a * x + b) / x
	//   M(a).F.S(b) = F.S(b / a).M(1 / a)
	//   M(a).S(b).F = S(a * b).M(a).F
	//   M(a).S(b).M(c) = S(a * b).M(a * c)
	//   M(a).M(b).F = M(a * b).F
}
