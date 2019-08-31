#ifndef QBOOT_SCALE_FACTOR_HPP_
#define QBOOT_SCALE_FACTOR_HPP_

#include <cstddef>  // for uint32_t

#include "matrix.hpp"      // for Vector
#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for real

namespace qboot
{
	// a common scale factor f(x) of elements in one constraint of a polynomial matrix programming
	// in a semidefiniteness
	//   \sum_{n = 0}^{N - 1} y[n] M[n](x) >= C(x) for all x >= 0,
	// M[n](x) for all n and C(x) must be a polynomial of x,
	// but in conformal bootstrap, they have a common factor f(x) like exp(-k x) / (polynomial of x),
	// which is positive in x >= 0.
	// we can divide the inequality by this factor, but natural scalings might be changed.
	// So, we allow M[n](x) and C(x) to be general functions of x,
	// but they must become polynomials when divided by a scale factor f(x).
	template <class Real>
	class ScaleFactor
	{
	public:
		ScaleFactor() = default;
		ScaleFactor(const ScaleFactor&) = default;
		ScaleFactor& operator=(const ScaleFactor&) = default;
		ScaleFactor(ScaleFactor&&) noexcept = default;
		ScaleFactor& operator=(ScaleFactor&&) noexcept = default;
		virtual ~ScaleFactor() = default;
		// the max degree D among M[n](x) / f(x), C(x) / f(x)
		[[nodiscard]] virtual uint32_t max_degree() const = 0;
		// evaluate f(x) at x = v
		// v must be non-negative
		[[nodiscard]] virtual Real eval(const Real& v) const = 0;
		// gives f(x_0), ..., f(x_{D})
		[[nodiscard]] virtual algebra::Vector<Real> sample_scalings() const& = 0;
		// gives f(x_0), ..., f(x_{D})
		[[nodiscard]] virtual algebra::Vector<Real> sample_scalings() && = 0;
		// gives x_k (0 <= k <= D)
		[[nodiscard]] virtual Real sample_point(uint32_t k) const = 0;
		// gives x_0, ..., x_{D}
		[[nodiscard]] virtual algebra::Vector<Real> sample_points() const& = 0;
		// gives x_0, ..., x_{D}
		[[nodiscard]] virtual algebra::Vector<Real> sample_points() && = 0;
		// gives q_m(x) (0 <= m <= D / 2)
		[[nodiscard]] virtual algebra::Polynomial<Real> bilinear_base(uint32_t m) const = 0;
		// gives q_0(x), ..., q_{D / 2}
		[[nodiscard]] virtual algebra::Vector<algebra::Polynomial<Real>> bilinear_bases() const& = 0;
		// gives q_0(x), ..., q_{D / 2}
		[[nodiscard]] virtual algebra::Vector<algebra::Polynomial<Real>> bilinear_bases() && = 0;
	};
}  // namespace qboot

#endif  // QBOOT_SCALE_FACTOR_HPP_
