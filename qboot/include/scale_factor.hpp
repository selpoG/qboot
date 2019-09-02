#ifndef QBOOT_SCALE_FACTOR_HPP_
#define QBOOT_SCALE_FACTOR_HPP_

#include <cstddef>  // for uint32_t

#include "matrix.hpp"      // for Vector
#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for real

namespace qboot
{
	// corresponds to \chi_j(x) discussed in README.md,
	// which is a function [0, \infty) -> (0, \infty)
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
		// the max degree D among M[n](x) / chi(x)
		[[nodiscard]] virtual uint32_t max_degree() const = 0;
		// evaluate chi(x) at x = v
		// v must be non-negative
		[[nodiscard]] virtual Real eval(const Real& v) const = 0;
		// gives f(x_0), ..., f(x_D)
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<Real> sample_scalings() = 0;
		// gives x_k (0 <= k <= D)
		// this function might store cache data
		[[nodiscard]] virtual Real sample_point(uint32_t k) = 0;
		// gives x_0, ..., x_D
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<Real> sample_points() = 0;
		// gives q_m(x) (0 <= m <= D / 2)
		// this function might store cache data
		[[nodiscard]] virtual algebra::Polynomial<Real> bilinear_base(uint32_t m) = 0;
		// gives q_0(x), ..., q_{D / 2}
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<algebra::Polynomial<Real>> bilinear_bases() = 0;
	};
}  // namespace qboot

#endif  // QBOOT_SCALE_FACTOR_HPP_
