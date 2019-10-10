#ifndef QBOOT_SCALE_FACTOR_HPP_
#define QBOOT_SCALE_FACTOR_HPP_

#include <cstddef>  // for uint32_t

#include "qboot/algebra/matrix.hpp"      // for Vector
#include "qboot/algebra/polynomial.hpp"  // for Polynomial
#include "qboot/mp/real.hpp"             // for real

namespace qboot
{
	// corresponds to \chi_j(x) discussed in README.md,
	// which is a function [0, \infty) -> (0, \infty)
	class ScaleFactor
	{
	public:
		ScaleFactor() = default;
		ScaleFactor(const ScaleFactor&) = default;
		ScaleFactor& operator=(const ScaleFactor&) = default;
		ScaleFactor(ScaleFactor&&) noexcept = default;
		ScaleFactor& operator=(ScaleFactor&&) noexcept = default;
		virtual ~ScaleFactor();
		// the max degree D among M[n](x) / chi(x)
		[[nodiscard]] virtual uint32_t max_degree() const = 0;
		// evaluate chi(x) at x = v
		// v must be non-negative
		[[nodiscard]] virtual mp::real eval(const mp::real& v) const = 0;
		// gives f(x_0), ..., f(x_D)
		[[nodiscard]] virtual algebra::Vector<mp::real> sample_scalings() const = 0;
		// gives x_k (0 <= k <= D)
		[[nodiscard]] virtual mp::real sample_point(uint32_t k) const = 0;
		// gives x_0, ..., x_D
		[[nodiscard]] virtual algebra::Vector<mp::real> sample_points() const = 0;
		// gives q_0(x), ..., q_{D / 2}
		[[nodiscard]] virtual algebra::Vector<algebra::Polynomial> bilinear_bases() const = 0;
	};
	class TrivialScale : public ScaleFactor
	{
	public:
		TrivialScale() = default;
		TrivialScale(const TrivialScale&) = delete;
		TrivialScale& operator=(const TrivialScale&) = delete;
		TrivialScale(TrivialScale&&) noexcept = default;
		TrivialScale& operator=(TrivialScale&&) noexcept = default;
		~TrivialScale() override;
		[[nodiscard]] uint32_t max_degree() const override { return 0; }
		[[nodiscard]] mp::real eval([[maybe_unused]] const mp::real& v) const override { return mp::real(1); }
		[[nodiscard]] algebra::Vector<mp::real> sample_scalings() const override
		{
			return algebra::Vector<mp::real>{mp::real(1)};
		}
		[[nodiscard]] mp::real sample_point([[maybe_unused]] uint32_t k) const override { return mp::real(0); }
		[[nodiscard]] algebra::Vector<mp::real> sample_points() const override
		{
			return algebra::Vector<mp::real>{mp::real(1)};
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial> bilinear_bases() const override
		{
			return algebra::Vector<algebra::Polynomial>{algebra::Polynomial(mp::real(1))};
		}
	};
}  // namespace qboot

#endif  // QBOOT_SCALE_FACTOR_HPP_
