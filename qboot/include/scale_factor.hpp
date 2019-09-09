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
		[[nodiscard]] virtual mpfr::real eval(const mpfr::real& v) const = 0;
		// gives f(x_0), ..., f(x_D)
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<mpfr::real> sample_scalings() = 0;
		// gives x_k (0 <= k <= D)
		// this function might store cache data
		[[nodiscard]] virtual mpfr::real sample_point(uint32_t k) = 0;
		// gives x_0, ..., x_D
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<mpfr::real> sample_points() = 0;
		// gives q_m(x) (0 <= m <= D / 2)
		// this function might store cache data
		[[nodiscard]] virtual algebra::Polynomial bilinear_base(uint32_t m) = 0;
		// gives q_0(x), ..., q_{D / 2}
		// this function might store cache data
		[[nodiscard]] virtual algebra::Vector<algebra::Polynomial> bilinear_bases() = 0;
	};
	class TrivialScale : public ScaleFactor
	{
	public:
		TrivialScale() = default;
		TrivialScale(const TrivialScale&) = delete;
		TrivialScale& operator=(const TrivialScale&) = delete;
		TrivialScale(TrivialScale&&) noexcept = default;
		TrivialScale& operator=(TrivialScale&&) noexcept = default;
		~TrivialScale() override = default;
		[[nodiscard]] uint32_t max_degree() const override { return 0; }
		[[nodiscard]] mpfr::real eval([[maybe_unused]] const mpfr::real& v) const override { return mpfr::real(1); }
		[[nodiscard]] algebra::Vector<mpfr::real> sample_scalings() override
		{
			return algebra::Vector<mpfr::real>{mpfr::real(1)};
		}
		[[nodiscard]] mpfr::real sample_point([[maybe_unused]] uint32_t k) override { return mpfr::real(0); }
		[[nodiscard]] algebra::Vector<mpfr::real> sample_points() override
		{
			return algebra::Vector<mpfr::real>{mpfr::real(1)};
		}
		[[nodiscard]] algebra::Polynomial bilinear_base([[maybe_unused]] uint32_t m) override
		{
			return {mpfr::real(1)};
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial> bilinear_bases() override
		{
			return algebra::Vector<algebra::Polynomial>{algebra::Polynomial(mpfr::real(1))};
		}
	};
}  // namespace qboot

#endif  // QBOOT_SCALE_FACTOR_HPP_
