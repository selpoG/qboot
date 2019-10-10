#ifndef QBOOT_CONFORMAL_SCALE_HPP_
#define QBOOT_CONFORMAL_SCALE_HPP_

#include <cassert>   // for assert
#include <cstdint>   // for uint32_t, int32_t
#include <optional>  // for optional
#include <utility>   // for move

#include "qboot/algebra/matrix.hpp"      // for Vector
#include "qboot/algebra/polynomial.hpp"  // for Polynomial
#include "qboot/block.hpp"               // for ConformalBlock
#include "qboot/context.hpp"             // for Context
#include "qboot/mp/real.hpp"             // for real, pow, log
#include "qboot/scale_factor.hpp"        // for ScaleFactor

namespace qboot
{
	inline bool include_odd(const mp::real& d1, const mp::real& d2, const mp::real& d3, const mp::real& d4)
	{
		return d1 != d2 && d3 != d4;
	}
	inline bool include_odd(const mp::real& d12, const mp::real& d34) { return d12 != 0 && d34 != 0; }
	class ConformalScale : public ScaleFactor
	{
		bool odd_included_;
		uint32_t spin_, lambda_;
		mp::rational epsilon_;
		mp::real rho_, gap_;
		std::optional<mp::real> end_;
		// pole at Delta = poles_[i]
		// 1 / (Delta - poles_[i])
		algebra::Vector<mp::real> poles_;
		algebra::Vector<algebra::Polynomial> bilinear_bases_{};
		void _set_bilinear_bases() &;

	public:
		void _reset() &&
		{
			odd_included_ = false;
			spin_ = lambda_ = 0;
			std::move(epsilon_)._reset();
			std::move(rho_)._reset();
			std::move(gap_)._reset();
			end_.reset();
			std::move(poles_)._reset();
			std::move(bilinear_bases_)._reset();
		}
		ConformalScale(const GeneralPrimaryOperator& op, const Context& context, bool include_odd)
		    : ConformalScale(op.num_poles(), op.spin(), context, include_odd, op.lower_bound(), op.upper_bound_safe())
		{
		}
		// constrain delta to be in the range gap <= delta < end
		// if end is nullopt (corresponds to infinity), gap <= delta
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, bool include_odd, const mp::real& gap,
		               const std::optional<mp::real>& end);
		ConformalScale(const ConformalScale&) = delete;
		ConformalScale& operator=(const ConformalScale&) = delete;
		ConformalScale(ConformalScale&&) noexcept = default;
		ConformalScale& operator=(ConformalScale&&) noexcept = default;
		~ConformalScale() override;
		[[nodiscard]] bool odd_included() const { return odd_included_; }
		[[nodiscard]] uint32_t max_degree() const override { return poles_.size() + lambda_; }
		[[nodiscard]] const algebra::Vector<mp::real>& get_poles() const& { return poles_; }
		[[nodiscard]] algebra::Vector<mp::real> get_poles() && { return std::move(poles_); }
		// evaluate at delta
		[[nodiscard]] mp::real eval_d(const mp::real& delta) const
		{
			mp::real ans = mp::pow(4 * rho_, delta);
			if (end_.has_value()) ans *= mp::pow(get_x(delta) + end_.value(), -int32_t(max_degree()));
			for (const auto& p : poles_) ans /= delta - p;
			return ans;
		}
		// convert x (in [0, \infty)) to delta
		[[nodiscard]] mp::real get_delta(const mp::real& x) const
		{
			if (end_.has_value()) return end_.value() * (x + gap_) / (x + end_.value());
			return x + gap_;
		}
		// convert delta to x (in [0, \infty))
		[[nodiscard]] mp::real get_x(const mp::real& delta) const
		{
			if (end_.has_value()) return end_.value() * (delta - gap_) / (end_.value() - gap_);
			return delta - gap_;
		}
		// evaluate at delta
		[[nodiscard]] mp::real eval(const mp::real& x) const override { return eval_d(get_delta(x)); }
		[[nodiscard]] algebra::Vector<mp::real> sample_scalings() const override
		{
			auto xs = sample_points();
			for (auto& x : xs) x = eval(x);
			return xs;
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial> bilinear_bases() const override
		{
			return bilinear_bases_.clone();
		}
		[[nodiscard]] mp::real sample_point(uint32_t k) const override
		{
			return ConformalScale::sample_point(max_degree(), k);
		}
		[[nodiscard]] algebra::Vector<mp::real> sample_points() const override
		{
			return ConformalScale::sample_points(max_degree());
		}
		static mp::real sample_point(uint32_t degree, uint32_t k)
		{
			return mp::pow(mp::const_pi() * (-1 + 4 * int32_t(k)), 2uL) /
			       ((-64 * int32_t(degree + 1)) * mp::log(3 - mp::sqrt(8)));
		}
		static algebra::Vector<mp::real> sample_points(uint32_t degree)
		{
			algebra::Vector<mp::real> v(degree + 1);
			for (uint32_t i = 0; i <= degree; ++i) v[i] = mp::pow(mp::const_pi() * (-1 + 4 * int32_t(i)), 2uL);
			v /= (-64 * int32_t(degree + 1)) * mp::log(3 - mp::sqrt(8));
			return v;
		}
	};
}  // namespace qboot

#endif  // QBOOT_CONFORMAL_SCALE_HPP_
