#ifndef QBOOT_CONFORMAL_SCALE_HPP_
#define QBOOT_CONFORMAL_SCALE_HPP_

#include <cassert>   // for assert
#include <cstdint>   // for uint32_t, int32_t
#include <memory>    // for make_unique, unique_ptr
#include <optional>  // for optional
#include <utility>   // for move

#include "block.hpp"         // for ConformalBlock
#include "context.hpp"       // for Context
#include "matrix.hpp"        // for Vector
#include "polynomial.hpp"    // for Polynomial
#include "real.hpp"          // for real, pow, log
#include "scale_factor.hpp"  // for ScaleFactor

namespace qboot
{
	inline mp::real sample_point(uint32_t degree, uint32_t k)
	{
		return mp::pow(mp::const_pi() * (-1 + 4 * int32_t(k)), 2uL) /
		       ((-64 * int32_t(degree + 1)) * mp::log(3 - mp::sqrt(8)));
	}
	inline algebra::Vector<mp::real> sample_points(uint32_t degree)
	{
		algebra::Vector<mp::real> v(degree + 1);
		for (uint32_t i = 0; i <= degree; ++i) v[i] = mp::pow(mp::const_pi() * (-1 + 4 * int32_t(i)), 2uL);
		v /= (-64 * int32_t(degree + 1)) * mp::log(3 - mp::sqrt(8));
		return v;
	}
	class PoleSequence
	{
		bool include_odd;
		uint32_t type, k, spin;
		mp::rational epsilon;

	public:
		PoleSequence(uint32_t type, uint32_t spin, const mp::rational& epsilon, bool include_odd = true)
		    : include_odd(include_odd), type(type), k(1), spin(spin), epsilon(epsilon)
		{
			assert(1 <= type && type <= 3);
			if (type != 2 && !include_odd) k = 2;
		}
		[[nodiscard]] bool valid() const { return type != 3 || k <= spin; }
		[[nodiscard]] mp::rational get() const
		{
			switch (type)
			{
			case 1: return mp::rational(-int32_t(spin + k - 1));
			case 2: return epsilon - (k - 1);
			default: return 2 * epsilon + (1 + spin - k);  // note: k <= spin -> 1 + spin - k >= 1
			}
		}
		void next() & { k += type == 2 || include_odd ? 1 : 2; }
	};
	template <class L, class R>
	class Merged
	{
		L seql;
		R seqr;
		bool next_l{};
		void update() & { next_l = !seqr.valid() || (seql.valid() && seql.get() >= seqr.get()); }

	public:
		Merged(L&& l, R&& r) : seql(std::move(l)), seqr(std::move(r)) { update(); }
		void next() &
		{
			if (next_l)
				seql.next();
			else
				seqr.next();
			update();
		}
		[[nodiscard]] bool valid() const { return seql.valid() || seqr.valid(); }
		[[nodiscard]] mp::rational get() const { return next_l ? seql.get() : seqr.get(); }
	};
	inline bool include_odd(const mp::real& d1, const mp::real& d2, const mp::real& d3, const mp::real& d4)
	{
		return d1 != d2 && d3 != d4;
	}
	inline bool include_odd(const mp::real& d12, const mp::real& d34) { return d12 != 0 && d34 != 0; }
	class ConformalScale : public ScaleFactor
	{
		bool odd_included_;
		uint32_t spin_, lambda_;
		const mp::rational* epsilon_;
		const mp::real* rho_;
		mp::real unitarity_bound_{}, gap_{};
		std::optional<mp::real> end_{};
		// pole at Delta = poles_[i]
		// 1 / (Delta - poles_[i])
		algebra::Vector<mp::real> poles_;
		std::optional<algebra::Vector<algebra::Polynomial>> bilinear_bases_{};
		void _set_bilinear_bases() &;

	public:
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, const mp::real& d1, const mp::real& d2,
		               const mp::real& d3, const mp::real& d4)
		    : ConformalScale(cutoff, spin, context, d1 - d2, d3 - d4)
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, const mp::real& delta12,
		               const mp::real& delta34)
		    : ConformalScale(cutoff, spin, context, include_odd(delta12, delta34))
		{
		}
		ConformalScale(const GeneralPrimaryOperator& op, const Context& context, bool include_odd)
		    : ConformalScale(op.num_poles(), op.spin(), context, include_odd)
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, bool include_odd);
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
			mp::real ans = mp::pow(4 * *rho_, delta);
			if (end_.has_value()) ans *= mp::pow(get_x(delta) + *end_, -int32_t(max_degree()));
			for (const auto& p : poles_) ans /= delta - p;
			return ans;
		}
		// convert x (in [0, \infty)) to delta
		[[nodiscard]] mp::real get_delta(const mp::real& x) const
		{
			if (end_.has_value()) return *end_ * (x + gap_) / (x + *end_);
			return x + gap_;
		}
		// convert delta to x (in [0, \infty))
		[[nodiscard]] mp::real get_x(const mp::real& delta) const
		{
			if (end_.has_value()) return *end_ * (delta - gap_) / (*end_ - gap_);
			return delta - gap_;
		}
		// evaluate at delta
		[[nodiscard]] mp::real eval(const mp::real& x) const override { return eval_d(get_delta(x)); }
		[[nodiscard]] algebra::Vector<mp::real> sample_scalings() override
		{
			auto xs = sample_points();
			for (auto& x : xs) x = eval(x);
			return xs;
		}
		// constrain delta to be in the range gap <= delta < end
		// if end is nullopt (corresponds to infinity), gap <= delta
		void set_gap(const mp::real& gap, const std::optional<mp::real>& end = {}) &
		{
			bilinear_bases_ = {};  // reset cache
			gap_ = gap;
			end_ = end;
		}
		[[nodiscard]] algebra::Polynomial bilinear_base(uint32_t m) override
		{
			_set_bilinear_bases();
			return bilinear_bases_->at(m).clone();
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial> bilinear_bases() override
		{
			_set_bilinear_bases();
			return bilinear_bases_.value().clone();
		}
		[[nodiscard]] mp::real sample_point(uint32_t k) override { return qboot::sample_point(max_degree(), k); }
		[[nodiscard]] algebra::Vector<mp::real> sample_points() override { return qboot::sample_points(max_degree()); }
	};
	inline std::unique_ptr<TrivialScale> get_scale([[maybe_unused]] const ConformalBlock<PrimaryOperator>& block,
	                                               [[maybe_unused]] uint32_t num_poles,
	                                               [[maybe_unused]] const Context& c)
	{
		return std::make_unique<TrivialScale>();
	}
	inline std::unique_ptr<ConformalScale> get_scale(const ConformalBlock<GeneralPrimaryOperator>& block,
	                                                 uint32_t num_poles, const Context& c)
	{
		auto scale = std::make_unique<ConformalScale>(num_poles, block.spin(), c, block.include_odd());
		const auto& op = block.get_op();
		scale->set_gap(op.lower_bound(), op.upper_bound_safe());
		return scale;
	}
}  // namespace qboot

#endif  // QBOOT_CONFORMAL_SCALE_HPP_
