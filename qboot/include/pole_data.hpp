#ifndef QBOOT_POLE_DATA_HPP_
#define QBOOT_POLE_DATA_HPP_

#include <cassert>   // for assert
#include <cstdint>   // for uint32_t, int32_t
#include <memory>    // for make_unique, unique_ptr
#include <optional>  // for optional
#include <utility>   // for move

#include "block.hpp"              // for ConformalBlock
#include "chol_and_inverse.hpp"   // for anti_band_to_inverse
#include "context_variables.hpp"  // for Context
#include "integral_decomp.hpp"    // for simple_pole_integral
#include "matrix.hpp"             // for Vector
#include "partial_fraction.hpp"   // for fast_partial_fraction
#include "polynomial.hpp"         // for Polynomial
#include "real.hpp"               // for real, pow, log
#include "scale_factor.hpp"       // for ScaleFactor

namespace qboot
{
	inline mpfr::real sample_point(uint32_t degree, uint32_t k)
	{
		return mpfr::pow(mpfr::const_pi() * (-1 + 4 * int32_t(k)), 2uL) /
		       ((-64 * int32_t(degree + 1)) * mpfr::log(3 - mpfr::sqrt(8)));
	}
	inline algebra::Vector<mpfr::real> sample_points(uint32_t degree)
	{
		algebra::Vector<mpfr::real> v(degree + 1);
		for (uint32_t i = 0; i <= degree; ++i) v[i] = mpfr::pow(mpfr::const_pi() * (-1 + 4 * int32_t(i)), 2uL);
		v /= (-64 * int32_t(degree + 1)) * mpfr::log(3 - mpfr::sqrt(8));
		return v;
	}
	class PoleSequence
	{
		bool include_odd;
		uint32_t type, k, spin;
		mpfr::real epsilon;

	public:
		PoleSequence(uint32_t type, uint32_t spin, const mpfr::real& epsilon, bool include_odd = true)
		    : include_odd(include_odd), type(type), k(1), spin(spin), epsilon(epsilon)
		{
			assert(1 <= type && type <= 3);
			if (type != 2 && !include_odd) k = 2;
		}
		[[nodiscard]] bool valid() const { return type != 3 || k <= spin; }
		[[nodiscard]] mpfr::real get() const
		{
			switch (type)
			{
			case 1: return mpfr::real(-int32_t(spin + k - 1));
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
		[[nodiscard]] mpfr::real get() const { return next_l ? seql.get() : seqr.get(); }
	};
	inline bool include_odd(const mpfr::real& d1, const mpfr::real& d2, const mpfr::real& d3, const mpfr::real& d4)
	{
		return d1 != d2 && d3 != d4;
	}
	inline bool include_odd(const mpfr::real& d12, const mpfr::real& d34) { return d12 != 0 && d34 != 0; }
	class ConformalScale : public ScaleFactor
	{
		bool odd_included_;
		uint32_t spin_, lambda_;
		const mpfr::real &epsilon_, &rho_;
		mpfr::real unitarity_bound_{}, gap_{};
		std::optional<mpfr::real> end_{};
		// pole at Delta = poles_[i]
		// 1 / (Delta - poles_[i])
		algebra::Vector<mpfr::real> poles_;
		std::optional<algebra::Vector<algebra::Polynomial>> bilinear_bases_{};
		void _set_bilinear_bases() &
		{
			if (bilinear_bases_.has_value()) return;
			if (end_.has_value())
			{
				bilinear_bases_ = algebra::Vector<algebra::Polynomial>{max_degree() / 2 + 1};
				for (uint32_t i = 0; i < bilinear_bases_->size(); ++i) bilinear_bases_->at(i) = algebra::Polynomial(i);
			}
			else
			{
				// orthogonal polynomial of weight function (4 rho) ^ {Delta} / \prod_i (Delta - poles[i])
				algebra::Vector<mpfr::real> shifted_poles(poles_.size());
				for (uint32_t i = 0; i < poles_.size(); ++i) shifted_poles[i] = poles_[i] - gap_;
				auto weight = fast_partial_fraction(shifted_poles);
				auto deg = max_degree() / 2;
				// inner_prods[i] = \int_{0}^{\infty} dx (4 rho) ^ x x ^ i / \prod_i (x - poles[i])
				algebra::Vector<mpfr::real> inner_prods(2 * deg + 1);
				for (uint32_t i = 0; i < poles_.size(); ++i)
					inner_prods += mul_scalar(weight[i], simple_pole_integral(2 * deg, 4 * rho_, shifted_poles[i]));
				inner_prods *= mpfr::pow(4 * rho_, gap_);
				auto mat = anti_band_to_inverse(inner_prods);
				bilinear_bases_ = algebra::Vector<algebra::Polynomial>{deg + 1};
				for (uint32_t i = 0; i <= deg; ++i)
				{
					algebra::Vector<mpfr::real> v(i + 1);
					for (uint32_t j = 0; j <= i; ++j) v[j] = mat.at(i, j);
					bilinear_bases_->at(i) = algebra::Polynomial(std::move(v));
				}
			}
		}

	public:
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, const mpfr::real& d1,
		               const mpfr::real& d2, const mpfr::real& d3, const mpfr::real& d4)
		    : ConformalScale(cutoff, spin, context, d1 - d2, d3 - d4)
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, const mpfr::real& delta12,
		               const mpfr::real& delta34)
		    : ConformalScale(cutoff, spin, context, include_odd(delta12, delta34))
		{
		}
		ConformalScale(const GeneralPrimaryOperator& op, const Context& context, bool include_odd)
		    : ConformalScale(op.num_poles(), op.spin(), context, include_odd)
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, bool include_odd)
		    : odd_included_(include_odd),
		      spin_(spin),
		      lambda_(context.lambda()),
		      epsilon_(context.epsilon()),
		      rho_(context.rho()),
		      poles_(cutoff)
		{
			// type 1 or 3 PoleData vanishes when delta12 == 0 or delta34 == 0 and k is odd
			auto get_pols = [this, include_odd](uint32_t type) {
				return PoleSequence(type, this->spin_, this->epsilon_, include_odd);
			};
			auto pole_seq = Merged(Merged(get_pols(1), get_pols(2)), get_pols(3));
			uint32_t pos = 0;
			while (pos < cutoff)
			{
				poles_[pos++] = pole_seq.get();
				pole_seq.next();
			}
		}
		ConformalScale(const ConformalScale&) = delete;
		ConformalScale& operator=(const ConformalScale&) = delete;
		ConformalScale(ConformalScale&&) noexcept = default;
		ConformalScale& operator=(ConformalScale&&) noexcept = default;
		~ConformalScale() override = default;
		[[nodiscard]] bool odd_included() const { return odd_included_; }
		[[nodiscard]] uint32_t max_degree() const override { return poles_.size() + lambda_; }
		[[nodiscard]] const algebra::Vector<mpfr::real>& get_poles() const& { return poles_; }
		[[nodiscard]] algebra::Vector<mpfr::real> get_poles() && { return std::move(poles_); }
		// evaluate at delta
		[[nodiscard]] mpfr::real eval_d(const mpfr::real& delta) const
		{
			mpfr::real ans = mpfr::pow(4 * rho_, delta);
			if (end_.has_value()) ans *= mpfr::pow(get_x(delta) + *end_, -int32_t(max_degree()));
			for (uint32_t i = 0; i < poles_.size(); ++i) ans /= delta - poles_[i];
			return ans;
		}
		// convert x (in [0, \infty)) to delta
		[[nodiscard]] mpfr::real get_delta(const mpfr::real& x) const
		{
			if (end_.has_value()) return *end_ * (x + gap_) / (x + *end_);
			return x + gap_;
		}
		// convert delta to x (in [0, \infty))
		[[nodiscard]] mpfr::real get_x(const mpfr::real& delta) const
		{
			if (end_.has_value()) return *end_ * (delta - gap_) / (*end_ - gap_);
			return delta - gap_;
		}
		// evaluate at delta
		[[nodiscard]] mpfr::real eval(const mpfr::real& x) const override { return eval_d(get_delta(x)); }
		[[nodiscard]] algebra::Vector<mpfr::real> sample_scalings() override
		{
			auto xs = sample_points();
			for (uint32_t i = 0; i < xs.size(); ++i) xs[i] = eval(xs[i]);
			return xs;
		}
		// constrain delta to be in the range gap <= delta < end
		// if end is nullopt (corresponds to infinity), gap <= delta
		void set_gap(const mpfr::real& gap, std::optional<mpfr::real> end = {}) &
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
		[[nodiscard]] mpfr::real sample_point(uint32_t k) override { return qboot::sample_point(max_degree(), k); }
		[[nodiscard]] algebra::Vector<mpfr::real> sample_points() override
		{
			return qboot::sample_points(max_degree());
		}
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

#endif  // QBOOT_CONTEXT_VARIABLES_HPP_
