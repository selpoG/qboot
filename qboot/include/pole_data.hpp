#ifndef QBOOT_POLE_DATA_HPP_
#define QBOOT_POLE_DATA_HPP_

#include <cassert>   // for assert
#include <cstdint>   // for uint32_t, int32_t
#include <optional>  // for optional
#include <utility>   // for move

#include "context_variables.hpp"  // for Context
#include "matrix.hpp"             // for Vector
#include "polynomial.hpp"         // for Polynomial
#include "real.hpp"               // for real, pow, log
#include "scale_factor.hpp"       // for ScaleFactor

namespace qboot
{
	template <class Real>
	Real sample_point(uint32_t degree, uint32_t k)
	{
		return mpfr::pow(Real::pi() * (-1 + 4 * int32_t(k)), 2uL) /
		       ((-64 * int32_t(degree + 1)) * mpfr::log(3 - Real::sqrt(8)));
	}
	template <class Real>
	algebra::Vector<Real> sample_points(uint32_t degree)
	{
		algebra::Vector<Real> v(degree + 1);
		for (uint32_t i = 0; i <= degree; ++i) v[i] = mpfr::pow(Real::pi() * (-1 + 4 * int32_t(i)), 2uL);
		v /= (-64 * int32_t(degree + 1)) * mpfr::log(3 - Real::sqrt(8));
		return v;
	}
	template <class Real>
	class PoleSequence
	{
		bool include_odd;
		uint32_t type, k, spin;
		Real epsilon;

	public:
		using value_type = Real;
		PoleSequence(uint32_t type, uint32_t spin, const Real& epsilon, bool include_odd = true)
		    : include_odd(include_odd), type(type), k(1), spin(spin), epsilon(epsilon)
		{
			assert(1 <= type && type <= 3);
			if (type != 2 && !include_odd) k = 2;
		}
		[[nodiscard]] bool valid() const { return type != 3 || k <= spin; }
		[[nodiscard]] Real get() const
		{
			switch (type)
			{
			case 1: return Real(-int32_t(spin + k - 1));
			case 2: return epsilon - (k - 1);
			default: return 2 * epsilon + (1 + spin - k);  // note: k <= spin -> 1 + spin - k >= 1
			}
		}
		void next() & { k += type == 2 || include_odd ? 1 : 2; }
	};
	template <class L, class R, class Real = typename L::value_type>
	class Merged
	{
		L seql;
		R seqr;
		bool next_l{};
		void update() & { next_l = !seqr.valid() || (seql.valid() && seql.get() >= seqr.get()); }

	public:
		using value_type = Real;
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
		[[nodiscard]] Real get() const { return next_l ? seql.get() : seqr.get(); }
	};
	template <class Real>
	bool include_odd(const Real& d1, const Real& d2, const Real& d3, const Real& d4)
	{
		return d1 != d2 && d3 != d4;
	}
	template <class Real>
	bool include_odd(const Real& d12, const Real& d34)
	{
		return d12 != 0 && d34 != 0;
	}
	template <class Real>
	class ConformalScale : public ScaleFactor<Real>
	{
		uint32_t spin_, lambda_;
		const Real &epsilon_, &rho_;
		Real unitarity_bound_{}, gap_{};
		std::optional<Real> end_{};
		// pole at Delta = poles_[i]
		// 1 / (Delta - poles_[i])
		algebra::Vector<Real> poles_;
		std::optional<algebra::Vector<algebra::Polynomial<Real>>> bilinear_bases_{};
		void _set_bilinear_bases() &
		{
			if (bilinear_bases_.has_value()) return;
			if (end_.has_value())
			{
				bilinear_bases_ = algebra::Vector<algebra::Polynomial<Real>>{max_degree() / 2 + 1};
				for (uint32_t i = 0; i < bilinear_bases_->size(); ++i)
					bilinear_bases_->at(i) = algebra::Polynomial<Real>(i);
			}
			else
			{
				// orthogonal polynomial of weight function (4 rho) ^ {Delta} / \prod_i (Delta - poles[i])
				algebra::Vector<Real> shifted_poles(poles_.size());
				for (uint32_t i = 0; i < poles_.size(); ++i) shifted_poles[i] = poles_[i] - gap_;
				auto weight = fast_partial_fraction(shifted_poles);
				auto deg = max_degree() / 2;
				// inner_prods[i] = \int_{0}^{\infty} dx (4 rho) ^ x x ^ i / \prod_i (x - poles[i])
				algebra::Vector<Real> inner_prods(2 * deg + 1);
				for (uint32_t i = 0; i < poles_.size(); ++i)
					inner_prods += mul_scalar(weight[i], simple_pole_integral(2 * deg, 4 * rho_, shifted_poles[i]));
				inner_prods *= mpfr::pow(4 * rho_, gap_);
				auto mat = anti_band_to_inverse(inner_prods);
				bilinear_bases_ = algebra::Vector<algebra::Polynomial<Real>>{deg + 1};
				for (uint32_t i = 0; i <= deg; ++i)
				{
					algebra::Vector<Real> v(i + 1);
					for (uint32_t j = 0; j <= i; ++j) v[j] = mat.at(i, j);
					bilinear_bases_->at(i) = algebra::Polynomial<Real>(std::move(v));
				}
			}
		}

	public:
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& d1, const Real& d2,
		               const Real& d3, const Real& d4)
		    : ConformalScale(cutoff, spin, context, d1 - d2, d3 - d4)
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& delta12,
		               const Real& delta34)
		    : ConformalScale(cutoff, spin, context, include_odd(delta12, delta34))
		{
		}
		ConformalScale(uint32_t cutoff, uint32_t spin, const Context<Real>& context, bool include_odd)
		    : spin_(spin), lambda_(context.lambda()), epsilon_(context.epsilon()), rho_(context.rho()), poles_(cutoff)
		{
			// type 1 or 3 PoleData vanishes when delta12 == 0 or delta34 == 0 and k is odd
			auto get_pols = [this, include_odd](uint32_t type) {
				return PoleSequence<Real>(type, this->spin_, this->epsilon_, include_odd);
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
		[[nodiscard]] uint32_t max_degree() const override { return poles_.size() + lambda_; }
		[[nodiscard]] const algebra::Vector<Real>& get_poles() const& { return poles_; }
		[[nodiscard]] algebra::Vector<Real> get_poles() && { return std::move(poles_); }
		// evaluate at delta
		[[nodiscard]] Real eval_d(const Real& delta) const
		{
			Real ans = mpfr::pow(4 * rho_, delta);
			if (end_.has_value()) ans *= mpfr::pow(get_x(delta) + *end_, -int32_t(max_degree()));
			for (uint32_t i = 0; i < poles_.size(); ++i) ans /= delta - poles_[i];
			return ans;
		}
		// convert x (in [0, \infty)) to delta
		[[nodiscard]] Real get_delta(const Real& x) const
		{
			if (end_.has_value()) return *end_ * (x + gap_) / (x + *end_);
			return x + gap_;
		}
		// convert delta to x (in [0, \infty))
		[[nodiscard]] Real get_x(const Real& delta) const
		{
			if (end_.has_value()) return *end_ * (delta - gap_) / (*end_ - gap_);
			return delta - gap_;
		}
		// evaluate at delta
		[[nodiscard]] Real eval(const Real& x) const override { return eval_d(get_delta(x)); }
		[[nodiscard]] algebra::Vector<Real> sample_scalings() override
		{
			auto xs = sample_points();
			for (uint32_t i = 0; i < xs.size(); i++) xs[i] = eval(xs[i]);
			return xs;
		}
		// constrain delta to be in the range gap <= delta < end
		// if end is nullopt (corresponds to infinity), gap <= delta
		void set_gap(const Real& gap, std::optional<Real> end = {}) &
		{
			bilinear_bases_ = {};  // reset cache
			gap_ = gap;
			end_ = end;
		}
		[[nodiscard]] algebra::Polynomial<Real> bilinear_base(uint32_t m) override
		{
			_set_bilinear_bases();
			return bilinear_bases_->at(m).clone();
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial<Real>> bilinear_bases() override
		{
			_set_bilinear_bases();
			return bilinear_bases_.value().clone();
		}
		[[nodiscard]] Real sample_point(uint32_t k) override { return qboot::sample_point<Real>(max_degree(), k); }
		[[nodiscard]] algebra::Vector<Real> sample_points() override
		{
			return qboot::sample_points<Real>(max_degree());
		}
	};
}  // namespace qboot

#endif  // QBOOT_CONTEXT_VARIABLES_HPP_
