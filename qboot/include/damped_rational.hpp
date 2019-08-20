#ifndef DAMPED_RATIONAL_
#define DAMPED_RATIONAL_

#include <cassert>
#include <cstddef>
#include <map>
#include <optional>

#include "real.hpp"

namespace qboot
{
	using std::size_t;
	// x -> (p * x + q) / (r * x + s)
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class MobiusTransformation
	{
		Real p_, q_, r_, s_;
		MobiusTransformation(const Real& p, const Real& q, const Real& r, const Real& s) : p_(p), q_(q), r_(r), s_(s)
		{
			assert(p * s - q * r != 0);
		}

	public:
		Real eval(const Real& x) const { return (p_ * x + q_) / (r_ * x + s_); }
		std::optional<Real> safe_eval(const std::optional<Real>& x) const
		{
			if (r_ == 0) return x.has_value() ? (p_ * x.value() + q_) / s_ : x;
			if (!x.has_value()) return p_ / r_;
			auto den = r_ * x.value() + s_;
			if (den == 0) return std::nullopt;
			return (p_ * x.value() + q_) / den;
		}
		std::optional<Real> get_shift() const noexcept
		{
			return p_ == s_ && p_ == 1 && r_ == 0 ? std::optional(q_) : std::nullopt;
		}
		std::optional<Real> get_scale() const noexcept
		{
			return q_ == r_ && q_ == 0 && s_ == 1 ? std::optional(p_) : std::nullopt;
		}
		bool is_flip() const noexcept { return p_ == s_ && p_ == 0 && q_ == r_ && q_ == 1; }
		// x -> x + a, [0, infty) -> [a, infty)
		static MobiusTransformation Shift(const Real& a) { return MobiusTransformation(Real(1), a, Real(0), Real(1)); }
		// x -> a * x, [0, infty) -> [0, infty)
		static MobiusTransformation Scale(const Real& a) { return MobiusTransformation(a, Real(0), Real(0), Real(1)); }
		// x -> 1 / x, [0, infty) -> (0, infty]
		static MobiusTransformation Flip() { return MobiusTransformation(Real(0), Real(1), Real(1), Real(0)); }
		// x -> b * (x + a) / (x + b), [0, infty) -> [a, b)
		static MobiusTransformation Range(const Real& a, const Real& b)
		{
			return MobiusTransformation(b, a * b, Real(1), b);
		}
		// x -> (a * x + b) / (c * x + d), [0, infty) -> [b / d, a / c) if c * d >= 0
		static MobiusTransformation General(const Real& a, const Real& b, const Real& c, const Real& d)
		{
			return MobiusTransformation(a, b, c, d);
		}
	};
	// function of x
	// pow(base_, k / (c - x)) (if c has value)
	// pow(base_, k * x) (otherwise)
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class ExpLike
	{
		Real base_, k_;
		std::optional<Real> c_;
	};
	// function of x, c * e(x) / (\\prod_p pow(x - p, poles[p]))
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class DampedRational
	{
		Real c_;
		ExpLike<Real> e_;
		std::map<Real, int> poles_;

	public:
		DampedRational(const Real& c, const ExpLike<Real>& e, const std::map<Real, int>& poles)
		    : c_(c), e_(e), poles_(poles)
		{
		}
	};
}  // namespace qboot

#endif  // DAMPED_RATIONAL_
