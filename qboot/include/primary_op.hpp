#ifndef PRIMARY_OP_HPP_
#define PRIMARY_OP_HPP_

#include <cstdint>  // for uint32_t
#include <string>   //for string

#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for real

namespace qboot
{
	template <class Real>
	class Context;
	// primary operator whose dimension is delta and spin is spin
	// delta may be a polynomial, especially just a monomial f(x) = x
	template <class Real = mpfr::real<1000, MPFR_RNDN>, class T = Real>
	class PrimaryOperator
	{
		T delta_;
		Real epsilon_;
		uint32_t spin_;
		const Context<Real>& context_;

	public:
		PrimaryOperator(const T& delta, uint32_t spin, const Context<Real>& context)
		    : delta_(delta.clone()), epsilon_(context.epsilon), spin_(spin), context_(context)
		{
		}
		const T& delta() const { return delta_; }
		uint32_t spin() const { return spin_; }
		const Real& epsilon() const { return epsilon_; }
		const Context<Real>& context() const { return context_; }
		T quadratic_casimir() const
		{
			return (2 * epsilon_ + spin_) * spin_ / 2 + (delta_ / 2 - epsilon_ - 1) * delta_;
		}
		T quartic_casimir() const
		{
			return spin_ * (spin_ + 2 * epsilon_) * (delta_ - 1) * (delta_ - 1 - 2 * epsilon_);
		}
		bool is_divergent_hor() const
		{
			if (spin_ > 0) return delta_ == spin_ + 2 * epsilon_ || (mpfr::is_integer(delta_) && delta_ <= 1);
			return delta_ == epsilon_ || delta_ == epsilon_ + Real(0.5);
		}
		std::string str() const
		{
			std::ostringstream os;
			os << "PrimaryOperator(delta=" << delta_ << ", spin=" << spin_ << ", context=" << context_.str() << ")";
			return os.str();
		}
	};
	template <class Real = mpfr::real<1000, MPFR_RNDN>, class T = Real>
	auto general_primary_operator(uint32_t spin, const Context<Real>& context)
	{
		return PrimaryOperator(algebra::Polynomial<Real>{Real(0), Real(1)}, spin, context);
	}
}  // namespace qboot

#endif  // PRIMARY_OP_HPP_
