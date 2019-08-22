#ifndef PRIMARY_OP_HPP_
#define PRIMARY_OP_HPP_

#include <cstdint>  // for uint32_t
#include <sstream>  // for ostringstream
#include <string>   // for string

#include "real.hpp"  // for real, is_integer

namespace qboot
{
	// primary operator whose dimension is delta and spin is spin
	// delta may be a polynomial, especially just a monomial f(x) = x
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class PrimaryOperator
	{
		Real delta_, epsilon_;
		uint32_t spin_;

	public:
		PrimaryOperator(const Real& delta, uint32_t spin, const Real& epsilon)
		    : delta_(delta), epsilon_(epsilon), spin_(spin)
		{
		}
		PrimaryOperator get_shifted(const Real& small) const
		{
			return PrimaryOperator(delta_ == 0 ? small : delta_ * (1 + small), spin_, epsilon_);
		}
		const Real& delta() const { return delta_; }
		uint32_t spin() const { return spin_; }
		const Real& epsilon() const { return epsilon_; }
		Real quadratic_casimir() const
		{
			return (2 * epsilon_ + spin_) * spin_ / 2 + (delta_ / 2 - epsilon_ - 1) * delta_;
		}
		Real quartic_casimir() const
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
			os << "PrimaryOperator(delta=" << delta_ << ", spin=" << spin_ << ", epsilon=" << epsilon_ << ")";
			return os.str();
		}
	};
}  // namespace qboot

#endif  // PRIMARY_OP_HPP_
