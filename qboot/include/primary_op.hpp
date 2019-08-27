#ifndef QBOOT_PRIMARY_OP_HPP_
#define QBOOT_PRIMARY_OP_HPP_

#include <cstdint>      // for uint32_t
#include <sstream>      // for ostringstream
#include <string>       // for string
#include <type_traits>  // for true_type

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
		[[nodiscard]] PrimaryOperator get_shifted(const Real& small) const
		{
			return PrimaryOperator(delta_ == 0 ? small : delta_ * (1 + small), spin_, epsilon_);
		}
		[[nodiscard]] const Real& delta() const { return delta_; }
		[[nodiscard]] uint32_t spin() const { return spin_; }
		[[nodiscard]] const Real& epsilon() const { return epsilon_; }
		[[nodiscard]] Real quadratic_casimir() const
		{
			return (2 * epsilon_ + spin_) * spin_ / 2 + (delta_ / 2 - epsilon_ - 1) * delta_;
		}
		[[nodiscard]] Real quartic_casimir() const
		{
			return spin_ * (spin_ + 2 * epsilon_) * (delta_ - 1) * (delta_ - 1 - 2 * epsilon_);
		}
		[[nodiscard]] bool is_divergent_hor() const
		{
			if (spin_ > 0) return delta_ == spin_ + 2 * epsilon_ || (mpfr::is_integer(delta_) && delta_ <= 1);
			return delta_ == epsilon_ || delta_ == epsilon_ + Real(0.5);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "PrimaryOperator(delta=" << delta_ << ", spin=" << spin_ << ", epsilon=" << epsilon_ << ")";
			return os.str();
		}
	};
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class GeneralPrimaryOperator
	{
		Real epsilon_;
		uint32_t spin_;

	public:
		GeneralPrimaryOperator(uint32_t spin, const Real& epsilon) : epsilon_(epsilon), spin_(spin) {}
		uint32_t spin() const { return spin_; }
		const Real& epsilon() const { return epsilon_; }
		[[nodiscard]] PrimaryOperator<Real> fix_delta(const Real& delta) const
		{
			return PrimaryOperator<Real>(delta, spin_, epsilon_);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "GeneralPrimaryOperator(spin=" << spin_ << ", epsilon=" << epsilon_ << ")";
			return os.str();
		}
	};
	template <class T>
	struct is_primary_operator;
	template <class T>
	constexpr bool is_primary_operator_v = is_primary_operator<T>::value;
	template <class Real>
	struct is_primary_operator<PrimaryOperator<Real>> : std::true_type
	{
	};
	template <class T>
	struct is_primary_operator : std::false_type
	{
	};
}  // namespace qboot

#endif  // QBOOT_PRIMARY_OP_HPP_
