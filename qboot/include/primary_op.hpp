#ifndef QBOOT_PRIMARY_OP_HPP_
#define QBOOT_PRIMARY_OP_HPP_

#include <cstdint>      // for uint32_t
#include <optional>     // for optional
#include <sstream>      // for ostringstream
#include <string>       // for string
#include <type_traits>  // for true_type
#include <variant>      // for variant

#include "real.hpp"  // for real, isinteger

namespace qboot
{
	template <class Real>
	Real unitarity_bound(const Real& epsilon, uint32_t spin)
	{
		return spin == 0 ? epsilon : spin + 2 * epsilon;
	}
	// primary operator whose dimension is delta and spin is spin
	template <class Real>
	class PrimaryOperator
	{
		Real delta_, epsilon_;
		uint32_t spin_;

	public:
		// unit operator
		explicit PrimaryOperator(const Real& epsilon) : delta_{}, epsilon_(epsilon), spin_(0) {}
		// on the unitarity bound
		PrimaryOperator(uint32_t spin, const Real& epsilon)
		    : delta_(unitarity_bound(epsilon, spin)), epsilon_(epsilon), spin_(spin)
		{
		}
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
			if (spin_ > 0) return delta_ == spin_ + 2 * epsilon_ || (mpfr::isinteger(delta_) && delta_ <= 1);
			return delta_.iszero() || delta_ == epsilon_ || delta_ == epsilon_ + Real(0.5);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "PrimaryOperator(delta=" << delta_ << ", spin=" << spin_ << ", epsilon=" << epsilon_ << ")";
			return os.str();
		}
	};
	template <class Real>
	class GeneralPrimaryOperator
	{
		Real epsilon_, from_;
		std::optional<Real> to_;
		uint32_t spin_;

	public:
		// primary operator whose dimension runs over unitarity bound
		GeneralPrimaryOperator(uint32_t spin, const Real& epsilon)
		    : epsilon_(epsilon), from_(unitarity_bound(epsilon, spin)), to_{}, spin_(spin)
		{
		}
		// primary operator whose dimension runs from lb to ub (empty ub means infinity)
		GeneralPrimaryOperator(uint32_t spin, const Real& epsilon, const Real& lb, std::optional<Real> ub = {})
		    : epsilon_(epsilon), from_(lb), to_(ub), spin_(spin)
		{
		}
		[[nodiscard]] uint32_t spin() const { return spin_; }
		[[nodiscard]] const Real& epsilon() const { return epsilon_; }
		[[nodiscard]] const Real& lower_bound() const { return from_; }
		[[nodiscard]] bool is_finite() const noexcept { return to_.has_value(); }
		// if not is_finite(), do not call this
		[[nodiscard]] const Real& upper_bound() const { return *to_; }
		// if not is_finite(), returns nullopt
		[[nodiscard]] const std::optional<Real>& upper_bound_safe() const { return to_; }
		[[nodiscard]] PrimaryOperator<Real> fix_delta(const Real& delta) const
		{
			assert(from_ <= delta && (!to_.has_value() || delta < *to_));
			return PrimaryOperator<Real>(delta, spin_, epsilon_);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "GeneralPrimaryOperator(spin=" << spin_ << ", epsilon=" << epsilon_ << ", delta in [" << from_ << ", "
			   << (to_.has_value() ? to_->str() : "infty") << "))";
			return os.str();
		}
	};
	template <class T>
	struct is_primary_operator;
	template <class T>
	inline constexpr bool is_primary_operator_v = is_primary_operator<T>::value;
	template <class Real>
	struct is_primary_operator<PrimaryOperator<Real>> : std::true_type
	{
	};
	template <class T>
	struct is_primary_operator : std::false_type
	{
	};
	template <class Real>
	using Operator = std::variant<PrimaryOperator<Real>, GeneralPrimaryOperator<Real>>;
}  // namespace qboot

#endif  // QBOOT_PRIMARY_OP_HPP_
