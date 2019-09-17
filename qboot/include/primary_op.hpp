#ifndef QBOOT_PRIMARY_OP_HPP_
#define QBOOT_PRIMARY_OP_HPP_

#include <cassert>      // for assert
#include <cstdint>      // for uint32_t
#include <optional>     // for optional
#include <sstream>      // for ostringstream
#include <string>       // for string
#include <type_traits>  // for true_type
#include <variant>      // for variant

#include "rational.hpp"  // for rational
#include "real.hpp"      // for real, isinteger

namespace qboot
{
	inline bool _is_positive_integer(const mpfr::real& x) { return x > 0 && mpfr::isinteger(x); }
	inline mpfr::rational unitarity_bound(const mpfr::rational& epsilon, uint32_t spin)
	{
		return spin == 0 ? epsilon : spin + 2 * epsilon;
	}
	class Context;
	// primary operator whose dimension is delta and spin is spin
	class PrimaryOperator
	{
		mpfr::real delta_{};
		mpfr::rational epsilon_;
		uint32_t spin_;

	public:
		// unit operator
		explicit PrimaryOperator(const Context& c);
		// unit operator
		explicit PrimaryOperator(const mpfr::rational& epsilon) : epsilon_(epsilon), spin_(0) {}
		// on the unitarity bound
		PrimaryOperator(uint32_t spin, const Context& c);
		// on the unitarity bound
		PrimaryOperator(uint32_t spin, const mpfr::rational& epsilon)
		    : delta_(unitarity_bound(epsilon, spin)), epsilon_(epsilon), spin_(spin)
		{
		}
		PrimaryOperator(const mpfr::real& delta, uint32_t spin, const Context& c);
		PrimaryOperator(const mpfr::real& delta, uint32_t spin, const mpfr::rational& epsilon)
		    : delta_(delta), epsilon_(epsilon), spin_(spin)
		{
		}
		[[nodiscard]] PrimaryOperator get_shifted(const mpfr::real& small) const
		{
			return PrimaryOperator(delta_ == 0 ? small : delta_ * (1 + small), spin_, epsilon_);
		}
		[[nodiscard]] const mpfr::real& delta() const { return delta_; }
		[[nodiscard]] uint32_t spin() const { return spin_; }
		[[nodiscard]] const mpfr::rational& epsilon() const { return epsilon_; }
		[[nodiscard]] mpfr::real quadratic_casimir() const
		{
			return (2 * epsilon_ + spin_) * spin_ / 2 + (delta_ / 2 - epsilon_ - 1) * delta_;
		}
		[[nodiscard]] mpfr::real quartic_casimir() const
		{
			return spin_ * (spin_ + 2 * epsilon_) * (delta_ - 1) * (delta_ - 1 - 2 * epsilon_);
		}
		[[nodiscard]] bool is_divergent_hor() const
		{
			// if spin > 0, exists n > 0
			// s.t. 0 = -n (n + 2 Delta - 2 epsilon - 2) (n + Delta - 2 epsilon - spin - 1) (n + Delta + spin - 1)
			// if spin = 0, exists n > 0
			// s.t. 0 = n (n + Delta - 1) (n + 2 Delta - 2 epsilon - 2)
			return _is_positive_integer((1 - int32_t(spin_)) - delta_) ||
			       _is_positive_integer(2 + 2 * epsilon_ - 2 * delta_) ||
			       (spin_ > 0 && _is_positive_integer((1 + spin_) + 2 * epsilon_ - delta_));
		}
		[[nodiscard]] std::string str() const;
	};
	class GeneralPrimaryOperator
	{
		mpfr::rational epsilon_;
		mpfr::real from_;
		std::optional<mpfr::real> to_{};
		uint32_t spin_, num_poles_;

	public:
		// primary operator whose dimension runs over unitarity bound
		GeneralPrimaryOperator(uint32_t spin, uint32_t num_poles, const mpfr::rational& epsilon)
		    : epsilon_(epsilon), from_(unitarity_bound(epsilon, spin)), spin_(spin), num_poles_(num_poles)
		{
		}
		// primary operator whose dimension runs from lb to ub (empty ub means infinity)
		GeneralPrimaryOperator(uint32_t spin, uint32_t num_poles, const mpfr::rational& epsilon, const mpfr::real& lb,
		                       const std::optional<mpfr::real>& ub = {})
		    : epsilon_(epsilon), from_(lb), to_(ub), spin_(spin), num_poles_(num_poles)
		{
		}
		[[nodiscard]] uint32_t spin() const { return spin_; }
		[[nodiscard]] uint32_t num_poles() const { return num_poles_; }
		[[nodiscard]] const mpfr::rational& epsilon() const { return epsilon_; }
		[[nodiscard]] const mpfr::real& lower_bound() const { return from_; }
		[[nodiscard]] bool is_finite() const noexcept { return to_.has_value(); }
		// if not is_finite(), do not call this
		[[nodiscard]] const mpfr::real& upper_bound() const { return *to_; }
		// if not is_finite(), returns nullopt
		[[nodiscard]] const std::optional<mpfr::real>& upper_bound_safe() const { return to_; }
		[[nodiscard]] PrimaryOperator fix_delta(const mpfr::real& delta) const
		{
			assert(from_ <= delta && (!to_.has_value() || delta < *to_));
			return PrimaryOperator(delta, spin_, epsilon_);
		}
		[[nodiscard]] std::string str() const;
	};
}  // namespace qboot

#endif  // QBOOT_PRIMARY_OP_HPP_
