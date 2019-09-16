#include "primary_op.hpp"
#include "context.hpp"

using std::string, std::ostringstream;

namespace qboot
{
	PrimaryOperator::PrimaryOperator(const Context& c) : PrimaryOperator(c.epsilon()) {}
	PrimaryOperator::PrimaryOperator(uint32_t spin, const Context& c) : PrimaryOperator(spin, c.epsilon()) {}
	PrimaryOperator::PrimaryOperator(const mpfr::real& delta, uint32_t spin, const Context& c)
	    : PrimaryOperator(delta, spin, c.epsilon())
	{
	}
	[[nodiscard]] string PrimaryOperator::str() const
	{
		ostringstream os;
		os << "PrimaryOperator(delta=" << delta_ << ", spin=" << spin_ << ", epsilon=" << epsilon_ << ")";
		return os.str();
	}
	[[nodiscard]] string GeneralPrimaryOperator::str() const
	{
		ostringstream os;
		os << "GeneralPrimaryOperator(spin=" << spin_ << ", epsilon=" << epsilon_ << ", delta in [" << from_ << ", "
		   << (to_.has_value() ? to_->str() : "infty") << "))";
		return os.str();
	}
}  // namespace qboot
