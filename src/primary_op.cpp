#include "qboot/primary_op.hpp"
#include "qboot/context.hpp"

using std::string, std::ostringstream;

namespace qboot
{
	PrimaryOperator::PrimaryOperator(const Context& c) : PrimaryOperator(c.epsilon()) {}
	PrimaryOperator::PrimaryOperator(uint32_t spin, const Context& c) : PrimaryOperator(spin, c.epsilon()) {}
	PrimaryOperator::PrimaryOperator(const mp::real& delta, uint32_t spin, const Context& c)
	    : PrimaryOperator(delta, spin, c.epsilon())
	{
	}
	[[nodiscard]] string PrimaryOperator::str() const
	{
		ostringstream os;
		os << "operator(spin=" << spin_ << ", dim=" << 2 * epsilon_ + 2 << ", delta=" << delta_ << ")";
		return os.str();
	}
	[[nodiscard]] string GeneralPrimaryOperator::str() const
	{
		ostringstream os;
		os << "operator(spin=" << spin_ << ", dim=" << 2 * epsilon_ + 2 << ", delta in [" << from_ << ", "
		   << (to_.has_value() ? to_.value().str() : "infty") << "))";
		return os.str();
	}
}  // namespace qboot
