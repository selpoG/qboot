#include "primary_op.hpp"

using std::string, std::ostringstream;

namespace qboot
{
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
