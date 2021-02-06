#include "qboot/block.hpp"

#include <sstream>  // for ostringstream

using qboot::algebra::FunctionSymmetry;
using std::string, std::ostringstream;

namespace qboot
{
	[[maybe_unused]] static void _instantiate(ConformalBlock<PrimaryOperator>, ConformalBlock<GeneralPrimaryOperator>)
	{
	}

	template <_operator Operator>
	[[nodiscard]] string ConformalBlock<Operator>::str() const
	{
		ostringstream os;
		os << "ConformalBlock(op=" << op_.str() << ", d12=" << d12_ << ", d34=" << d34_ << ", d23h=" << d23h_
		   << ", type=" << (sym_ == FunctionSymmetry::Odd ? 'F' : 'H') << ")";
		return os.str();
	}

	[[nodiscard]] string GeneralConformalBlock::str() const
	{
		ostringstream os;
		os << "ConformalBlock(d12=" << d12_ << ", d34=" << d34_ << ", d23h=" << d23h_
		   << ", type=" << (sym_ == FunctionSymmetry::Odd ? 'F' : 'H') << ")";
		return os.str();
	}
}  // namespace qboot
