#include "qboot/algebra/real_function.hpp"

using qboot::mp::real;

namespace qboot::algebra
{
	template class RealFunction<real>;
	template class RealFunction<Polynomial>;

	RealFunction<real> power_function(const real& a, const real& b, const real& p, uint32_t lambda)
	{
		assert(a != 0);
		RealFunction<real> f(lambda);
		f.at(0) = mp::pow(a, p);
		real tmp = b / a;
		for (uint32_t k = 1; k <= lambda; ++k) { f.at(k) = f.at(k - 1) * tmp * (p - (k - 1)) / k; }
		return f;
	}

	RealConverter::RealConverter(const RealFunction<real>& func)
	    : lambda_(func.lambda()), mat_(func.lambda() + 1, func.lambda() + 1)
	{
		// assert(func.at(0).iszero());
		mat_.at(0, 0) = 1;
		auto pf = func.clone();
		for (uint32_t n = 1; n <= lambda_; ++n)
		{
			for (uint32_t i = n; i <= lambda_; ++i) mat_.at(n, i) = pf.at(i);
			pf = mul(pf, func);
		}
	}
	[[nodiscard]] RealConverter RealConverter::inverse() const
	{
		auto invmat = mat_.clone();
		invmat.transpose();
		invmat = algebra::lower_triangular_inverse(invmat);
		RealFunction<real> f(lambda_);
		for (uint32_t k = 1; k <= lambda_; ++k) f.at(k) = invmat.at(k, 1);
		return RealConverter(f);
	}
}  // namespace qboot::algebra
