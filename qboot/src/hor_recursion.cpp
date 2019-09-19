#include "hor_recursion.hpp"

#include <utility>  // for move

#include "context.hpp"      // for Context
#include "hor_formula.hpp"  // for _get_rec_coeffs

using algebra::RealFunction, algebra::RealFunctionWithPower, algebra::ComplexFunction;
using mp::real;
using std::move;

namespace qboot
{
	RealFunction<real> hBlock_shifted(const PrimaryOperator& op, const real& S, const real& P, uint32_t n_Max)
	{
		RealFunction<real> b(n_Max);
		b.at(0) = mp::pochhammer(op.epsilon() * 2, op.spin()) /
		          (mp::pow(2, op.spin()) * mp::pochhammer(op.epsilon(), op.spin()));
		if (op.spin() % 2 == 1) b.at(0).negate();
		if (op.delta().iszero()) return b;
		if (op.is_divergent_hor())
		{
			real small(1ul, -(3 * mp::global_prec) / 8 + 15);
			b = hBlock_shifted(op.get_shifted(small), S, P, n_Max);
			b += hBlock_shifted(op.get_shifted(-small), S, P, n_Max);
			b /= 2;
			return b;
		}
		real sum;
		auto p = _get_rec_coeffs(op, S, P);
		for (uint32_t n = 1; n <= n_Max; ++n)
		{
			sum = 0;
			for (uint32_t i = 1; i < p.size(); ++i)
				if (i <= n) sum += p[i].eval(n) * b.at(n - i);
			b.at(n) = -sum / p[0].eval(n);
		}
		return b;
	}
	RealFunction<real> gBlock_real(const PrimaryOperator& op, const real& S, const real& P, const Context& context)
	{
		auto h_at_0 = hBlock_shifted(op, S, P, context.n_Max());
		h_at_0 *= mp::pow(4, op.delta());
		const auto& rho = context.rho();
		RealFunctionWithPower f_at_0(move(h_at_0), op.delta());
		RealFunction<real> f_of_rho(context.lambda());
		f_of_rho.at(0) = f_at_0.approximate(rho);
		real tmp(1);
		for (uint32_t k = 1; k <= f_of_rho.lambda(); ++k)
		{
			tmp *= k;
			f_at_0.derivate();
			f_of_rho.at(k) = f_at_0.approximate(rho) / tmp;
		}
		return context.rho_to_z().convert(f_of_rho);
	}
	ComplexFunction<real> gBlock(const PrimaryOperator& op, const real& S, const real& P, const Context& context)
	{
		return context.expand_off_diagonal(gBlock_real(op, S, P, context), op, S, P);
	}
	ComplexFunction<real> gBlock(const PrimaryOperator& op, const real& d1, const real& d2, const real& d3,
	                             const real& d4, const Context& context)
	{
		real d12 = d1 - d2, d34 = d3 - d4;
		return gBlock(op, (d34 - d12) / 2, -d12 * d34 / 2, context);
	}
}  // namespace qboot
