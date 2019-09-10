#include "hor_recursion.hpp"

#include <utility>  // for move

#include "context_variables.hpp"

using algebra::RealFunction, algebra::RealFunctionWithPower, algebra::ComplexFunction;
using mpfr::real;
using std::move;

namespace qboot
{
	RealFunction<real> hBlock_shifted(const PrimaryOperator& op, const real& S, const real& P, uint32_t n_Max)
	{
		RealFunction<real> b(n_Max);
		if (op.is_divergent_hor())
		{
			real small(1ul, -(3 * mpfr::global_prec) / 8 + 15);
			b = hBlock_shifted(op.get_shifted(small), S, P, n_Max);
			b += hBlock_shifted(op.get_shifted(-small), S, P, n_Max);
			b /= 2;
			return b;
		}
		real sum;
		auto p = _get_rec_coeffs(op, S, P);
		b.at(0) = 1;
		for (uint32_t n = 1; n <= n_Max; ++n)
		{
			sum = 0;
			for (uint32_t i = 1; i < p.size(); ++i)
				if (i <= n) sum += p[i].eval(n) * b.at(n - i);
			b.at(n) = -sum / p[0].eval(n);
		}
		return b;
	}
	RealFunction<real> hBlock_powered(const real& exp, const PrimaryOperator& op, const real& S, const real& P,
	                                  const Context& context)
	{
		auto h_at_0 = hBlock_shifted(op, S, P, context.n_Max());
		h_at_0 *= mpfr::pow(4, exp);
		const auto& rho = context.rho();
		RealFunctionWithPower<real> f_at_0(move(h_at_0), exp);
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
	RealFunction<real> gBlock_real(const PrimaryOperator& op, const real& S, const real& P, const Context& context)
	{
		return hBlock_powered(op.delta(), op, S, P, context);
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
	RealFunction<real> h_asymptotic(const real& S, const Context& context)
	{
		// calculate (1 + sgn r) ^ {-epsilon - 1 + 2 sgn S} as a function of r - rho
		auto getFactor = [&context, &S](int sign) {
			return algebra::power_function(real(1) + sign * context.rho(), real(sign),
			                               (2 * sign) * S - 1 - context.epsilon(), context.lambda());
		};
		return context.rho_to_z().convert(mul(getFactor(1), getFactor(-1)));
	}
}  // namespace qboot
