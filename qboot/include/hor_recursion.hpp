#ifndef QBOOT_HOR_RECURSION_HPP_
#define QBOOT_HOR_RECURSION_HPP_

#include <cstdint>  // for uint32_t

#include "complex_function.hpp"   // for ComplexFunction
#include "context_variables.hpp"  // for Context, cb_context
#include "hor_formula.hpp"        // for _get_rec_coeffs
#include "primary_op.hpp"         // for PrimaryOperator
#include "real.hpp"               // for real, mpfr_prec_t, mpfr_rnd_t, mpfr_t, pow
#include "real_function.hpp"      // for RealFunction, RealFunctionWithPower, power_function

namespace qboot
{
	// a function f(rho) of rho at rho = 0 (not crossing symmetric point)
	// f(rho) = h_{Delta, spin}^{d12, d34}(z, z),
	// g_{Delta, spin}^{d12, d34}(z, z) = (4 rho) ^ {Delta} f(rho)
	// if p[0] may be 0, we use continuity of conformal block.
	template <class Real>
	algebra::RealFunction<Real> hBlock_shifted(const PrimaryOperator<Real>& op, const Real& S, const Real& P,
	                                           uint32_t n_Max)
	{
		algebra::RealFunction<Real> b(n_Max);
		if (op.is_divergent_hor())
		{
			Real small(1ul, -(3 * Real::prec) / 8 + 15);
			b = hBlock_shifted(op.get_shifted(small), S, P, n_Max);
			b += hBlock_shifted(op.get_shifted(-small), S, P, n_Max);
			b /= 2;
			return b;
		}
		Real sum;
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

	// a function of z - 1 / 2 expanded at z = 1 / 2,
	// (4 * rho) ^ {exp} * h_{Delta, spin}^{d12, d34}(z, z)
	// = (4 * rho) ^ {exp - Delta} * g_{Delta, spin}^{d12, d34}(z, z)
	template <class Real>
	algebra::RealFunction<Real> hBlock_powered(const Real& exp, const PrimaryOperator<Real>& op, const Real& S,
	                                           const Real& P, const Context<Real>& context)
	{
		auto h_at_0 = hBlock_shifted(op, S, P, context.n_Max());
		h_at_0 *= mpfr::pow(4, exp);
		const auto& rho = context.rho();
		algebra::RealFunctionWithPower<Real> f_at_0(std::move(h_at_0), exp);
		algebra::RealFunction<Real> f_of_rho(context.lambda());
		f_of_rho.at(0) = f_at_0.approximate(rho);
		Real tmp(1);
		for (uint32_t k = 1; k <= f_of_rho.lambda(); ++k)
		{
			tmp *= k;
			f_at_0.derivate();
			f_of_rho.at(k) = f_at_0.approximate(rho) / tmp;
		}
		return context.rho_to_z().convert(std::move(f_of_rho));
	}

	// a function of z - 1 / 2 expanded at z = 1 / 2, g_{Delta, spin}^{d12, d34}(z, z)
	template <class Real>
	algebra::RealFunction<Real> gBlock_real(const PrimaryOperator<Real>& op, const Real& S, const Real& P,
	                                        const Context<Real>& context)
	{
		return hBlock_powered(op.delta(), op, S, P, context);
	}

	template <class Real>
	algebra::ComplexFunction<Real> gBlock(const PrimaryOperator<Real>& op, const Real& S, const Real& P,
	                                      const Context<Real>& context)
	{
		return context.expand_off_diagonal(gBlock_real(op, S, P, context), op, S, P);
	}

	template <class Real>
	algebra::ComplexFunction<Real> gBlock(const PrimaryOperator<Real>& op, const Real& d1, const Real& d2,
	                                      const Real& d3, const Real& d4, const Context<Real>& context)
	{
		Real d12 = d1 - d2, d34 = d3 - d4;
		return gBlock(op, (d34 - d12) / 2, -d12 * d34 / 2, context);
	}

	// calculate \tilde{h}(r, 1) as a function of z - 1 / 2 eq (4.6) in arXiv:1406:4858
	// \tilde{h}(r, 1)
	//   = (1 - r ^ 2) ^ {-epsilon} (1 + r) ^ {-1 - d12 + d34} (1 - r) ^ {-1 + d12 - d34}
	//   = (1 + r) ^ {-epsilon - 1 + 2 S} (1 - r) ^ {-epsilon - 1 - 2 S}
	template <class Real>
	algebra::RealFunction<Real> h_asymptotic(const Real& S, const Context<Real>& context)
	{
		// calculate (1 + sgn r) ^ {-epsilon - 1 + 2 sgn S} as a function of r - rho
		auto getFactor = [&context, &S](int sign) {
			return algebra::power_function(Real(1) + sign * context.rho(), Real(sign),
			                               (2 * sign) * S - 1 - context.epsilon(), context.lambda());
		};
		return context.rho_to_z().convert(mul(getFactor(1), getFactor(-1)));
	}
}  // namespace qboot

#endif  // QBOOT_HOR_RECURSION_HPP_
