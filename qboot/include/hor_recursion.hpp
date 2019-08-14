#ifndef HOR_RECURSION_HPP_
#define HOR_RECURSION_HPP_

#include <cstdint>  // for uint32_t, int32_t
#include <memory>   // for unique_ptr

#include "context_variables.hpp"  // for Context, cb_context
#include "hor_formula.hpp"        // for _get_rec_coeffs, _evaluate_at_n
#include "matrix.hpp"             // for Vector
#include "primary_op.hpp"         // for PrimaryOperator
#include "real.hpp"               // for real, mpfr_prec_t, mpfr_rnd_t, mpfr_t, is_integer, factorial, pow

namespace qboot2
{
	std::unique_ptr<mpfr_t[]> recursionNonZeroVector(uint32_t nMax, const mpfr_t& epsilon, const mpfr_t& ell,
	                                                 const mpfr_t& Delta, const mpfr_t& S, const mpfr_t& P,
	                                                 mpfr_prec_t prec, mpfr_rnd_t rnd);
	std::unique_ptr<mpfr_t[]> recursionSpinZeroVector(uint32_t nMax, const mpfr_t& epsilon, mpfr_t& Delta,
	                                                  const mpfr_t& S, const mpfr_t& P, mpfr_prec_t prec,
	                                                  mpfr_rnd_t rnd);
	std::unique_ptr<mpfr_t[]> real_axis_result(const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta, const mpfr_t& S,
	                                           const mpfr_t& P, const cb_context& context);
	uint32_t indexOfConformalBlock(const cb_context& context, int32_t n, int32_t m);
	void element_helper(const cb_context& context, const std::unique_ptr<mpfr_t[]>& array, mpfr_t& r, int32_t m,
	                    int32_t n);
	std::unique_ptr<mpfr_t[]> casimirExpander(mpfr_t* realAxisResult, const mpfr_t& epsilon, const mpfr_t& ell,
	                                          const mpfr_t& Delta, const mpfr_t& S, const mpfr_t& P,
	                                          const cb_context& context);
	std::unique_ptr<mpfr_t[]> gBlock_full(const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta, const mpfr_t& S,
	                                      const mpfr_t& P, const cb_context& context);
	std::unique_ptr<mpfr_t[]> hBlock_times_rho_n(uint32_t n, const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta,
	                                             const mpfr_t& S, const mpfr_t& P, const cb_context& context);
	std::unique_ptr<mpfr_t[]> h_asymptotic(const mpfr_t& epsilon, const mpfr_t& S, const cb_context& context);
}  // namespace qboot2

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> power_series_in_rho(const PrimaryOperator<Real>& op, const Real& S, const Real& P);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> hBlock_powered(const Real& exp, const PrimaryOperator<Real>& op, const Real& S,
	                                     const Real& P);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> gBlock_full(const PrimaryOperator<Real>& op, const Real& S, const Real& P);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> h_asymptotic(const Real& S, const Context<Real>& context);

	template <class Real>
	PrimaryOperator<Real> _shift_op(const PrimaryOperator<Real>& op, const Real& small)
	{
		return op.context().get_primary(op.delta() == 0 ? small : op.delta() * (1 + small), op.spin());
	}

	// calculate b[n].
	// g_{\Delta, spin}^{d12, d34}(z, z) ~ (4 * \rho) ^ \Delta * \sum_{n = 0}^{nMax} b[n] * \rho ^ n
	// if p[0] may be 0, we use continuity of conformal block.
	template <class Real>
	algebra::Vector<Real> power_series_in_rho(const PrimaryOperator<Real>& op, const Real& S, const Real& P)
	{
		algebra::Vector<Real> b(op.context().n_Max + 1);
		if (op.is_divergent_hor())
		{
			Real small(1ul, -(3 * Real::prec) / 8 + 15);
			b = power_series_in_rho(_shift_op(op, small), S, P);
			b += power_series_in_rho(_shift_op(op, -small), S, P);
			b /= Real(2);
			return b;
		}
		Real sum;
		auto p = _get_rec_coeffs(op, S, P);
		b[0] = 1;
		for (uint32_t n = 1; n < b.size(); ++n)
		{
			sum = 0;
			for (uint32_t i = 1; i < p.size(); ++i)
				if (i <= n) sum += p[i].eval(n) * b[n - i];
			b[n] = -sum / p[0].eval(n);
		}
		return b;
	}

	// let f(z) = g_{\Delta, spin}^{d12, d34}(z, z), which is a function of a on real axis.
	// hBlock_powered gives (f(z), f'(z), ..., f^{(lambda)}(z)).
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> hBlock_powered(const PrimaryOperator<Real>& op, const Real& S, const Real& P)
	{
		return hBlock_powered(op.delta(), op, S, P);
	}

	// let f(z) = (4 * rho) ^ {exp} * h_{\Delta, spin}^{d12, d34}(z, z)
	//   = (4 * rho) ^ {exp - \Delta} * g_{\Delta, spin}^{d12, d34}(z, z), which is a function of a on real axis.
	// hBlock_powered gives (f(z), f'(z), ..., f^{(lambda)}(z)).
	template <class Real>
	algebra::Vector<Real> hBlock_powered(const Real& exp, const PrimaryOperator<Real>& op, const Real& S, const Real& P)
	{
		auto b = power_series_in_rho(op, S, P);
		// f(z) ~ (4 * \rho) ^ {exp} * \sum_{n = 0}^{nMax} b[n] * \rho ^ n
		// k-th derivative of f(z) in \rho is
		//   (4 * \rho) ^ {exp} *
		//     \sum_{n = 0}^{nMax} (exp + n) * (exp + n - 1) * ... * (exp + n - k + 1) b[n] * \rho ^ {n - k}
		algebra::Vector<Real> result_in_rho(op.context().lambda + 1);
		Real tmp;
		const auto& rho = op.context().rho;
		result_in_rho[0] = 0;
		tmp = mpfr::pow(4 * rho, exp);
		for (uint32_t n = 0; n < b.size(); ++n)
		{
			b[n] *= tmp;
			result_in_rho[0] += b[n];
			tmp *= rho;
		}
		// b[n] -> (4 * \rho) ^ {exp} b[n] * \rho ^ n
		tmp = 1;
		for (int32_t k = 1; uint32_t(k) < result_in_rho.size(); ++k)
		{
			result_in_rho[uint32_t(k)] = 0;
			for (int32_t n = 0; uint32_t(n) < b.size(); ++n)
			{
				b[uint32_t(n)] *= exp + (n - k + 1);
				result_in_rho[uint32_t(k)] += b[uint32_t(n)];
			}
			tmp *= k;
			tmp *= rho;
			result_in_rho[uint32_t(k)] /= tmp;
		}
		// convert derivatives in rho to derivatives in z
		return op.context().rho_to_z_matrix * result_in_rho;
	}

	template <class Real>
	algebra::Vector<Real> gBlock_full(const PrimaryOperator<Real>& op, const Real& S, const Real& P)
	{
		return op.context().expand_off_diagonal(hBlock_powered(op, S, P), op, S, P);
	}

	template <class Real>
	algebra::Vector<Real> h_asymptotic(const Real& S, const Context<Real>& context)
	{
		auto getFactor = [&context, &S](int sign) {
			Real temp1, temp2;
			temp1 = (2 * sign) * S - 1 - context.epsilon;
			algebra::Vector<Real> factor(context.lambda + 1);
			temp2 = 1 + sign * context.rho;
			factor[0] = mpfr::pow(temp2, temp1);
			temp2 = sign / temp2;
			for (uint32_t j = 1; j <= context.lambda; ++j)
			{
				factor[j] = factor[j - 1] * (temp1 - (j - 1));
				factor[j] *= temp2;
				factor[j] /= j;
			}
			return factor;
		};
		auto firstFactor = getFactor(1);
		auto secondFactor = getFactor(-1);
		auto result_in_rho = firstFactor.convolve(secondFactor);
		return context.rho_to_z_matrix * result_in_rho;
	}
}  // namespace qboot

#endif  // HOR_RECURSION_HPP_
