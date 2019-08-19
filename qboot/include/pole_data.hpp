#ifndef POLE_DATA_HPP_
#define POLE_DATA_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t, int32_t
#include <utility>  // for move
#include <vector>   // for vector

#include "complex_function.hpp"   // for ComplexFunction
#include "context_variables.hpp"  // for Context
#include "matrix.hpp"             // for Vector, Matrix
#include "polynomial.hpp"         // for Polynomial
#include "primary_op.hpp"         // for PrimaryOperator
#include "real.hpp"               // for real, factorial, pochhammer, pow
#include "real_function.hpp"      // for RealFunction

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class PoleData
	{
		uint32_t type, k, spin;
		Real epsilon, a, b;
		const Context<Real>& context;

	public:
		PoleData(uint32_t type, uint32_t k, uint32_t spin, const Context<Real>& context, const Real& a, const Real& b)
		    : type(type), k(k), spin(spin), epsilon(context.epsilon), a(a), b(b), context(context)
		{
			assert(1 <= type && type <= 3);
			assert(k >= 1);
			assert(type != 3 || k <= spin);
		}
		Real S() const { return a + b; }
		Real P() const { return 2 * a * b; }
		uint32_t descendant_level() const { return type == 2 ? 2 * k : k; }
		Real pole_position() const
		{
			switch (type)
			{
			case 1: return Real(-int32_t(spin + k - 1));
			case 2: return epsilon - (k - 1);
			default: return 2 * epsilon + (1 + spin - k);  // note: k <= spin -> 1 + spin - k >= 1
			}
		}
		Real residueDelta() const
		{
			switch (type)
			{
			case 1: return Real(1 - int32_t(spin));
			case 2: return epsilon + (k + 1);
			default: return 2 * epsilon + (spin + 1);
			}
		}
		uint32_t residueEll() const
		{
			switch (type)
			{
			case 1: return spin + k;
			case 2: return spin;
			default: return spin - k;
			}
		}
		Real coeff() const
		{
			Real p0 = mpfr::factorial<Real::prec, Real::rnd>(k);
			Real val = (k % 2 == 0 ? -1 : 1) * int32_t(k) / (p0 * p0);
			auto km1 = int32_t(k) - 1;
			switch (type)
			{
			case 1:
				val *= mpfr::pochhammer(spin + 2 * epsilon, k);
				val *= mpfr::pochhammer((2 * a - km1) / 2, k);
				val *= mpfr::pochhammer((2 * b - km1) / 2, k);
				return val /= mpfr::pochhammer(spin + epsilon, k);
			case 2:
			{
				auto x = int32_t(spin) - int32_t(k);
				val *= mpfr::pochhammer(epsilon - k, 2 * k);
				val /= mpfr::pochhammer(epsilon + x, 2 * k);
				val /= mpfr::pochhammer(epsilon + (x + 1), 2 * k);
				Real tmp = (epsilon + (x + 1)) / 2;
				val *= mpfr::pochhammer(tmp + a, k);
				val *= mpfr::pochhammer(tmp - a, k);
				val *= mpfr::pochhammer(tmp + b, k);
				return val *= mpfr::pochhammer(tmp - b, k);
			}
			default:
				val *= mpfr::pochhammer(Real(1 + spin - k), k);
				val *= mpfr::pochhammer((2 * a - km1) / 2, k);
				val *= mpfr::pochhammer((2 * b - km1) / 2, k);
				return val /= mpfr::pochhammer(epsilon + (1 + spin - k), k);
			}
		}
		PrimaryOperator<Real> get_op() const { return PrimaryOperator(residueDelta(), residueEll(), context); }
		RealFunction<Real> residue_of_h() const
		{
			return coeff() * context.h_times_rho_k(descendant_level(), get_op(), S(), P());
		}
	};
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class RationalApproxData
	{
		uint32_t cutoff, spin;
		Real epsilon, a, b, S{}, P{}, unitarity_bound{};
		std::vector<PoleData<Real>> poles{}, approx_poles{};
		algebra::Matrix<Real> approx_matrix_inv{0, 0};
		const Context<Real>& context;

	public:
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& d1, const Real& d2,
		                   const Real& d3, const Real& d4)
		    : RationalApproxData(cutoff, spin, context, d1 - d2, d3 - d4)
		{
		}
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& delta12,
		                   const Real& delta34)
		    : cutoff(cutoff), spin(spin), epsilon(context.epsilon), a(-delta12 / 2), b(delta34 / 2), context(context)
		{
			// type 1 or 3 PoleData vanishes when delta12 == 0 or delta34 == 0 and k is odd
			auto include_odd = delta12 != 0 && delta34 != 0;
			S = a + b;
			P = 2 * a * b;
			for (uint32_t x = 1; x <= cutoff; ++x)
				if (include_odd || x % 2 == 0) poles.emplace_back(1, x, spin, context, a, b);
			for (uint32_t x = 1; x <= cutoff / 2; ++x) poles.emplace_back(2, x, spin, context, a, b);
			for (uint32_t x = 1; x <= std::min(spin, cutoff / 2); ++x)
				if (include_odd || x % 2 == 0) poles.emplace_back(3, x, spin, context, a, b);
			auto N = uint32_t(poles.size());
			unitarity_bound = context.unitary_bound(spin) + Real(1L, -10);
			{
				algebra::Matrix<Real> approx_matrix{N, N};
				for (uint32_t i = 0; i < N; ++i)
				{
					auto vec = approx_row(poles[i].pole_position());
					for (uint32_t j = 0; j < N; ++j) approx_matrix.get(i, j) = vec.get(j);
				}
				approx_matrix_inv = approx_matrix.inverse();
			}
			for (uint32_t x = cutoff + 1; x <= 2 * cutoff; ++x)
				if (include_odd || x % 2 == 0) approx_poles.emplace_back(1, x, spin, context, a, b);
			for (uint32_t x = cutoff / 2 + 1; x <= 2 * (cutoff / 2); ++x)
				approx_poles.emplace_back(2, x, spin, context, a, b);
			for (uint32_t x = std::min(spin, cutoff / 2) + 1; x <= spin; ++x)
				if (include_odd || x % 2 == 0) approx_poles.emplace_back(3, x, spin, context, a, b);
		}
		// converts 1 / (x - p) into some vector
		algebra::Vector<Real> approx_row(const Real& p)
		{
			// regard 1 / (x - p) as an operator acting on a space of functions of x,
			// and take a basis ket(f_i) which spans the space spanned by this.poles.
			// we define a projector Pi = \sum_{i=0}^{N - 1} ket(f_i) bra(f_i),
			// then 1 / (x - p) are projected to \sum_{i=0}^{N - 1} braket(f_i | 1 / (x - p)) ket(f_i).
			// what we need is ket(1 / (x - p)) -> \sum_{i=0}^{N - 1} C_i(p) ket(1 / (x - poles[i]))
			// where C_i(p) = braket(1 / (x - poles[i]), 1 / (x - p)).
			// C_i(p) = \sum_{j = 0}^{N - 1} braket(1 / (x - poles[i]), f_j) braket(f_j, 1 / (x - p))
			auto N = uint32_t(poles.size());
			uint32_t j = 0;
			algebra::Vector<Real> v(N);
			// evaluate the (i - 1)-th taylor coeff at unitarity bound of 1 / (x - p)
			// braket(f_j | g) = (-1) ^ {i - 1} (der x) ^ {i - 1} g(x) / (i - 1)! | x = unitarity_bound
			for (int32_t i = 1; uint32_t(i) <= N / 2 + 1; ++i) v[j++] = mpfr::pow(unitarity_bound - p, -i);
			// evaluate the (i + 1)-th taylor coeff at infty of 1 / (x - p)
			// braket(f_j | g) = (der t) ^ {i + 1} g(1 / t) / (i + 1)! | t = 0
			for (uint32_t i = 0; i < (N + 1) / 2 - 1; ++i) v[j++] = mpfr::pow(p, i);
			return v;
		}
		RealFunction<algebra::Polynomial<Real>> approx_h()
		{
			auto N = uint32_t(poles.size());
			std::vector<algebra::Polynomial<Real>> pole_pol;
			for (uint32_t i = 0; i < N; ++i) pole_pol.push_back({-poles[i].pole_position(), Real(1)});
			// prod_left[i] = pole_pol[0] ... pole_pol[i - 1]
			// prod_right[i] = pole_pol[N - i] ... pole_pol[N - 1]
			std::vector<algebra::Polynomial<Real>> prod_left, prod_right;
			prod_left.emplace_back(0);
			prod_right.emplace_back(0);
			for (uint32_t i = 0; i < N; ++i) prod_left.push_back(prod_left[i] * pole_pol[i]);
			for (uint32_t i = 0; i < N; ++i) prod_right.push_back(prod_right[i] * pole_pol[N - i - 1]);

			auto res = RealFunction<algebra::Polynomial<Real>>(context.h_asymptotic_form(S));
			// multipluy (x - poles[i]) for all i
			res *= prod_left[N];

			// residues[i] = residue of 1 / (x - poles[i])
			std::vector<RealFunction<Real>> residues;
			for (uint32_t i = 0; i < N; ++i) residues.push_back(poles[i].residue_of_h());
			// for each ap_pole, approximate 1 / (x - ap_pole) as a sum of 1 / (x - poles[i])
			for (const auto& ap_pole : approx_poles)
			{
				auto v = approx_row(ap_pole.pole_position()) * approx_matrix_inv;
				auto r = ap_pole.residue_of_h();
				for (uint32_t i = 0; i < N; ++i) residues[i] += v[i] * r;
			}

			for (uint32_t i = 0; i < N; ++i)
			{
				auto f = RealFunction<algebra::Polynomial<Real>>(residues[i]);
				// multiply (x - poles[j]) for all j but i
				f *= prod_left[i] * prod_right[N - i - 1];
				res += f;
			}
			return res;
		}
		// gives reduced gBlock
		// Let returned value to be f(z), g = (4 rho) ^ {Delta} f(z) / \prod_{i = 0}^{N - 1} (Delta - poles[i])
		ComplexFunction<algebra::Polynomial<Real>> approx_g()
		{
			return context.expand_off_diagonal(context.rho_to_delta * approx_h(), spin, S, P);
		}
		ComplexFunction<Real> my_approx(const Real& delta)
		{
			auto N = uint32_t(poles.size());

			auto sum = context.h_asymptotic_form(S);
			// auto res = RealFunction<Real>(context.h_asymptotic_form(S));
			// std::cout << context.expand_off_diagonal(context.rho_to_delta * res, spin, S, P) << std::endl;

			// residues[i] = residue of 1 / (x - poles[i])
			std::vector<RealFunction<Real>> residues;
			for (uint32_t i = 0; i < N; ++i) residues.push_back(poles[i].residue_of_h());
			// for each ap_pole, approximate 1 / (x - ap_pole) as a sum of 1 / (x - poles[i])
			for (const auto& ap_pole : approx_poles)
			{
				auto v = approx_row(ap_pole.pole_position()) * approx_matrix_inv;
				auto r = ap_pole.residue_of_h();
				for (uint32_t i = 0; i < N; ++i) residues[i] += v[i] * r;
			}

			for (uint32_t i = 0; i < N; ++i)
			{
				auto f = RealFunction<algebra::Polynomial<Real>>(residues[i]);
				RealFunction<Real> g(f.lambda());
				for (uint32_t k = 0; k <= f.lambda(); ++k)
					g.get(k) = f.get(k).eval(delta) / (delta - poles[i].pole_position());
				sum += g;
			}
			RealFunction<Real> rtd(context.rho_to_delta.lambda());
			for (uint32_t k = 0; k <= rtd.lambda(); ++k) rtd.get(k) = context.rho_to_delta.get(k).eval(delta);
			rtd *= mpfr::pow(4 * context.rho, delta);
			return context.expand_off_diagonal(rtd * sum, context.get_primary(delta, spin), S, P);
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_HPP_
