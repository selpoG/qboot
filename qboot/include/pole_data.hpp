#ifndef POLE_DATA_HPP_
#define POLE_DATA_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t, int32_t
#include <vector>   // for vector

#include "context_variables.hpp"  // for Context
#include "matrix.hpp"             // for Vector, Matrix
#include "polynomial.hpp"         // for Polynomial
#include "real.hpp"               // for real, factorial, pochhammer. pow

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class PoleData
	{
		uint32_t type, k;
		Real ell, epsilon, a, b;

	public:
		PoleData(uint32_t type, uint32_t k, uint32_t spin, const Real& eps, const Real& a, const Real& b)
		    : type(type), k(k), ell(Real(spin)), epsilon(eps), a(a), b(b)
		{
			assert(1 <= type && type <= 3);
			assert(type != 3 || k <= spin);
		}
		Real S() const { return a + b; }
		Real P() const { return 2 * a * b; }
		uint32_t descendant_level() const { return type == 2 ? 2 * k : k; }
		Real pole_position() const
		{
			switch (type)
			{
			case 1: return 1 - ell - k;
			case 2: return 1 + epsilon - k;
			default: return 1 + ell + 2 * epsilon - k;
			}
		}
		Real residueDelta() const
		{
			switch (type)
			{
			case 1: return 1 - ell;
			case 2: return 1 + epsilon + k;
			default: return 1 + ell + 2 * epsilon;
			}
		}
		Real residueEll() const
		{
			switch (type)
			{
			case 1: return ell + k;
			case 2: return ell;
			default: return ell - k;
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
				val *= mpfr::pochhammer(ell + 2 * epsilon, k);
				val *= mpfr::pochhammer((2 * a - km1) / 2, k);
				val *= mpfr::pochhammer((2 * b - km1) / 2, k);
				return val /= mpfr::pochhammer(ell + epsilon, k);
			case 2:
			{
				val *= mpfr::pochhammer(epsilon - k, 2 * k);
				val /= mpfr::pochhammer(ell + epsilon - k, 2 * k);
				val /= mpfr::pochhammer(ell + epsilon - km1, 2 * k);
				Real tmp = (ell + epsilon - km1) / 2;
				val *= mpfr::pochhammer(tmp + a, k);
				val *= mpfr::pochhammer(tmp - a, k);
				val *= mpfr::pochhammer(tmp + b, k);
				return val *= mpfr::pochhammer(tmp - b, k);
			}
			default:
				val *= mpfr::pochhammer(ell + 1 - k, k);
				val *= mpfr::pochhammer((2 * a - km1) / 2, k);
				val *= mpfr::pochhammer((2 * b - km1) / 2, k);
				return val /= mpfr::pochhammer(1 + ell + epsilon - k, k);
			}
		}
		algebra::Vector<Real> residue_of_h(const Context<Real>& c) const
		{
			return coeff() * c.h_times_rho_k(descendant_level(), residueEll(), residueDelta(), S(), P());
		}
	};
	template <class Real>
	class RationalApproxData
	{
		uint32_t cutoff;
		Real epsilon, ell, a, b, S{}, P{}, unitarity_bound{};
		std::vector<PoleData<Real>> poles{}, approx_poles{};
		algebra::Matrix<Real> approx_matrix_inv{0, 0};

	public:
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Real& eps, const Real& delta12, const Real& delta34)
		    : cutoff(cutoff), epsilon(eps), ell(Real(spin)), a(-delta12 / 2), b(delta34 / 2)
		{
			bool is_correlator_multiple = delta12 != 0 || delta34 != 0;
			S = a + b;
			P = 2 * a * b;
			for (uint32_t x = 1; x <= cutoff; x++)
				if (is_correlator_multiple || x % 2 == 0) poles.emplace_back(1, x, spin, epsilon, a, b);
			for (uint32_t x = 1; x <= cutoff / 2; x++) poles.emplace_back(2, x, spin, epsilon, a, b);
			for (uint32_t x = 1; x <= std::min(spin, cutoff / 2); x++)
				if (is_correlator_multiple || x % 2 == 0) poles.emplace_back(3, x, spin, epsilon, a, b);
			auto N = uint32_t(poles.size());
			unitarity_bound = spin == 0 ? epsilon : ell + 2 * epsilon;
			unitarity_bound += Real(1L, -10);
			algebra::Matrix<Real> approx_matrix{N, N};
			for (uint32_t i = 0; i < poles.size(); i++)
			{
				auto vec = approx_column(poles[i]);
				for (uint32_t j = 0; j < poles.size(); j++) approx_matrix.get(i, j) = vec.get(j);
			}
			approx_matrix.transpose();
			approx_matrix_inv = approx_matrix.inverse();
			for (uint32_t x = cutoff + 1; x <= 2 * cutoff; x++)
				if (is_correlator_multiple || x % 2 == 0) approx_poles.emplace_back(1, x, spin, epsilon, a, b);
			for (uint32_t x = cutoff / 2 + 1; x <= 2 * (cutoff / 2); x++)
				approx_poles.emplace_back(2, x, spin, epsilon, a, b);
			for (uint32_t x = std::min(spin, cutoff / 2) + 1; x <= spin; x++)
				if (is_correlator_multiple || x % 2 == 0) approx_poles.emplace_back(3, x, spin, epsilon, a, b);
		}
		algebra::Vector<Real> approx_column(const PoleData<Real>& p)
		{
			auto dim_approx_base = uint32_t(poles.size());
			uint32_t j = 0;
			algebra::Vector<Real> res(dim_approx_base);
			for (uint32_t i = 1; i <= dim_approx_base / 2 + 1; i++)
				res[j++] = 1 / mpfr::pow(unitarity_bound - p.pole_position(), i);
			for (uint32_t i = 0; i < (dim_approx_base + 1) / 2 - 1; i++) res[j++] = mpfr::pow(p.pole_position(), i);
			return res;
		}
		algebra::Vector<algebra::Polynomial<Real>> approx_h(const Context<Real>& c)
		{
			auto N = uint32_t(poles.size());
			std::vector<algebra::Polynomial<Real>> pole_pol;
			for (uint32_t i = 0; i < N; i++)
				pole_pol.push_back(algebra::Polynomial<Real>{-poles[i].pole_position(), Real(1)});
			auto res = algebra::Vector<algebra::Polynomial<Real>>(c.h_asymptotic_form(S));
			for (uint32_t i = 0; i < N; i++) res *= pole_pol[i];
			// prod_left[i] = pole_pol[0] * ... * pole_pol[i - 1]
			// prod_right[i] = pole_pol[N - i] * ... * pole_pol[N - 1]
			std::vector<algebra::Polynomial<Real>> prod_left, prod_right;
			prod_left.emplace_back(0);
			prod_right.emplace_back(0);
			for (uint32_t i = 0; i < N; i++) prod_left.push_back(prod_left[i] * pole_pol[i]);
			for (uint32_t i = 0; i < N; i++) prod_right.push_back(prod_right[i] * pole_pol[N - i - 1]);
			for (uint32_t i = 0; i < N; i++) res += poles[i].residue_of_h(c) * prod_left[i] * prod_right[N - i - 1];
			if (!approx_poles.empty())
			{
				auto col = uint32_t(approx_poles.size());
				algebra::Matrix<Real> approx_target(col, N);
				for (uint32_t i = 0; i < col; i++)
				{
					auto v = approx_column(approx_poles[i]);
					for (uint32_t j = 0; j < N; j++) approx_target.get(i, j) = v.get(j);
				}
				approx_target.transpose();
				algebra::Vector<algebra::Polynomial<Real>> polys(N);
				for (uint32_t i = 0; i < N; i++) polys[i] = prod_left[i] * prod_right[N - i - 1];
				auto approx_polys = (polys * approx_matrix_inv) * approx_target;
				for (uint32_t i = 0; i < N; i++) res += approx_poles[i].residue_of_h(c) * approx_polys[i];
			}
			return res;
		}
		algebra::Vector<algebra::Polynomial<Real>> approx_g(const Context<Real>& c)
		{
			const auto& v1 = c.rho_to_delta;
			auto v2 = approx_h(c);
			algebra::Vector<algebra::Polynomial<Real>> v(c.lambda + 1);
			for (uint32_t i = 0; i <= c.lambda; i++)
			{
				v[i] = {};
				for (uint32_t j = 0; j <= i; j++) v[i] += v1[j] * v2[i - j];
			}
			return c.c2_expand(std::move(v), ell, S, P);
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_HPP_
