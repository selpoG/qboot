#ifndef CONTEXT_VARIABLES_HPP_
#define CONTEXT_VARIABLES_HPP_

#include <array>    // for array
#include <cstdint>  // for uint32_t, int32_t
#include <memory>   // for unique_ptr
#include <sstream>  // for ostringstream
#include <string>   // for string
#include <vector>   // for vector

#include "matrix.hpp"      // for Vector, Matrix
#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for mpfr_prec_t, mpfr_rnd_t, mpfr_t, real, sqrt, pow

namespace qboot2
{
	typedef struct _context
	{
		uint32_t n_Max{};
		mpfr_prec_t prec{};
		mpfr_rnd_t rnd{};
		uint32_t lambda{};
		std::unique_ptr<mpfr_t[]> rho_to_z_matrix{};
		mpfr_t rho;
	} cb_context;

	/* basic constructor for cb_context */
	cb_context context_construct(uint32_t n_Max, mpfr_prec_t prec, uint32_t lambda);
	void clear_cb_context(cb_context& context);
}  // namespace qboot2

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class RationalApproxData;
	constexpr uint32_t _triangle_num(uint32_t n) noexcept { return n * (n + 1) / 2; }
	constexpr uint32_t get_dimG(uint32_t lambda) noexcept { return (lambda + 2) * (lambda + 2) / 4; }
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Matrix<Real> z_zbar_derivative_to_x_y_derivative_Matrix(uint32_t lambda)
	{
		auto dimG = get_dimG(lambda);
		algebra::Matrix result(dimG, dimG);
		auto to_index = [lambda](uint32_t a, uint32_t b) { return (lambda + 2 - a) * a + b; };
		std::vector<algebra::Polynomial<Real>> p, q;
		for (uint32_t i = 0; i <= lambda; i++)
			p.push_back(algebra::Polynomial<Real>::linear_power(Real(1), Real(1), i));
		for (uint32_t i = 0; i <= lambda; i++)
			q.push_back(algebra::Polynomial<Real>::linear_power(Real(1), Real(-1), i));
		for (uint32_t i = 0; i <= lambda / 2; i++)
			for (uint32_t j = i; j <= lambda - i; j++)
			{
				auto coeff = i == j ? p[i] * q[i] : p[i] * q[j] + p[j] * q[i];
				if (coeff.iszero()) continue;
				auto d = uint32_t(coeff.degree());
				for (uint32_t k = (i + j) % 2; k <= d; k += 2)
					result.get(to_index((i + j - k) / 2, k), to_index(i, j - i)) = coeff[k];
			}
		return result;
	}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Context
	{
	public:
		using type = Real;
		static constexpr mpfr_prec_t prec = Real::prec;
		static constexpr mpfr_rnd_t rnd = Real::rnd;
		const uint32_t n_Max, lambda, dimension;
		const Real epsilon, rho;
		algebra::Matrix<Real> rho_to_z_matrix, zzbar_to_xy_marix;
		std::vector<std::array<uint32_t, 2>> index_list;
		algebra::Vector<algebra::Polynomial<Real>> rho_to_delta;
		std::string str() const
		{
			std::ostringstream os;
			os << "Context(nMax=" << n_Max << ", lambda=" << lambda << ", dim=" << dimension << ")";
			return os.str();
		}
		Context(uint32_t n_Max, uint32_t lambda, uint32_t dim)
		    : n_Max(n_Max),
		      lambda(lambda),
		      dimension(dim),
		      epsilon(Real(dim - 2) / 2),
		      rho(3 - mpfr::sqrt(Real(8))),
		      rho_to_z_matrix(lambda + 1, lambda + 1),
		      zzbar_to_xy_marix(0, 0),
		      index_list(),
		      rho_to_delta(lambda + 1)
		{
			assert(dim >= 3 && dim % 2 == 1);
			for (uint32_t i = 0; 2 * i <= lambda; i++)
				for (uint32_t j = 0; j <= lambda - 2 * i; j++) index_list.push_back(std::array{i, j});
			auto Lambda = lambda + 1;
			zzbar_to_xy_marix = z_zbar_derivative_to_x_y_derivative_Matrix<Real>(lambda);
			algebra::Vector<Real> tmps(Lambda);
			tmps[0] = -mpfr::sqrt(Real(8));
			for (uint32_t j = 1; j < Lambda; j++) tmps[j] = (tmps[j - 1] * (2 * Real(j) - 3)) / Real(j);
			tmps[1] -= 2;
			tmps[0] += 3;

			Real tmp;
			rho_to_z_matrix.get(0, 0) = 1;
			for (uint32_t j = 1; j < Lambda; j++)
			{
				tmp = 1;
				for (uint32_t k = 0; k <= j; k++)
				{
					rho_to_z_matrix.get(1, j) += tmps[j - k] * tmp;
					tmp *= -2;
				}
			}
			for (uint32_t i = 2; i < Lambda; i++)
				for (uint32_t j = 1; j < Lambda; j++)
					for (uint32_t k = i - 1; k < Lambda - j; k++)
						rho_to_z_matrix.get(i, k + j) += rho_to_z_matrix.get(i - 1, k) * rho_to_z_matrix.get(1, j);

			rho_to_z_matrix.transpose();

			rho_to_delta[0] = {Real(1)};
			for (uint32_t i = 1; i <= lambda; i++)
			{
				tmp = 1 / (rho * i);
				rho_to_delta[i] = rho_to_delta[i - 1] * algebra::Polynomial<Real>({Real(1 - int32_t(i)) * tmp, tmp});
			}
			rho_to_delta = rho_to_z_matrix * rho_to_delta;
		}
		uint32_t dim_f() const noexcept { return _triangle_num((lambda + 1) / 2); }
		uint32_t dim_h() const noexcept { return _triangle_num((lambda + 2) / 2); }
		template <class T>
		Real operator()(const T& val) const
		{
			return Real(val);
		}
		// compute the table of derivative of v = (z z_bar) ^ d in the x_y basis
		algebra::Vector<Real> v_to_d(const Real& d) const
		{
			algebra::Vector<Real> table(lambda + 1);
			table[0] = Real(1);
			for (uint32_t i = 1; i <= lambda; i++) table[i] = table[i - 1] * (-2 * (d - Real(i - 1))) / Real(i);
			uint32_t dimG = get_dimG(lambda);
			algebra::Vector<Real> res(dimG);
			uint32_t pos = 0;
			for (uint32_t i = 0; i <= lambda; i++)
				for (uint32_t j = i; j <= lambda - i; j++) res[pos++] = table[i] * table[j];
			return zzbar_to_xy_marix * res;
		}
		// parity = 0 or 1
		algebra::Matrix<Real> F_matrix_impl(const Real& d, uint32_t parity) const
		{
			uint32_t row = parity > 0 ? dim_f() : dim_h(), col = get_dimG(lambda);
			algebra::Matrix<Real> ans(row, col);
			auto v = v_to_d(d);
			auto k = mpfr::pow(Real(0.25), d);
			uint32_t r = 0;
			for (const auto& i : index_list)
			{
				if (i[1] % 2 != parity) continue;
				uint32_t c = 0;
				for (const auto& m : index_list)
				{
					if (i[0] >= m[0] && i[1] >= m[1])
					{
						uint32_t d0 = i[0] - m[0], d1 = i[1] - m[1];
						ans.get(r, c) = k * v[(lambda + 2 - d0) * d0 + d1];
					}
					else
						ans.get(r, c) = {};
					++c;
				}
				++r;
			}
			return ans;
		}
		algebra::Matrix<Real> F_minus_matrix(const Real& d) const { return F_matrix_impl(d, 1); }
		algebra::Matrix<Real> F_plus_matrix(const Real& d) const { return F_matrix_impl(d, 0); }
		algebra::Vector<Real> h_times_rho_k(uint32_t k, const Real& ell, const Real& Delta, const Real& S,
		                                    const Real& P) const
		{
			return hBlock_times_rho_n(k, epsilon, ell, Delta, S, P, *this);
		}
		uint32_t aligned_index(uint32_t y, uint32_t x) const noexcept { return (lambda + 2 - y) * y + x; }
		algebra::Vector<Real> h_asymptotic_form(const Real& S) const { return h_asymptotic(epsilon, S, *this); }
		algebra::Vector<algebra::Polynomial<Real>> c2_expand(algebra::Vector<algebra::Polynomial<Real>>&& array_real,
		                                                     const Real& ell, const Real& S, const Real& P) const
		{
			algebra::Polynomial<Real> local_c2{ell * (ell + 2 * epsilon), -2 - 2 * epsilon, Real(1)};
			algebra::Vector<algebra::Polynomial<Real>> ans(get_dimG(lambda));
			for (uint32_t i = 0; i <= lambda; i++) ans[i] = std::move(array_real[i]);
			for (uint32_t i = 1; i <= lambda / 2; i++)
			{
				Real ri = Real(i);
				for (uint32_t j = 0; j + 2 * i <= lambda; j++)
				{
					Real rj = Real(j);
					algebra::Polynomial<Real> val{};
					Real common_factor = 2 * epsilon + ri * 2 - 1;
					if (j >= 3) val += ans[aligned_index(i, j - 3)] * Real(16);
					if (j >= 2) val += ans[aligned_index(i, j - 2)] * Real(8);
					if (j >= 1) val -= ans[aligned_index(i, j - 1)] * Real(4);
					val *= common_factor;
					val += (-(rj + 1) * (rj + 2) / ri) * ans[aligned_index(i - 1, j + 2)];
					val += (2 * (rj + 1) * (2 * S + 2 * epsilon - 4 * ri - rj + 6) / ri) *
					       ans[aligned_index(i - 1, j + 1)];
					val += (Real(4) *
					        (2 * P + 8 * S * ri + 4 * S * rj - 8 * S + Real(2) * local_c2 + 4 * epsilon * ri +
					         4 * epsilon * rj - 4 * epsilon + 4 * ri * ri + 8 * ri * rj - 2 * ri + rj * rj - 5 * rj -
					         Real(2)) /
					        ri) *
					       ans[aligned_index(i - 1, j)];
					if (j >= 1)
						val += (8 *
						        (2 * P + 8 * S * ri + 2 * S * rj - 10 * S - 4 * epsilon * ri + 2 * epsilon * rj +
						         2 * epsilon + 12 * ri * ri + 12 * ri * rj - 34 * ri + rj * rj - 13 * rj + 22) /
						        ri) *
						       ans[aligned_index(i - 1, j - 1)];
					if (i >= 2)
					{
						val += (4 * ((rj + 1) * (rj + 2)) / ri) * ans[aligned_index(i - 2, j + 2)];
						val += (8 * (rj + 1) * (2 * S - 2 * epsilon + 4 * ri + 3 * rj - 6) / ri) *
						       ans[aligned_index(i - 2, j + 1)];
					}
					ans[aligned_index(i, j)] = val / (2 * common_factor);
				}
			}
			return ans;
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_HPP_
