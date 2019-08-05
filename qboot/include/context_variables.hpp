#ifndef CONTEXT_VARIABLES_H
#define CONTEXT_VARIABLES_H

#include <array>
#include <cstddef>
#include <vector>

#include "matrix.hpp"
#include "polynomial.hpp"
#include "real.hpp"

namespace qboot2
{
	typedef struct _context
	{
		long n_Max;
		mpfr_prec_t prec;
		mpfr_rnd_t rnd;
		int lambda;
		mpfr_t* rho_to_z_matrix;
		mpfr_t rho;
	} cb_context;

	/* basic constructor for cb_context */
	cb_context context_construct(long n_Max, mpfr_prec_t prec, int lambda);
	void clear_cb_context(cb_context& context);
}  // namespace qboot2

namespace qboot
{
	constexpr size_t _triangle_num(size_t n) noexcept { return n * (n + 1) / 2; }
	constexpr size_t get_dimG(size_t lambda) noexcept { return (lambda + 2) * (lambda + 2) / 4; }
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Matrix<Real> z_zbar_derivative_to_x_y_derivative_Matrix(size_t lambda)
	{
		auto dimG = get_dimG(lambda);
		algebra::Matrix result(dimG, dimG);
		auto to_index = [lambda](size_t a, size_t b) { return (lambda + 2 - a) * a + b; };
		std::vector<algebra::Polynomial<Real>> p, q;
		for (size_t i = 0; i <= lambda; i++) p.push_back(algebra::Polynomial<Real>::linear_power(Real(1), Real(1), i));
		for (size_t i = 0; i <= lambda; i++) q.push_back(algebra::Polynomial<Real>::linear_power(Real(1), Real(-1), i));
		for (size_t i = 0; i <= lambda / 2; i++)
			for (size_t j = i; j <= lambda - i; j++)
			{
				auto coeff = i == j ? p[i] * q[i] : p[i] * q[j] + p[j] * q[i];
				for (size_t p = (i + j) % 2; p <= coeff.degree(); p += 2)
					result.get(to_index((i + j - p) / 2, p), to_index(i, j - i)) = coeff[p];
			}
		return result;
	}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	struct cb_context
	{
		using type = Real;
		static constexpr mpfr_prec_t prec = Real::prec;
		static constexpr mpfr_rnd_t rnd = Real::rnd;
		const size_t n_Max, lambda;
		const Real rho;
		algebra::Matrix<Real> rho_to_z_matrix, zzbar_to_xy_marix;
		std::vector<std::array<size_t, 2>> index_list;
		std::vector<algebra::Polynomial<Real>> rho_to_delta;
		cb_context(size_t n_Max, size_t lambda)
		    : n_Max(n_Max),
		      lambda(lambda),
		      rho(3 - mpfr::sqrt(Real(8))),
		      rho_to_z_matrix(lambda + 1, lambda + 1),
		      zzbar_to_xy_marix(0, 0),
		      index_list(),
		      rho_to_delta(lambda + 1)
		{
			for (size_t i = 0; 2 * i <= lambda; i++)
				for (size_t j = 0; j <= lambda - 2 * i; j++) index_list.push_back(std::array{i, j});
			auto Lambda = lambda + 1;
			zzbar_to_xy_marix = z_zbar_derivative_to_x_y_derivative_Matrix<Real>(lambda);
			algebra::Vector<Real> tmps(Lambda);
			tmps[0] = -mpfr::sqrt(Real(8));
			for (size_t j = 1; j < Lambda; j++) tmps[j] = (tmps[j - 1] * (2 * Real(j) - 3)) / Real(j);
			tmps[1] -= 2;
			tmps[0] += 3;

			rho_to_z_matrix.get(0, 0) = 1;
			for (size_t j = 1; j < Lambda; j++)
			{
				Real tmp(1);
				for (size_t k = 0; k <= j; k++)
				{
					rho_to_z_matrix.get(1, j) += tmps[j - k] * tmp;
					tmp *= -2;
				}
			}
			for (size_t i = 2; i < Lambda; i++)
				for (size_t j = 1; j < Lambda; j++)
					for (size_t k = i - 1; k < Lambda - j; k++)
						rho_to_z_matrix.get(i, k + j) += rho_to_z_matrix.get(i - 1, k) * rho_to_z_matrix.get(1, j);

			rho_to_z_matrix.transpose();

			rho_to_delta[0] = {1};
			for (size_t i = 1; i <= lambda; i++)
			{
				auto tmp = Real(1) / (rho * Real(i));
				rho_to_delta[i] = rho_to_delta[i - 1] * algebra::Polynomial<Real>({Real(1 - i) * tmp, tmp});
			}
		}
		size_t dim_f() const noexcept { return _triangle_num((lambda + 1) / 2); }
		size_t dim_g() const noexcept { return _triangle_num((lambda + 2) / 2); }
		template <class T>
		Real operator()(const T& val)
		{
			return Real(val);
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_H
