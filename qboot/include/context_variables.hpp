#ifndef CONTEXT_VARIABLES_H
#define CONTEXT_VARIABLES_H

#include <cstddef>

#include "matrix.hpp"
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
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	struct cb_context
	{
		using type = Real;
		static constexpr mpfr_prec_t prec = Real::prec;
		static constexpr mpfr_rnd_t rnd = Real::rnd;
		const size_t n_Max, lambda;
		const Real rho;
		algebra::Matrix<Real> rho_to_z_matrix;
		cb_context(size_t n_Max, size_t lambda)
		    : n_Max(n_Max), lambda(lambda), rho(3 - mpfr::sqrt(Real(8))), rho_to_z_matrix(lambda + 1, lambda + 1)
		{
			auto Lambda = lambda + 1;
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
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_H
