#ifndef HOR_RECURSION_H
#define HOR_RECURSION_H

#include "context_variables.hpp"
#include "hor_formula.hpp"
#include "matrix.hpp"

namespace qboot2
{
	mpfr_t* recursionNonZeroVector(unsigned long nMax, const mpfr_t& epsilon, const mpfr_t& ell, const mpfr_t& Delta,
	                               const mpfr_t& S, const mpfr_t& P, mpfr_prec_t prec, mpfr_rnd_t rnd);
	mpfr_t* recursionSpinZeroVector(unsigned long nMax, const mpfr_t& epsilon, mpfr_t& Delta, const mpfr_t& S,
	                                const mpfr_t& P, mpfr_prec_t prec, mpfr_rnd_t rnd);
	mpfr_t* real_axis_result(const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta, const mpfr_t& S, const mpfr_t& P,
	                         const cb_context& context);
	long indexOfConformalBlock(const cb_context& context, int n, int m);
	void element_helper(const cb_context& context, mpfr_t* array, mpfr_t& r, int m, int n);
	mpfr_t* casimirExpander(mpfr_t* realAxisResult, const mpfr_t& epsilon, const mpfr_t& ell, const mpfr_t& Delta,
	                        const mpfr_t& S, const mpfr_t& P, const cb_context& context);
	mpfr_t* gBlock_full(const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta, const mpfr_t& S, const mpfr_t& P,
	                    const cb_context& context);
	mpfr_t* hBlock_times_rho_n(unsigned long n, const mpfr_t& epsilon, const mpfr_t& ell, mpfr_t& Delta,
	                           const mpfr_t& S, const mpfr_t& P, const cb_context& context);
	mpfr_t* h_asymptotic(const mpfr_t& epsilon, const mpfr_t& S, const cb_context& context);
}  // namespace qboot2

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> _recursion_vector(size_t nMax, const Real& epsilon, const Real& ell, const Real& Delta,
	                                        const Real& S, const Real& P);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> _real_axis_result(const Real& epsilon, const Real& ell, const Real& Delta, const Real& S,
	                                        const Real& P, const cb_context<Real>& context);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> _casimir_expander(const algebra::Vector<Real>& realAxisResult, const Real& epsilon,
	                                        const Real& ell, const Real& Delta, const Real& S, const Real& P,
	                                        size_t lambda);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> gBlock_full(const Real& epsilon, const Real& ell, const Real& Delta, const Real& S,
	                                  const Real& P, const cb_context<Real>& context);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> hBlock_times_rho_n(size_t n, const Real& epsilon, const Real& ell, const Real& Delta,
	                                         const Real& S, const Real& P, const cb_context<Real>& context);
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> h_asymptotic(const Real& epsilon, const Real& S, const cb_context<Real>& context);

	template <class Real>
	algebra::Vector<Real> _recursion_vector(size_t nMax, const Real& epsilon, const Real& ell, const Real& Delta,
	                                        const Real& S, const Real& P)
	{
		size_t minNMax = ell > 0 ? 8 : 6;
		assert(nMax >= minNMax);
		auto order = nMax + 1;
		algebra::Vector<Real> result(order);
		// if delta >= 2 * epsilon + ell and i >= 1 and epsilon >= 0 and ell >= 1,
		// av[0] == 0 iff. i == 1 and delta == 2 * epsilon + ell
		// if delta >= epsilon and i >= 1 and epsilon >= 0 and ell == 0
		// av[0] == 0 iff. (i == 1 and delta == ell) or (i == 2 and delta == ell + 1 / 2)
		auto av_divergent = (ell > 0 && Delta == 2 * epsilon + ell) ||
		                    (ell == 0 && (Delta == epsilon || Delta == epsilon + Real("0.5")));
		if (av_divergent)
		{
			Real small(1ul, -(3 * Real::prec) / 8 + 15);
			result = _recursion_vector(nMax, epsilon, ell, Delta * (1 + small), S, P);
			result += _recursion_vector(nMax, epsilon, ell, Delta * (1 - small), S, P);
			result /= Real(2);
			return result;
		}
		Real sum;
		auto recCoeffs = _get_rec_coeffs(epsilon, ell, Delta, S, P);
		result[0] = 1;
		for (size_t i = 1; i <= nMax; i++)
		{
			auto av = _evaluate_at_n(recCoeffs, i);
			sum = 0;
			auto jmax = i < minNMax ? i : minNMax - 1;
			for (size_t j = 1; j <= jmax; j++) sum += av[j] * result[i - j];
			result[i] = -sum / av[0];
		}
		return result;
	}

	template <class Real>
	algebra::Vector<Real> _real_axis_result(const Real& epsilon, const Real& ell, const Real& Delta, const Real& S,
	                                        const Real& P, const cb_context<Real>& context)
	{
		auto hBlock = _recursion_vector(context.n_Max, epsilon, ell, Delta, S, P);
		algebra::Vector<Real> result_in_rho(context.lambda + 1);
		Real tmp;
		result_in_rho[0] = 0;
		tmp = mpfr::pow(4 * context.rho, Delta);
		for (size_t j = 0; j <= context.n_Max; j++)
		{
			hBlock[j] *= tmp;
			result_in_rho[0] += hBlock[j];
			if (j < context.n_Max) tmp *= context.rho;
		}
		for (int64_t i = 1; size_t(i) <= context.lambda; i++)
		{
			result_in_rho[size_t(i)] = 0;
			for (int64_t j = 0; size_t(j) <= context.n_Max; j++)
			{
				hBlock[size_t(j)] *= Delta + Real(j - i + 1);
				result_in_rho[size_t(i)] += hBlock[size_t(j)];
			}
			result_in_rho[size_t(i)] /= mpfr::factorial(uint32_t(i)) * mpfr::pow(context.rho, uint32_t(i));
		}
		return context.rho_to_z_matrix * result_in_rho;
	}

	template <class Real>
	algebra::Vector<Real> _casimir_expander(const algebra::Vector<Real>& realAxisResult, const Real& epsilon,
	                                        const Real& ell, const Real& Delta, const Real& S, const Real& P,
	                                        size_t lambda)
	{
		auto result = algebra::Vector<Real>((lambda + 2) * (lambda + 2) / 4);
		for (size_t i = 0; i <= lambda; ++i) result[i] = realAxisResult[i];

		Real Casimir, val, term;

		/* computing quadratic Casimir times 2 */
		Casimir = (2 * epsilon + ell) * ell + (Delta - 2 * epsilon - 2) * Delta;

		auto get_index = [lambda](int64_t n, int64_t m) { return size_t((int64_t(lambda) + 2 - n) * n + m); };

		for (int64_t i = 1; size_t(i) <= lambda / 2; i++)
		{
			for (int64_t j = 0; size_t(j + 2 * i) <= lambda; ++j)
			{
				val = 0;
				/* The first line in arXiv 1203.6064 (C.1)
				 * note: n(there) = i (here), and m(there) = j (here)
				 * (a / 2) = x
				 * (b / 2) = sqrt(y)
				 * d/da ^m d/db ^n (g) = h_{n, m} = (1 / 2) ^ (m + n) d/dx ^m d/dy ^n g
				 * and h_{m, n} = (1 / 2) ^ (i + j) * i! j! h_{i, j}
				 * */
				term = 0;
				if (j >= 3) term += 8 * result[get_index(i, j - 3)];
				if (j >= 2) term += 4 * result[get_index(i, j - 2)];
				if (j >= 1) term -= 2 * result[get_index(i, j - 1)];
				term *= 4 * epsilon + Real(4 * i - 2);
				val += term;

				/* The second line */
				if (i >= 1) val += result[get_index(i - 1, j + 2)] * Real(-(j + 2) * (j + 1)) / Real(i);

				if (i >= 1)
				{
					term = 2 * epsilon - Real(j + 4 * i) + 6 + 2 * S;
					term *= result[get_index(i - 1, j + 1)];
					term *= Real(2 * (j + 1));
					term /= Real(i);
					val += term;
				}

				/* The third line */
				if (i >= 1)
				{
					term = epsilon * Real(4 * (j + i - 1));
					term += Real(j * (j + 8 * i - 5) + i * (4 * i - 2) - 2);
					term += 2 * P;
					term += S * Real(8 * i + 4 * j - 8);
					term += 2 * Casimir;
					term *= result[get_index(i - 1, j)];
					term *= 4;
					term /= Real(i);
					val += term;
				}

				/* The fourth line */
				if (i >= 1 && j >= 1)
				{
					term = epsilon * Real(2 * (j + 1 - 2 * i));
					term += Real(j * (j + 12 * i - 13) + i * (12 * i - 34) + 22);
					term += 2 * P;
					term += S * Real(8 * i + 2 * j - 10);
					term *= result[get_index(i - 1, j - 1)];
					term *= 8;
					term /= Real(i);
					val += term;
				}

				/* The last line */
				if (i >= 2)
				{
					term = 2 * epsilon;
					term += Real(6 - 3 * j - 4 * i);
					term += -2 * S;
					term *= result[get_index(i - 2, j + 1)];
					term *= Real(8 * (j + 1));
					term -= result[get_index(i - 2, j + 2)] * Real(4 * (j + 1) * (j + 2));
					term /= Real(i);
					val -= term;
				}

				/* finally division by 2 (D + 2 n - 3) */
				result[get_index(i, j)] = val / (4 * epsilon + Real(4 * i - 2));
			}
		}
		return result;
	}

	template <class Real>
	algebra::Vector<Real> gBlock_full(const Real& epsilon, const Real& ell, const Real& Delta, const Real& S,
	                                  const Real& P, const cb_context<Real>& context)
	{
		auto realAxisResult = _real_axis_result(epsilon, ell, Delta, S, P, context);
		return _casimir_expander(realAxisResult, epsilon, ell, Delta, S, P, context.lambda);
	}

	template <class Real>
	algebra::Vector<Real> hBlock_times_rho_n(size_t n, const Real& epsilon, const Real& ell, const Real& Delta,
	                                         const Real& S, const Real& P, const cb_context<Real>& context)
	{
		/* *
		 * gives (4 * rho) ^ {n} * h(\Delta, l,...) = (4 * rho) ^ {n - \Delta} g * (\Delta, l, ...)
		 * , evaluated in x-y coordinate
		 * */
		auto hBlock = _recursion_vector(context.n_Max, epsilon, ell, Delta, S, P);
		algebra::Vector<Real> result_in_rho(context.lambda + 1);
		Real tmp;

		result_in_rho[0] = 0;
		tmp = mpfr::pow(4 * context.rho, int32_t(n));
		for (size_t j = 0; j <= context.n_Max; j++)
		{
			hBlock[j] *= tmp;
			result_in_rho[0] += hBlock[j];
			if (j < context.n_Max) tmp *= context.rho;
		}
		for (int64_t i = 1; size_t(i) <= context.lambda; i++)
		{
			result_in_rho[size_t(i)] = 0;
			for (int64_t j = 0; size_t(j) <= context.n_Max; j++)
			{
				hBlock[size_t(j)] *= Real(int64_t(n) + j - i + 1);
				result_in_rho[size_t(i)] += hBlock[size_t(j)];
			}
			result_in_rho[size_t(i)] /= mpfr::pow(context.rho, int32_t(i)) * mpfr::factorial(uint32_t(i));
		}

		return context.rho_to_z_matrix * result_in_rho;
	}

	template <class Real>
	algebra::Vector<Real> h_asymptotic(const Real& epsilon, const Real& S, const cb_context<Real>& context)
	{
		auto getFactor = [&epsilon, &context, &S](int sign) {
			Real temp1, temp2;
			temp1 = (2 * sign) * S - 1 - epsilon;
			algebra::Vector<Real> factor(context.lambda + 1);
			temp2 = 1 + sign * context.rho;
			factor[0] = mpfr::pow(temp2, temp1);
			temp2 = sign / temp2;
			for (size_t j = 1; j <= context.lambda; j++)
			{
				factor[j] = factor[j - 1] * (temp1 - Real(j - 1));
				factor[j] *= temp2;
				factor[j] /= Real(j);
			}
			return factor;
		};
		auto firstFactor = getFactor(1);
		auto secondFactor = getFactor(-1);
		auto result_in_rho = firstFactor.convolve(secondFactor);
		return context.rho_to_z_matrix * result_in_rho;
	}
}  // namespace qboot

#endif  // HOR_RECURSION_H
