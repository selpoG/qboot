#ifndef CONTEXT_VARIABLES_HPP_
#define CONTEXT_VARIABLES_HPP_

#include <array>    // for array
#include <cstdint>  // for uint32_t, int32_t
#include <memory>   // for unique_ptr
#include <sstream>  // for ostringstream
#include <string>   // for string
#include <vector>   // for vector

#include "complex_function.hpp"  // for ComplexFunction
#include "matrix.hpp"            // for Vector, Matrix
#include "polynomial.hpp"        // for Polynomial
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for mpfr_prec_t, mpfr_rnd_t, mpfr_t, real, sqrt, pow

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

	// basic constructor for cb_context
	cb_context context_construct(uint32_t n_Max, mpfr_prec_t prec, uint32_t lambda);
	void clear_cb_context(cb_context& context);
}  // namespace qboot2

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class RationalApproxData;
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Context
	{
	public:
		using type = Real;
		static constexpr mpfr_prec_t prec = Real::prec;
		static constexpr mpfr_rnd_t rnd = Real::rnd;
		const uint32_t n_Max, lambda, dimension;
		const Real epsilon, rho;
		algebra::Matrix<Real> rho_to_z_matrix;
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
		      rho_to_delta(lambda + 1)
		{
			assert(dim >= 3 && dim % 2 == 1);
			auto Lambda = lambda + 1;

			// z = 4 * rho / pow(1 + rho, 2)
			// (d/dz)^n / n! = sum(rho_to_z_matrix[n, i] * (d/drho)^i / i!, 0 <= i <= n)
			algebra::Vector<Real> tmps(Lambda);
			tmps[0] = -mpfr::sqrt(Real(8));
			for (uint32_t j = 1; j < Lambda; ++j) tmps[j] = (tmps[j - 1] * (2 * int32_t(j) - 3)) / j;
			tmps[1] -= 2;
			tmps[0] += 3;

			Real tmp;
			rho_to_z_matrix.get(0, 0) = 1;
			for (uint32_t j = 1; j < Lambda; ++j)
			{
				tmp = 1;
				for (uint32_t k = 0; k <= j; ++k)
				{
					rho_to_z_matrix.get(1, j) += tmps[j - k] * tmp;
					tmp *= -2;
				}
			}
			for (uint32_t i = 2; i < Lambda; ++i)
				for (uint32_t j = 1; j < Lambda; ++j)
					for (uint32_t k = i - 1; k < Lambda - j; ++k)
						rho_to_z_matrix.get(i, k + j) += rho_to_z_matrix.get(i - 1, k) * rho_to_z_matrix.get(1, j);

			rho_to_z_matrix.transpose();

			rho_to_delta[0] = {Real(1)};
			for (uint32_t i = 1; i <= lambda; ++i)
			{
				tmp = 1 / (rho * i);
				rho_to_delta[i] = rho_to_delta[i - 1] * algebra::Polynomial<Real>({(1 - int32_t(i)) * tmp, tmp});
			}
			rho_to_delta = rho_to_z_matrix * rho_to_delta;
		}
		Context(Context&&) = default;
		Context& operator=(Context&&) = default;
		Context(const Context&) = delete;
		Context& operator=(const Context&) = delete;
		~Context() = default;
		template <class T>
		Real operator()(const T& val) const
		{
			return Real(val);
		}
		Real unitary_bound(uint32_t spin) const { return spin == 0 ? epsilon : spin + 2 * epsilon; }
		auto get_primary(const Real& delta, uint32_t spin) const { return PrimaryOperator(delta, spin, *this); }
		auto get_general_primary(uint32_t spin) const { return general_primary_operator(spin, *this); }
		// calculate v ^ d gBlock and project to sym-symmetric part
		// F-type corresponds to sym = Odd, H-type to Even
		ComplexFunction<Real> F_block(const Real& d, const ComplexFunction<Real>& gBlock, FunctionSymmetry sym) const
		{
			return (v_to_d(d, lambda) * gBlock).proj(sym);
		}
		algebra::Vector<Real> h_times_rho_k(uint32_t k, const PrimaryOperator<Real>& op, const Real& S,
		                                    const Real& P) const
		{
			return hBlock_powered(Real(k), op, S, P);
		}
		algebra::Vector<Real> h_asymptotic_form(const Real& S) const { return h_asymptotic(S, *this); }
		ComplexFunction<algebra::Polynomial<Real>> expand_off_diagonal(
		    algebra::Vector<algebra::Polynomial<Real>>&& realAxisResult, uint32_t spin, const Real& S,
		    const Real& P) const
		{
			return expand_off_diagonal(std::move(realAxisResult), get_general_primary(spin), S, P);
		}
		// recover off-diagonal derivatives of conformal block from the diagonal derivatives.
		// use recursion relation of eq (C.1) in arXiv:1203.6064.
		// note:
		//   z   = x + sqrt(y) = (a + sqrt(b)) / 2
		//   z^* = x - sqrt(y) = (a - sqrt(b)) / 2
		//   h_{m, n}(there) := (der a) ^ m (der b) ^ n g_{\Delta, spin}
		//                    = (1 / 2) ^ {m + n} (der x) ^ m (der y) ^ n g_{\Delta, spin}
		//   and h_{m, n}(there) = m! n! h_{m, n}(here) / 2 ^ {m + n}
		//   h_{m, n}(here) = (der x) ^ m (der y) ^ n g_{\Delta, spin} / (m! n!)
		//   h_{m, n}(here) satisfies g_{\Delta, spin} = \sum_{m, n >= 0} h_{m, n}(here) (x - 1 / 2) ^ m y ^ n
		//   ((x, y) = (1 / 2, 0) is the crossing symmetric point)
		//   h_{m, n}(here) is returned as a vector
		//   (h_{0, 0}, ..., h_{lambda, 0}, h_{0, 1}, ..., h_{lambda - 2, 1}, h_{0, 2}, ...)
		template <class T>
		ComplexFunction<T> expand_off_diagonal(algebra::Vector<T>&& realAxisResult, const PrimaryOperator<Real, T>& op,
		                                       const Real& S, const Real& P) const
		{
			ComplexFunction<T> f(lambda);
			for (uint32_t m = 0; m <= lambda; ++m) f.get(m, 0u) = std::move(realAxisResult[m]);

			T val{}, term{}, quad_casimir = op.quadratic_casimir();

			for (int32_t n = 1; uint32_t(n) <= lambda / 2; ++n)
			{
				Real common_factor = (4 * n) * epsilon + (4 * n - 2) * n;
				for (int32_t m = 0; uint32_t(m + 2 * n) <= lambda; ++m)
				{
					// The second line
					val = (-(m + 1) * (m + 2)) * f.get(m + 2, n - 1);

					term = T(2 * (epsilon + S) - (m + 4 * n - 6));
					term *= f.get(m + 1, n - 1);
					term *= 2 * (m + 1);
					val += term;

					// The third line
					term = T(epsilon * (4 * (m + n - 1)));
					term += m * (m + 8 * n - 5) + n * (4 * n - 2) - 2;
					term += 2 * P;
					term += S * (8 * n + 4 * m - 8);
					term += 4 * quad_casimir;
					term *= f.get(m, n - 1);
					term *= 4;
					val += term;

					// The fourth line
					if (m >= 1)
					{
						term = T(epsilon * (2 * (m + 1 - 2 * n)));
						term += m * (m + 12 * n - 13) + n * (12 * n - 34) + 22;
						term += 2 * P;
						term += S * (8 * n + 2 * m - 10);
						term *= f.get(m - 1, n - 1);
						term *= 8;
						val += term;
					}

					// The last line
					if (n >= 2)
					{
						term = T(2 * epsilon);
						term += 6 - 3 * m - 4 * n;
						term += -2 * S;
						term *= f.get(m + 1, n - 2);
						term *= 8 * (m + 1);
						term -= f.get(m + 2, n - 2) * (4 * (m + 1) * (m + 2));
						val -= term;
					}

					// finally divide by 2 * (D + 2 n - 3)
					val /= common_factor;

					// The first line
					term = {};
					if (m >= 3) term += 8 * f.get(m - 3, n);
					if (m >= 2) term += 4 * f.get(m - 2, n);
					if (m >= 1) term -= 2 * f.get(m - 1, n);

					f.get(m, n) = val + term;
				}
			}
			return f;
		}
	};
}  // namespace qboot

#endif  // CONTEXT_VARIABLES_HPP_
