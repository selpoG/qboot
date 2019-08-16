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
#include "real_function.hpp"     // for RealFunction, RealConverter

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
	// z = 4 r / (1 + r) ^ 2
	// calculate z - 1 / 2 as a function of r' (= r - 3 + 2 sqrt(2)) upto r' ^ {lambda}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	RealFunction<Real> z_as_func_rho(uint32_t lambda)
	{
		// z - 1 / 2 = (-r' ^ 2 / 2 + 2 sqrt(2) r') (4 - 2sqt(2) + r') ^ {-2}
		// f = (4 - 2sqt(2) + r') ^ {-2}
		auto f = qboot::power_function<Real>(4 - mpfr::sqrt(Real(8)), Real(1), Real(-2), lambda);
		// g1 = 2 sqrt(2) r' f
		auto g1 = f.clone();
		g1.shift(1);
		g1 *= mpfr::sqrt(Real(8));
		// g2 = -r' ^ 2 f / 2
		auto g2 = f.clone();
		g2.shift(2);
		g2 /= -2;
		// z - 1 / 2 = g1 + g2
		return g1 + g2;
	}
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
		// convert a function of rho - (3 - 2 sqrt(2)) to a function of z - 1 / 2
		RealConverter<Real> rho_to_z;
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
		      rho_to_z(RealConverter<Real>(z_as_func_rho<Real>(lambda)).inverse()),
		      rho_to_delta(lambda + 1)
		{
			assert(dim >= 3 && dim % 2 == 1);

			// (rho + r') ^ Delta = rho ^ Delta \sum_{i = 0}^{lambda} rho_to_delta[i] r' ^ i + O(r' ^ {lambda + 1})
			rho_to_delta[0] = {Real(1)};
			for (uint32_t i = 1; i <= lambda; ++i)
			{
				Real tmp = 1 / (rho * i);
				rho_to_delta[i] = rho_to_delta[i - 1] * algebra::Polynomial<Real>({(1 - int32_t(i)) * tmp, tmp});
			}
			rho_to_delta = rho_to_delta * rho_to_z.matrix();
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
		RealFunction<Real> h_times_rho_k(uint32_t k, const PrimaryOperator<Real>& op, const Real& S,
		                                 const Real& P) const
		{
			return hBlock_powered(Real(k), op, S, P);
		}
		RealFunction<Real> h_asymptotic_form(const Real& S) const { return h_asymptotic(S, *this); }
		ComplexFunction<algebra::Polynomial<Real>> expand_off_diagonal(
		    RealFunction<algebra::Polynomial<Real>>&& realAxisResult, uint32_t spin, const Real& S, const Real& P) const
		{
			return expand_off_diagonal(std::move(realAxisResult), get_general_primary(spin), S, P);
		}
		// recover off-diagonal derivatives of conformal block from the diagonal derivatives.
		// use recursion relation of eq (2.17) in arXiv:1602.02810 (generalized ver. of eq (C.1) in arXiv:1203.6064).
		// note:
		//   z   = x + sqrt(y) = (a + sqrt(b)) / 2
		//   z^* = x - sqrt(y) = (a - sqrt(b)) / 2
		//   a = 2 x, b = 4 y
		//   h_{m, n} := (der a) ^ m (der b) ^ n g_{\Delta, spin}
		//             = (1 / 2) ^ {m + 2 n} (der x) ^ m (der y) ^ n g_{\Delta, spin}
		//   and h_{m, n} = m! n! f[m, n] / 2 ^ {m + 2 n}
		template <class T>
		ComplexFunction<T> expand_off_diagonal(RealFunction<T>&& realAxisResult, const PrimaryOperator<Real, T>& op,
		                                       const Real& S, const Real& P) const
		{
			assert(realAxisResult.lambda() == lambda);
			ComplexFunction<T> f(lambda);
			for (uint32_t m = 0; m <= lambda; ++m) f.get(m, 0u) = std::move(realAxisResult.get(m));

			T val{}, term{}, quad_casimir = op.quadratic_casimir();

			for (int32_t n = 1; uint32_t(n) <= lambda / 2; ++n)
			{
				// multiply n * (4 * epsilon + 4 * n - 2) and divide
				Real common_factor = (4 * n) * epsilon + (4 * n - 2) * n;
				for (int32_t m = 0; uint32_t(m + 2 * n) <= lambda; ++m)
				{
					// h[m + 2, n - 1]
					val = (-(m + 1) * (m + 2)) * f.get(m + 2, n - 1);

					// h[m + 1, n - 1]
					term = T(2 * (epsilon + S) - (m + 4 * n - 6));
					term *= 2 * (m + 1);
					term *= f.get(m + 1, n - 1);
					val += term;

					// h[m, n - 1]
					term = T(epsilon * (4 * (m + n - 1)));
					term += m * (m + 8 * n - 5) + n * (4 * n - 2) - 2;
					term += 2 * P;
					term += S * (8 * n + 4 * m - 8);
					term += 4 * quad_casimir;
					term *= 4;
					term *= f.get(m, n - 1);
					val += term;

					// h[m - 1, n - 1]
					if (m >= 1)
					{
						term = T(epsilon * (2 * (m + 1 - 2 * n)));
						term += m * (m + 12 * n - 13) + n * (12 * n - 34) + 22;
						term += 2 * P;
						term += S * (8 * n + 2 * m - 10);
						term *= 8;
						term *= f.get(m - 1, n - 1);
						val += term;
					}

					// h[m + 1, n - 2] and h[m + 2, n - 2]
					if (n >= 2)
					{
						term = T(-2 * epsilon);
						term += 3 * m + 4 * n - 6;
						term += 2 * S;
						term *= 8 * (m + 1);
						term *= f.get(m + 1, n - 2);
						term += f.get(m + 2, n - 2) * (4 * (m + 1) * (m + 2));
						val += term;
					}

					val /= common_factor;

					// h[m - 1, n], h[m - 2, n] and h[m - 3, n]
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