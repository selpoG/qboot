#include "context.hpp"

using algebra::RealFunction, algebra::ComplexFunction, algebra::RealConverter;
using mpfr::real;

namespace qboot
{
	RealFunction<real> z_as_func_rho(uint32_t lambda)
	{
		// z - 1 / 2 = (-r' ^ 2 / 2 + 2 sqrt(2) r') (4 - 2sqt(2) + r') ^ {-2}
		// f = (4 - 2sqt(2) + r') ^ {-2}
		auto f = algebra::power_function(4 - mpfr::sqrt(8), real(1), real(-2), lambda);
		// g1 = 2 sqrt(2) r' f
		auto g1 = f.clone();
		g1.shift(1);
		g1 *= mpfr::sqrt(8);
		// g2 = -r' ^ 2 f / 2
		auto g2 = f.clone();
		g2.shift(2);
		g2 /= -2;
		// z - 1 / 2 = g1 + g2
		return g1 + g2;
	}
	Context::Context(uint32_t n_Max, uint32_t lambda, uint32_t dim)
	    : n_Max_(n_Max),
	      lambda_(lambda),
	      dim_(dim),
	      epsilon_(real(dim) / 2 - 1),
	      rho_(3 - mpfr::sqrt(8)),
	      rho_to_z_(RealConverter(z_as_func_rho(lambda)).inverse())
	{
		assert(dim >= 3 && dim % 2 == 1);
	}
	ComplexFunction<real> Context::expand_off_diagonal(RealFunction<real>&& realAxisResult, const PrimaryOperator& op,
	                                                   const real& S, const real& P) const
	{
		assert(realAxisResult.lambda() == lambda_);
		ComplexFunction<real> f(lambda_);
		for (uint32_t m = 0; m <= lambda_; ++m) f.at(m, 0u).swap(realAxisResult.at(m));

		real val{}, term{}, quad_casimir = op.quadratic_casimir();

		for (int32_t n = 1; uint32_t(n) <= lambda_ / 2; ++n)
		{
			// multiply n (4 epsilon + 4 n - 2) and divide
			real common_factor = (4 * n) * epsilon_ + (4 * n - 2) * n;
			for (int32_t m = 0; uint32_t(m + 2 * n) <= lambda_; ++m)
			{
				// h[m + 2, n - 1]
				val = (-(m + 1) * (m + 2)) * f.at(uint32_t(m + 2), uint32_t(n - 1));

				// h[m + 1, n - 1]
				term = 2 * (epsilon_ + S) - (m + 4 * n - 6);
				term *= 2 * (m + 1);
				term *= f.at(uint32_t(m + 1), uint32_t(n - 1));
				val += term;

				// h[m, n - 1]
				term = epsilon_ * (4 * (m + n - 1));
				term += m * (m + 8 * n - 5) + n * (4 * n - 2) - 2;
				term += 2 * P;
				term += S * (8 * n + 4 * m - 8);
				term += 4 * quad_casimir;
				term *= 4;
				term *= f.at(uint32_t(m), uint32_t(n - 1));
				val += term;

				// h[m - 1, n - 1]
				if (m >= 1)
				{
					term = epsilon_ * (2 * (m + 1 - 2 * n));
					term += m * (m + 12 * n - 13) + n * (12 * n - 34) + 22;
					term += 2 * P;
					term += S * (8 * n + 2 * m - 10);
					term *= 8;
					term *= f.at(uint32_t(m - 1), uint32_t(n - 1));
					val += term;
				}

				// h[m + 1, n - 2] and h[m + 2, n - 2]
				if (n >= 2)
				{
					term = -2 * epsilon_;
					term += 3 * m + 4 * n - 6;
					term += 2 * S;
					term *= 8 * (m + 1);
					term *= f.at(uint32_t(m + 1), uint32_t(n - 2));
					term += f.at(uint32_t(m + 2), uint32_t(n - 2)) * (4 * (m + 1) * (m + 2));
					val += term;
				}

				val /= common_factor;

				// h[m - 1, n], h[m - 2, n] and h[m - 3, n]
				term = {};
				if (m >= 3) term += 8 * f.at(uint32_t(m - 3), uint32_t(n));
				if (m >= 2) term += 4 * f.at(uint32_t(m - 2), uint32_t(n));
				if (m >= 1) term -= 2 * f.at(uint32_t(m - 1), uint32_t(n));

				f.at(uint32_t(m), uint32_t(n)) = val + term;
			}
		}
		return f;
	}
}  // namespace qboot
