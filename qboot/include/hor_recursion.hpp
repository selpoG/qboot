#ifndef QBOOT_HOR_RECURSION_HPP_
#define QBOOT_HOR_RECURSION_HPP_

#include <cstdint>  // for uint32_t

#include "complex_function.hpp"  // for ComplexFunction
#include "hor_formula.hpp"       // for _get_rec_coeffs
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for real, mpfr_prec_t, mpfr_rnd_t, mpfr_t, pow
#include "real_function.hpp"     // for RealFunction, RealFunctionWithPower, power_function

namespace qboot
{
	class Context;
	// a function f(rho) of rho at rho = 0 (not crossing symmetric point)
	// f(rho) = h_{Delta, spin}^{d12, d34}(z, z),
	// g_{Delta, spin}^{d12, d34}(z, z) = (4 rho) ^ {Delta} f(rho)
	// if p[0] may be 0, we use continuity of conformal block.
	algebra::RealFunction<mpfr::real> hBlock_shifted(const PrimaryOperator& op, const mpfr::real& S,
	                                                 const mpfr::real& P, uint32_t n_Max);

	// a function of z - 1 / 2 expanded at z = 1 / 2,
	// (4 * rho) ^ {exp} * h_{Delta, spin}^{d12, d34}(z, z)
	// = (4 * rho) ^ {exp - Delta} * g_{Delta, spin}^{d12, d34}(z, z)
	algebra::RealFunction<mpfr::real> hBlock_powered(const mpfr::real& exp, const PrimaryOperator& op,
	                                                 const mpfr::real& S, const mpfr::real& P, const Context& context);

	// a function of z - 1 / 2 expanded at z = 1 / 2, g_{Delta, spin}^{d12, d34}(z, z)
	algebra::RealFunction<mpfr::real> gBlock_real(const PrimaryOperator& op, const mpfr::real& S, const mpfr::real& P,
	                                              const Context& context);

	algebra::ComplexFunction<mpfr::real> gBlock(const PrimaryOperator& op, const mpfr::real& S, const mpfr::real& P,
	                                            const Context& context);

	algebra::ComplexFunction<mpfr::real> gBlock(const PrimaryOperator& op, const mpfr::real& d1, const mpfr::real& d2,
	                                            const mpfr::real& d3, const mpfr::real& d4, const Context& context);

	// calculate \tilde{h}(r, 1) as a function of z - 1 / 2 eq (4.6) in arXiv:1406:4858
	// \tilde{h}(r, 1)
	//   = (1 - r ^ 2) ^ {-epsilon} (1 + r) ^ {-1 - d12 + d34} (1 - r) ^ {-1 + d12 - d34}
	//   = (1 + r) ^ {-epsilon - 1 + 2 S} (1 - r) ^ {-epsilon - 1 - 2 S}
	algebra::RealFunction<mpfr::real> h_asymptotic(const mpfr::real& S, const Context& context);
}  // namespace qboot

#endif  // QBOOT_HOR_RECURSION_HPP_
