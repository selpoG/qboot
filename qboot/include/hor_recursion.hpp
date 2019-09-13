#ifndef QBOOT_HOR_RECURSION_HPP_
#define QBOOT_HOR_RECURSION_HPP_

#include <cstdint>  // for uint32_t

#include "complex_function.hpp"  // for ComplexFunction
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for real, pow
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

	// a function of z - 1 / 2 expanded at z = 1 / 2, g_{Delta, spin}^{d12, d34}(z, z)
	algebra::RealFunction<mpfr::real> gBlock_real(const PrimaryOperator& op, const mpfr::real& S, const mpfr::real& P,
	                                              const Context& context);

	algebra::ComplexFunction<mpfr::real> gBlock(const PrimaryOperator& op, const mpfr::real& S, const mpfr::real& P,
	                                            const Context& context);

	algebra::ComplexFunction<mpfr::real> gBlock(const PrimaryOperator& op, const mpfr::real& d1, const mpfr::real& d2,
	                                            const mpfr::real& d3, const mpfr::real& d4, const Context& context);
}  // namespace qboot

#endif  // QBOOT_HOR_RECURSION_HPP_
