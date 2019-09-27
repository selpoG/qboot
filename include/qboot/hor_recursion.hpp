#ifndef QBOOT_HOR_RECURSION_HPP_
#define QBOOT_HOR_RECURSION_HPP_

#include <cstdint>  // for uint32_t

#include "qboot/algebra/complex_function.hpp"  // for ComplexFunction
#include "qboot/algebra/real_function.hpp"     // for RealFunction, RealFunctionWithPower, power_function
#include "qboot/mp/real.hpp"                   // for real, pow
#include "qboot/primary_op.hpp"                // for PrimaryOperator

namespace qboot
{
	class Context;

	// a function f(rho) of rho at rho = 0 (not crossing symmetric point)
	// f(rho) = h_{Delta, spin}^{d12, d34}(z, z),
	// g_{Delta, spin}^{d12, d34}(z, z) = (4 rho) ^ {Delta} f(rho)
	// if p[0] may be 0, we use continuity of conformal block.
	algebra::RealFunction<mp::real> hBlock_shifted(const PrimaryOperator& op, const mp::real& S, const mp::real& P,
	                                               uint32_t n_Max);

	// a function of z - 1 / 2 expanded at z = 1 / 2, g_{Delta, spin}^{d12, d34}(z, z)
	algebra::RealFunction<mp::real> gBlock_real(const PrimaryOperator& op, const mp::real& S, const mp::real& P,
	                                            const Context& context);

	algebra::ComplexFunction<mp::real> gBlock(const PrimaryOperator& op, const mp::real& S, const mp::real& P,
	                                          const Context& context);

	algebra::ComplexFunction<mp::real> gBlock(const PrimaryOperator& op, const mp::real& d1, const mp::real& d2,
	                                          const mp::real& d3, const mp::real& d4, const Context& context);
}  // namespace qboot

#endif  // QBOOT_HOR_RECURSION_HPP_
