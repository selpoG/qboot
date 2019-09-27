#ifndef QBOOT_HOR_FORMULA_HPP_
#define QBOOT_HOR_FORMULA_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t

#include "qboot/algebra/matrix.hpp"      // for Vector
#include "qboot/algebra/polynomial.hpp"  // for Polynomial
#include "qboot/mp/real.hpp"             // for real
#include "qboot/primary_op.hpp"          // for PrimaryOperator

// Hogervorst-Osborn-Rycykov recursion relation (derived from eq (4.9) in arXiv:1305.1321) takes the form
//
//   \sum_{i = 0}^{7} b[n - i] p[i] (n, parameters (epsilon, spin, Delta, S, P)) = 0,
//   \sum_{i = 0}^{5} b[n - i] p[i] (n, parameters (epsilon, spin = 0, Delta, S, P)) = 0
//
// and each coefficients p[i] is a polynomial of n whose coefficients are functions of other params.
//
// Let f(z) be the diagonal conformal block g_{Delta, spin}^{d12, d34}(z, z).
// f(z) can be expanded at z = 0 as
//   f(z) = (4 rho) ^ {Delta} \sum_{n = 0}^{\infty} b[n] rho ^ n, b[0] = 1
//   z = 4 rho / (1 + rho) ^ 2

namespace qboot
{
	algebra::Vector<algebra::Polynomial> _get_nonzero_spin_rec_coeffs(const PrimaryOperator& op, const mp::real& S,
	                                                                  const mp::real& P);

	algebra::Vector<algebra::Polynomial> _get_zero_spin_rec_coeffs(const PrimaryOperator& op, const mp::real& S,
	                                                               const mp::real& P);

	// gives coefficient polynomials p[i] of recursion equation of b[i].
	// \sum_{i = 0}^{p.size() - 1} b[n - i] p[i] = 0
	inline algebra::Vector<algebra::Polynomial> _get_rec_coeffs(const PrimaryOperator& op, const mp::real& S,
	                                                            const mp::real& P)
	{
		if (op.spin() == 0) return _get_zero_spin_rec_coeffs(op, S, P);
		return _get_nonzero_spin_rec_coeffs(op, S, P);
	}
}  // namespace qboot

#endif  // QBOOT_HOR_FORMULA_HPP_
