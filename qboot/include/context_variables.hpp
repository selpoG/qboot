#ifndef QBOOT_CONTEXT_VARIABLES_HPP_
#define QBOOT_CONTEXT_VARIABLES_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t, int32_t
#include <sstream>  // for ostringstream
#include <string>   // for string

#include "block.hpp"             // for ConformalBlock
#include "complex_function.hpp"  // for ComplexFunction, FunctionSymmetry
#include "hor_recursion.hpp"     // for gBlock
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for mpfr_prec_t, mpfr_rnd_t, mpfr_t, real, sqrt
#include "real_function.hpp"     // for RealFunction, RealConverter

namespace qboot
{
	// z = 4 r / (1 + r) ^ 2
	// calculate z - 1 / 2 as a function of r' (= r - 3 + 2 sqrt(2)) upto r' ^ {lambda}
	algebra::RealFunction<mpfr::real> z_as_func_rho(uint32_t lambda);

	// controls n_Max, lambda and dim
	class Context
	{
		uint32_t n_Max_, lambda_, dim_;
		mpfr::real epsilon_, rho_;
		// convert a function of rho - (3 - 2 sqrt(2)) to a function of z - 1 / 2
		algebra::RealConverter<mpfr::real> rho_to_z_;

	public:
		// cutoff of the power series expansion of conformal blocks at rho = 0
		[[nodiscard]] uint32_t n_Max() const { return n_Max_; }
		// passed to the constructor of ComplexFunction<Real>
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		// the dimension of the spacetime
		[[nodiscard]] uint32_t dimension() const { return dim_; }
		// (dim - 2) / 2
		[[nodiscard]] const mpfr::real& epsilon() const { return epsilon_; }
		// 3 - 2 sqrt(2)
		[[nodiscard]] const mpfr::real& rho() const { return rho_; }
		// convert a function of rho - (3 - 2 sqrt(2)) to a function of z - 1 / 2
		[[nodiscard]] const algebra::RealConverter<mpfr::real>& rho_to_z() const { return rho_to_z_; }
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "Context(nMax=" << n_Max_ << ", lambda=" << lambda_ << ", dim=" << dim_ << ")";
			return os.str();
		}
		Context(uint32_t n_Max, uint32_t lambda, uint32_t dim);
		Context(Context&&) noexcept = default;
		Context& operator=(Context&&) noexcept = default;
		Context(const Context&) = delete;
		Context& operator=(const Context&) = delete;
		~Context() = default;
		[[nodiscard]] mpfr::real unitarity_bound(uint32_t spin) const { return qboot::unitarity_bound(epsilon_, spin); }
		// calculate v ^ d gBlock and project to sym-symmetric part
		// F-type corresponds to sym = Odd, H-type to Even
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> F_block(const mpfr::real& d,
		                                                           const algebra::ComplexFunction<mpfr::real>& gBlock,
		                                                           algebra::FunctionSymmetry sym) const
		{
			return mul(algebra::v_to_d(d, lambda_), gBlock).proj(sym);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> evaluate(const ConformalBlock<PrimaryOperator>& block) const
		{
			return F_block(block.delta_half(), gBlock(block.get_op(), block.S(), block.P(), *this), block.symmetry());
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> evaluate(const ConformalBlock<GeneralPrimaryOperator>& block,
		                                                            const mpfr::real& delta) const
		{
			return F_block(block.delta_half(), gBlock(block.get_op(delta), block.S(), block.P(), *this),
			               block.symmetry());
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> F_block(const PrimaryOperator& op, const mpfr::real& d1,
		                                                           const mpfr::real& d2, const mpfr::real& d3,
		                                                           const mpfr::real& d4,
		                                                           algebra::FunctionSymmetry sym) const
		{
			return F_block((d2 + d3) / 2, gBlock(op, d1, d2, d3, d4, *this), sym);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> F_block(
		    const mpfr::real& d2, const mpfr::real& d3, const algebra::ComplexFunction<mpfr::real>& gBlock) const
		{
			return F_block((d2 + d3) / 2, gBlock);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> F_block(
		    const mpfr::real& d, const algebra::ComplexFunction<mpfr::real>& gBlock) const
		{
			return F_block(d, gBlock, algebra::FunctionSymmetry::Odd);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> F_block(const PrimaryOperator& op, const mpfr::real& d1,
		                                                           const mpfr::real& d2, const mpfr::real& d3,
		                                                           const mpfr::real& d4) const
		{
			return F_block(op, d1, d2, d3, d4, algebra::FunctionSymmetry::Odd);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> H_block(
		    const mpfr::real& d2, const mpfr::real& d3, const algebra::ComplexFunction<mpfr::real>& gBlock) const
		{
			return H_block((d2 + d3) / 2, gBlock);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> H_block(
		    const mpfr::real& d, const algebra::ComplexFunction<mpfr::real>& gBlock) const
		{
			return F_block(d, gBlock, algebra::FunctionSymmetry::Even);
		}
		[[nodiscard]] algebra::ComplexFunction<mpfr::real> H_block(const PrimaryOperator& op, const mpfr::real& d1,
		                                                           const mpfr::real& d2, const mpfr::real& d3,
		                                                           const mpfr::real& d4) const
		{
			return F_block(op, d1, d2, d3, d4, algebra::FunctionSymmetry::Even);
		}
		// recover off-diagonal derivatives of conformal block from the diagonal derivatives.
		// use recursion relation of eq (2.17) in arXiv:1602.02810 (generalized ver. of eq (C.1) in arXiv:1203.6064).
		// note:
		//   z   = x + sqrt(y) = (a + sqrt(b)) / 2
		//   z^* = x - sqrt(y) = (a - sqrt(b)) / 2
		//   a = 2 x, b = 4 y
		//   h_{m, n} := (der a) ^ m (der b) ^ n g_{Delta, spin}
		//             = (1 / 2) ^ {m + 2 n} (der x) ^ m (der y) ^ n g_{Delta, spin}
		//   and h_{m, n} = m! n! f[m, n] / 2 ^ {m + 2 n}
		algebra::ComplexFunction<mpfr::real> expand_off_diagonal(algebra::RealFunction<mpfr::real>&& realAxisResult,
		                                                         const PrimaryOperator& op, const mpfr::real& S,
		                                                         const mpfr::real& P) const;
	};
}  // namespace qboot

#endif  // QBOOT_CONTEXT_VARIABLES_HPP_
