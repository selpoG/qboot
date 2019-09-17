#ifndef QBOOT_BLOCK_HPP_
#define QBOOT_BLOCK_HPP_

#include <cstdint>      // for uint32_t
#include <string>       // for string
#include <type_traits>  // enable_if, is_same_v
#include <variant>      // for varinat

#include "complex_function.hpp"  // for FunctionSymmetry
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for real

namespace qboot
{
	// F_{\mp, op}^{d1 d2, d3 d4}
	template <class Operator>
	class ConformalBlock
	{
		Operator op_;
		mpfr::real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		// sym must be Even or Odd (Mixed is not allowed)
		ConformalBlock(const Operator& op, const mpfr::real& d12, const mpfr::real& d34, const mpfr::real& d23h,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : op_(op), d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = (d34_ - d12_) / 2;
			P_ = d12_ * d34_ / (-2);
		}
		ConformalBlock(const Operator& op, const mpfr::real& d1, const mpfr::real& d2, const mpfr::real& d3,
		               const mpfr::real& d4, algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, d1 - d2, d3 - d4, (d2 + d3) / 2, sym)
		{
		}
		ConformalBlock(const Operator& op, const PrimaryOperator& op1, const PrimaryOperator& op2,
		               const PrimaryOperator& op3, const PrimaryOperator& op4,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, op1.delta(), op2.delta(), op3.delta(), op4.delta(), sym)
		{
		}
		template <class = std::enable_if<std::is_same_v<Operator, PrimaryOperator>>>
		[[nodiscard]] const mpfr::real& delta() const
		{
			return op_.delta();
		}
		[[nodiscard]] uint32_t spin() const { return op_.spin(); }
		[[nodiscard]] const mpfr::rational& epsilon() const { return op_.epsilon(); }
		[[nodiscard]] const mpfr::real& delta_half() const { return d23h_; }
		[[nodiscard]] const mpfr::real& S() const { return S_; }
		[[nodiscard]] const mpfr::real& P() const { return P_; }
		[[nodiscard]] bool include_odd() const { return d12_ != 0 && d34_ != 0; }
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] const Operator& get_op() const { return op_; }
		template <class = std::enable_if<!std::is_same_v<Operator, PrimaryOperator>>>
		[[nodiscard]] PrimaryOperator get_op(const mpfr::real& delta) const
		{
			return op_.fix_delta(delta);
		}
		template <class = std::enable_if<!std::is_same_v<Operator, PrimaryOperator>>>
		[[nodiscard]] ConformalBlock<PrimaryOperator> fix_delta(const mpfr::real& delta) const
		{
			return ConformalBlock<PrimaryOperator>(op_.fix_delta(delta), d12_, d34_, d23h_, sym_);
		}
		[[nodiscard]] std::string str() const;
	};
	class GeneralConformalBlock
	{
		mpfr::real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		// sym must be Even or Odd (Mixed is not allowed)
		GeneralConformalBlock(const mpfr::real& d12, const mpfr::real& d34, const mpfr::real& d23h,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = (d34_ - d12_) / 2;
			P_ = d12_ * d34_ / (-2);
		}
		GeneralConformalBlock(const mpfr::real& d1, const mpfr::real& d2, const mpfr::real& d3, const mpfr::real& d4,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : GeneralConformalBlock(d1 - d2, d3 - d4, (d2 + d3) / 2, sym)
		{
		}
		GeneralConformalBlock(const PrimaryOperator& op1, const PrimaryOperator& op2, const PrimaryOperator& op3,
		                      const PrimaryOperator& op4,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : GeneralConformalBlock(op1.delta(), op2.delta(), op3.delta(), op4.delta(), sym)
		{
		}
		[[nodiscard]] const mpfr::real& delta_half() const { return d23h_; }
		[[nodiscard]] const mpfr::real& S() const { return S_; }
		[[nodiscard]] const mpfr::real& P() const { return P_; }
		[[nodiscard]] bool include_odd() const { return d12_ != 0 && d34_ != 0; }
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		template <class Operator>
		[[nodiscard]] ConformalBlock<Operator> fix_op(const Operator& op) const
		{
			return ConformalBlock<Operator>(op, d12_, d34_, d23h_, sym_);
		}
		[[nodiscard]] std::string str() const;
	};
	using Block = std::variant<ConformalBlock<PrimaryOperator>, GeneralConformalBlock>;
}  // namespace qboot

#endif  // QBOOT_BLOCK_HPP_
