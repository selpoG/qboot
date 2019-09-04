#ifndef QBOOT_BLOCK_HPP_
#define QBOOT_BLOCK_HPP_

#include <cstdint>      // for uint32_t
#include <sstream>      // for ostringstream
#include <string>       // for string
#include <type_traits>  // enable_if

#include "complex_function.hpp"  // for FunctionSymmetry
#include "primary_op.hpp"        // for PrimaryOperator
#include "real.hpp"              // for real

namespace qboot
{
	// F_{\mp, op}^{d1 d2, d3 d4}
	template <class Real, class Operator>
	class ConformalBlock
	{
		Operator op_;
		Real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		// sym must be Even or Odd (Mixed is not allowed)
		ConformalBlock(const Operator& op, const Real& d12, const Real& d34, const Real& d23h,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : op_(op), d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = (d34_ - d12_) / 2;
			P_ = d12_ * d34_ / (-2);
		}
		ConformalBlock(const Operator& op, const Real& d1, const Real& d2, const Real& d3, const Real& d4,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, d1 - d2, d3 - d4, (d2 + d3) / 2, sym)
		{
		}
		ConformalBlock(const Operator& op, const PrimaryOperator<Real>& op1, const PrimaryOperator<Real>& op2,
		               const PrimaryOperator<Real>& op3, const PrimaryOperator<Real>& op4,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, op1.delta(), op2.delta(), op3.delta(), op4.delta(), sym)
		{
		}
		template <class = std::enable_if<is_primary_operator_v<Operator>>>
		[[nodiscard]] const Real& delta() const
		{
			return op_.delta();
		}
		[[nodiscard]] uint32_t spin() const { return op_.spin(); }
		[[nodiscard]] const Real& epsilon() const { return op_.epsilon(); }
		[[nodiscard]] const Real& delta_half() const { return d23h_; }
		[[nodiscard]] const Real& S() const { return S_; }
		[[nodiscard]] const Real& P() const { return P_; }
		[[nodiscard]] bool include_odd() const { return d12_ != 0 && d34_ != 0; }
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] const Operator& get_op() const { return op_; }
		template <class = std::enable_if<!is_primary_operator_v<Operator>>>
		[[nodiscard]] PrimaryOperator<Real> get_op(const Real& delta) const
		{
			return op_.fix_delta(delta);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "ConformalBlock(op=" << op_.str() << ", d12=" << d12_ << ", d34=" << d34_ << ", d23h=" << d23h_
			   << ", type=" << (sym_ == algebra::FunctionSymmetry::Odd ? 'F' : 'H') << ")";
			return os.str();
		}
	};
	template <class Real>
	class GeneralConformalBlock
	{
		Real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		// sym must be Even or Odd (Mixed is not allowed)
		GeneralConformalBlock(const Real& d12, const Real& d34, const Real& d23h,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = (d34_ - d12_) / 2;
			P_ = d12_ * d34_ / (-2);
		}
		GeneralConformalBlock(const Real& d1, const Real& d2, const Real& d3, const Real& d4,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : GeneralConformalBlock(d1 - d2, d3 - d4, (d2 + d3) / 2, sym)
		{
		}
		GeneralConformalBlock(const PrimaryOperator<Real>& op1, const PrimaryOperator<Real>& op2,
		                      const PrimaryOperator<Real>& op3, const PrimaryOperator<Real>& op4,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : GeneralConformalBlock(op1.delta(), op2.delta(), op3.delta(), op4.delta(), sym)
		{
		}
		[[nodiscard]] const Real& delta_half() const { return d23h_; }
		[[nodiscard]] const Real& S() const { return S_; }
		[[nodiscard]] const Real& P() const { return P_; }
		[[nodiscard]] bool include_odd() const { return d12_ != 0 && d34_ != 0; }
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		template <class Operator>
		[[nodiscard]] ConformalBlock<Real, Operator> fix_op(const Operator& op) const
		{
			return ConformalBlock<Real, Operator>(op, d12_, d34_, d23h_, sym_);
		}
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "ConformalBlock(d12=" << d12_ << ", d34=" << d34_ << ", d23h=" << d23h_
			   << ", type=" << (sym_ == algebra::FunctionSymmetry::Odd ? 'F' : 'H') << ")";
			return os.str();
		}
	};
	template <class Real>
	using Block = std::variant<ConformalBlock<Real, PrimaryOperator<Real>>, GeneralConformalBlock<Real>>;
}  // namespace qboot

#endif  // QBOOT_BLOCK_HPP_
