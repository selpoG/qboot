#ifndef QBOOT_BLOCK_HPP_
#define QBOOT_BLOCK_HPP_

#include <concepts>     // for same_as
#include <cstdint>      // for uint32_t
#include <string>       // for string
#include <type_traits>  // enable_if, is_same_v
#include <variant>      // for varinat

#include "qboot/algebra/complex_function.hpp"  // for FunctionSymmetry
#include "qboot/mp/real.hpp"                   // for real
#include "qboot/primary_op.hpp"                // for PrimaryOperator

namespace qboot
{
	inline mp::real _delta_S(const mp::real& d12, const mp::real& d34) { return (d34 - d12) / 2; }
	inline mp::real _delta_S(const mp::real& d1, const mp::real& d2, const mp::real& d3, const mp::real& d4)
	{
		return _delta_S(d1 - d2, d3 - d4);
	}
	inline mp::real _delta_P(const mp::real& d12, const mp::real& d34) { return d12 * d34 / (-2); }
	inline mp::real _delta_P(const mp::real& d1, const mp::real& d2, const mp::real& d3, const mp::real& d4)
	{
		return _delta_P(d1 - d2, d3 - d4);
	}
	class PrimaryOperator;
	class GeneralPrimaryOperator;
	template <class T>
	concept _operator = std::same_as<T, PrimaryOperator> || std::same_as<T, GeneralPrimaryOperator>;
	// F_{\mp, op}^{d1 d2, d3 d4}
	template <_operator Operator>
	class ConformalBlock
	{
		Operator op_;
		mp::real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		void _reset() &&
		{
			std::move(op_)._reset();
			std::move(d12_)._reset();
			std::move(d34_)._reset();
			std::move(d23h_)._reset();
			std::move(S_)._reset();
			std::move(P_)._reset();
			sym_ = algebra::FunctionSymmetry::Mixed;
		}
		// sym must be Even or Odd (Mixed is not allowed)
		ConformalBlock(const Operator& op, const mp::real& d12, const mp::real& d34, const mp::real& d23h,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : op_(op), d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = _delta_S(d12_, d34_);
			P_ = _delta_P(d12_, d34_);
		}
		ConformalBlock(const Operator& op, const mp::real& d1, const mp::real& d2, const mp::real& d3,
		               const mp::real& d4, algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, d1 - d2, d3 - d4, (d2 + d3) / 2, sym)
		{
		}
		ConformalBlock(const Operator& op, const PrimaryOperator& op1, const PrimaryOperator& op2,
		               const PrimaryOperator& op3, const PrimaryOperator& op4,
		               algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : ConformalBlock(op, op1.delta(), op2.delta(), op3.delta(), op4.delta(), sym)
		{
		}
		[[nodiscard]] const mp::real& delta() const requires std::same_as<Operator, PrimaryOperator>
		{
			return op_.delta();
		}
		[[nodiscard]] uint32_t spin() const { return op_.spin(); }
		[[nodiscard]] const mp::rational& epsilon() const { return op_.epsilon(); }
		[[nodiscard]] const mp::real& delta_half() const { return d23h_; }
		[[nodiscard]] const mp::real& S() const { return S_; }
		[[nodiscard]] const mp::real& P() const { return P_; }
		[[nodiscard]] bool include_odd() const { return d12_ != 0 && d34_ != 0; }
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] const Operator& get_op() const { return op_; }
		[[nodiscard]] PrimaryOperator get_op(const mp::real& delta) const
		    requires std::same_as<Operator, GeneralPrimaryOperator>
		{
			return op_.fix_delta(delta);
		}
		[[nodiscard]] ConformalBlock<PrimaryOperator> fix_delta(const mp::real& delta) const
		    requires std::same_as<Operator, GeneralPrimaryOperator>
		{
			return ConformalBlock<PrimaryOperator>(op_.fix_delta(delta), d12_, d34_, d23h_, sym_);
		}
		[[nodiscard]] std::string str() const;
	};
	class GeneralConformalBlock
	{
		mp::real d12_, d34_, d23h_, S_{}, P_{};
		// Odd  => type F (F_{-})
		// Even => type H (F_{+})
		algebra::FunctionSymmetry sym_;

	public:
		void _reset() &&
		{
			std::move(d12_)._reset();
			std::move(d34_)._reset();
			std::move(d23h_)._reset();
			std::move(S_)._reset();
			std::move(P_)._reset();
			sym_ = algebra::FunctionSymmetry::Mixed;
		}
		// sym must be Even or Odd (Mixed is not allowed)
		GeneralConformalBlock(const mp::real& d12, const mp::real& d34, const mp::real& d23h,
		                      algebra::FunctionSymmetry sym = algebra::FunctionSymmetry::Odd)
		    : d12_(d12), d34_(d34), d23h_(d23h), sym_(sym)
		{
			assert(sym == algebra::FunctionSymmetry::Even || sym == algebra::FunctionSymmetry::Odd);
			S_ = _delta_S(d12_, d34_);
			P_ = _delta_P(d12_, d34_);
		}
		GeneralConformalBlock(const mp::real& d1, const mp::real& d2, const mp::real& d3, const mp::real& d4,
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
		[[nodiscard]] const mp::real& delta_half() const { return d23h_; }
		[[nodiscard]] const mp::real& S() const { return S_; }
		[[nodiscard]] const mp::real& P() const { return P_; }
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
