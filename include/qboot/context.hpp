#ifndef QBOOT_CONTEXT_HPP_
#define QBOOT_CONTEXT_HPP_

#include <cassert>     // for assert
#include <cstdint>     // for uint32_t, int32_t
#include <functional>  // for function
#include <future>      // for shared_future
#include <map>         // for map
#include <mutex>       // for mutex, lock_guard
#include <sstream>     // for ostringstream
#include <string>      // for string
#include <tuple>       // for tuple, apply

#include "qboot/algebra/complex_function.hpp"  // for ComplexFunction, FunctionSymmetry
#include "qboot/algebra/real_function.hpp"     // for RealFunction, RealConverter
#include "qboot/block.hpp"                     // for ConformalBlock
#include "qboot/hor_recursion.hpp"             // for gBlock
#include "qboot/mp/rational.hpp"               // for rational
#include "qboot/mp/real.hpp"                   // for real, sqrt
#include "qboot/primary_op.hpp"                // for PrimaryOperator
#include "qboot/task_queue.hpp"                // for _task_queue

namespace qboot
{
	template <class T>
	class _memoized;
	// create a thread-safe momoized function
	// all types of arguments must be comparable
	template <class T, class... TArgs>
	class _memoized<T(TArgs...)>
	{
		mutable std::mutex mutex_{};
		std::function<T(TArgs...)> f_;
		mutable std::map<std::tuple<TArgs...>, std::shared_future<T>> memo_{};
		mutable _task_queue tq_;

	public:
		explicit _memoized(std::function<T(TArgs...)> f, uint32_t p = 1) : f_(f), tq_(p) {}
		const T& operator()(const TArgs&... args) const
		{
			auto key = std::tuple<TArgs...>(args...);
			{
				std::lock_guard<std::mutex> guard(mutex_);
				if (memo_.find(key) == memo_.end())
					memo_.emplace(key, tq_.push([this, key] { return std::apply(f_, key); }).share());
			}
			return memo_.at(key).get();
		}
		[[nodiscard]] uint32_t _total_memory() const
		{
			std::lock_guard<std::mutex> guard(mutex_);
			uint32_t sum = 0;
			for (const auto& t : memo_) sum += t.second.get().size();
			return sum;
		}
	};
	template <class T, class... TArgs>
	_memoized(T (*)(TArgs...))->_memoized<T(TArgs...)>;

	// controls n_Max, lambda and dim
	// all gBlocks are evaluated in parallel and memoized
	class Context
	{
		uint32_t n_Max_, lambda_, parallel_;
		mp::rational dim_, epsilon_;
		mp::real rho_;
		// convert a function of rho - (3 - 2 sqrt(2)) to a function of z - 1 / 2
		algebra::RealConverter rho_to_z_;
		_memoized<algebra::ComplexFunction<mp::real>(mp::real)> v_to_d_;
		_memoized<algebra::ComplexFunction<mp::real>(PrimaryOperator, mp::real, mp::real)> gBlock_;
		[[nodiscard]] const algebra::ComplexFunction<mp::real>& v_to_d(const mp::real& d) const { return v_to_d_(d); }
		[[nodiscard]] const algebra::ComplexFunction<mp::real>& gBlock(const PrimaryOperator& op, const mp::real& S,
		                                                               const mp::real& P) const
		{
			return gBlock_(op, S, P);
		}

	public:
		void _reset() &&
		{
			n_Max_ = lambda_ = parallel_ = 0;
			std::move(dim_)._reset();
			std::move(epsilon_)._reset();
			std::move(rho_)._reset();
			std::move(rho_to_z_)._reset();
		}
		// cutoff of the power series expansion of conformal blocks at rho = 0
		[[nodiscard]] uint32_t n_Max() const { return n_Max_; }
		// passed to the constructor of ComplexFunction<Ring>
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		// the dimension of the spacetime
		[[nodiscard]] const mp::rational& dimension() const { return dim_; }
		// the dimension of the spacetime
		[[nodiscard]] uint32_t parallel() const { return parallel_; }
		// (dim - 2) / 2
		[[nodiscard]] const mp::rational& epsilon() const { return epsilon_; }
		// 3 - 2 sqrt(2)
		[[nodiscard]] const mp::real& rho() const { return rho_; }
		// convert a function of rho - (3 - 2 sqrt(2)) to a function of z - 1 / 2
		[[nodiscard]] const algebra::RealConverter& rho_to_z() const { return rho_to_z_; }
		[[nodiscard]] std::string str() const
		{
			std::ostringstream os;
			os << "Context(nMax=" << n_Max_ << ", lambda=" << lambda_ << ", dim=" << dim_ << ")";
			return os.str();
		}
		Context(uint32_t n_Max, uint32_t lambda, const mp::rational& dim, uint32_t p = 1);
		Context(Context&&) noexcept = delete;
		Context& operator=(Context&&) noexcept = delete;
		Context(const Context&) = delete;
		Context& operator=(const Context&) = delete;
		~Context() = default;
		[[nodiscard]] uint32_t _total_memory() const
		{
			return 1 + v_to_d_._total_memory() + gBlock_._total_memory() + rho_to_z_._total_memory();
		}
		[[nodiscard]] mp::rational unitarity_bound(uint32_t spin) const
		{
			return qboot::unitarity_bound(epsilon_, spin);
		}
		// calculate v ^ d gBlock and project to sym-symmetric part
		// F-type corresponds to sym = Odd, H-type to Even
		[[nodiscard]] algebra::ComplexFunction<mp::real> F_block(const mp::real& d,
		                                                         const algebra::ComplexFunction<mp::real>& gBlock,
		                                                         algebra::FunctionSymmetry sym) const
		{
			return mul_scalar(2, mul(v_to_d(d), gBlock).proj(sym));
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> evaluate(const ConformalBlock<PrimaryOperator>& block) const
		{
			assert(block.symmetry() != algebra::FunctionSymmetry::Mixed);
			return F_block(block.delta_half(), gBlock(block.get_op(), block.S(), block.P()), block.symmetry());
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> evaluate(const ConformalBlock<GeneralPrimaryOperator>& block,
		                                                          const mp::real& delta) const
		{
			assert(block.symmetry() != algebra::FunctionSymmetry::Mixed);
			return F_block(block.delta_half(), gBlock(block.get_op(delta), block.S(), block.P()), block.symmetry());
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> F_block(const PrimaryOperator& op, const mp::real& d1,
		                                                         const mp::real& d2, const mp::real& d3,
		                                                         const mp::real& d4,
		                                                         algebra::FunctionSymmetry sym) const
		{
			return F_block((d2 + d3) / 2, gBlock(op, _delta_S(d1, d2, d3, d4), _delta_P(d1, d2, d3, d4)), sym);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> F_block(const mp::real& d2, const mp::real& d3,
		                                                         const algebra::ComplexFunction<mp::real>& gBlock) const
		{
			return F_block((d2 + d3) / 2, gBlock);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> F_block(const mp::real& d,
		                                                         const algebra::ComplexFunction<mp::real>& gBlock) const
		{
			return F_block(d, gBlock, algebra::FunctionSymmetry::Odd);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> F_block(const PrimaryOperator& op, const mp::real& d1,
		                                                         const mp::real& d2, const mp::real& d3,
		                                                         const mp::real& d4) const
		{
			return F_block(op, d1, d2, d3, d4, algebra::FunctionSymmetry::Odd);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> H_block(const mp::real& d2, const mp::real& d3,
		                                                         const algebra::ComplexFunction<mp::real>& gBlock) const
		{
			return H_block((d2 + d3) / 2, gBlock);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> H_block(const mp::real& d,
		                                                         const algebra::ComplexFunction<mp::real>& gBlock) const
		{
			return F_block(d, gBlock, algebra::FunctionSymmetry::Even);
		}
		[[nodiscard]] algebra::ComplexFunction<mp::real> H_block(const PrimaryOperator& op, const mp::real& d1,
		                                                         const mp::real& d2, const mp::real& d3,
		                                                         const mp::real& d4) const
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
		algebra::ComplexFunction<mp::real> expand_off_diagonal(algebra::RealFunction<mp::real>&& realAxisResult,
		                                                       const PrimaryOperator& op, const mp::real& S,
		                                                       const mp::real& P) const;
	};
}  // namespace qboot

#endif  // QBOOT_CONTEXT_HPP_
