#ifndef QBOOT_ALGEBRA_REAL_FUNCTION_HPP_
#define QBOOT_ALGEBRA_REAL_FUNCTION_HPP_

#include <cassert>  // for assert
#include <cstddef>  // for uint32_t
#include <ostream>  // for ostream
#include <utility>  // for swap, move

#include "qboot/algebra/matrix.hpp"      // for Vector, Matrix
#include "qboot/algebra/polynomial.hpp"  // for _polynomialize_t, to_pol
#include "qboot/mp/real.hpp"             // for real, pow

namespace qboot::algebra
{
	// real function of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->at(k) x ^ k + O(x ^ {lambda + 1})
	// R must be mp::real or Polynomial
	template <Ring R>
	class RealFunction
	{
		template <Ring R2>
		friend class RealFunction;
		uint32_t lambda_;
		Vector<R> coeffs_;
		explicit RealFunction(Vector<R>&& vec) : lambda_(vec.size() - 1), coeffs_(std::move(vec)) {}

	public:
		friend class RealConverter;
		void _reset() &&
		{
			lambda_ = 0;
			std::move(coeffs_)._reset();
		}
		// make an uninitialized object
		RealFunction() : lambda_(0), coeffs_{} {}
		explicit RealFunction(uint32_t lambda) : lambda_(lambda), coeffs_(lambda + 1) {}
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		[[nodiscard]] uint32_t size() const { return coeffs_.size(); }
		void swap(RealFunction& other) &
		{
			coeffs_.swap(other.coeffs_);
			std::swap(lambda_, other.lambda_);
		}
		void negate() & { coeffs_.negate(); }
		[[nodiscard]] bool iszero() const { return coeffs_.iszero(); }
		[[nodiscard]] RealFunction clone() const { return RealFunction(coeffs_.clone()); }
		// get coefficient of the term x ^ k
		// 0 <= k <= lambda
		[[nodiscard]] R& at(uint32_t k) & { return coeffs_[k]; }
		[[nodiscard]] const R& at(uint32_t k) const& { return coeffs_[k]; }

		[[nodiscard]] auto abs() const { return coeffs_.abs(); }
		[[nodiscard]] auto norm() const { return coeffs_.norm(); }
		// multiply x ^ p
		void shift(uint32_t p) &
		{
			for (uint32_t i = lambda_; i >= p; i--) at(i) = std::move(at(i - p));
			for (uint32_t i = 0; i < p; ++i) at(i) = {};
		}

		RealFunction& operator+=(const RealFunction& v) &
		{
			coeffs_ += v.coeffs_;
			return *this;
		}
		RealFunction& operator-=(const RealFunction& v) &
		{
			coeffs_ -= v.coeffs_;
			return *this;
		}
		template <class S>
		RealFunction& operator*=(const S& r) &
		{
			coeffs_ *= r;
			return *this;
		}
		template <class S>
		RealFunction& operator/=(const S& r) &
		{
			coeffs_ /= r;
			return *this;
		}
		RealFunction operator+() const& { return clone(); }
		RealFunction operator+() && { return std::move(*this); }
		RealFunction operator-() const& { return RealFunction(-coeffs_); }
		RealFunction operator-() &&
		{
			negate();
			return std::move(*this);
		}
		friend RealFunction operator+(const RealFunction& x, const RealFunction& y)
		{
			return RealFunction(x.coeffs_ + y.coeffs_);
		}
		friend RealFunction operator+(RealFunction&& x, const RealFunction& y) { return std::move(x += y); }
		friend RealFunction operator+(const RealFunction& x, RealFunction&& y) { return std::move(y += x); }
		friend RealFunction operator+(RealFunction&& x, RealFunction&& y)
		{
			x += y;
			std::move(y)._reset();
			return std::move(x);
		}
		friend RealFunction operator-(const RealFunction& x, const RealFunction& y)
		{
			return RealFunction(x.coeffs_ - y.coeffs_);
		}
		friend RealFunction operator-(RealFunction&& x, const RealFunction& y) { return std::move(x -= y); }
		friend RealFunction operator-(const RealFunction& x, RealFunction&& y)
		{
			y.coeffs_ = x.coeffs_ - std::move(y.coeffs_);
			return std::move(y);
		}
		friend RealFunction operator-(RealFunction&& x, RealFunction&& y)
		{
			x -= y;
			std::move(y)._reset();
			return std::move(x);
		}
		friend RealFunction mul(const RealFunction& x, const RealFunction& y)
		{
			assert(x.lambda_ == y.lambda_);
			RealFunction z(x.lambda_);
			for (uint32_t k1 = 0; k1 <= x.lambda_; ++k1)
				for (uint32_t k2 = 0; k1 + k2 <= x.lambda_; ++k2) z.at(k1 + k2) += mul(x.at(k1), y.at(k2));
			return z;
		}
		template <class S>
		friend RealFunction mul_scalar(const S& r, const RealFunction& x)
		{
			return RealFunction(mul_scalar(r, x.coeffs_));
		}
		template <class S>
		friend RealFunction mul_scalar(const S& r, RealFunction&& x)
		{
			return std::move(x *= r);
		}
		template <class S>
		friend RealFunction operator/(const RealFunction& x, const S& r)
		{
			return RealFunction(x.coeffs_ / r);
		}
		template <class S>
		friend RealFunction operator/(RealFunction&& x, const S& r)
		{
			return std::move(x /= r);
		}
		friend bool operator==(const RealFunction& x, const RealFunction& y)
		{
			return x.lambda_ == y.lambda_ && x.coeffs_ == y.coeffs_;
		}
		friend bool operator!=(const RealFunction& x, const RealFunction& y) { return !(x == y); }
		[[nodiscard]] RealFunction<_evaluated_t<R>> eval(const mp::real& x) const
		{
			return RealFunction<_evaluated_t<R>>(coeffs_.eval(x));
		}
	};
	template <Ring R>
	RealFunction<_polynomialize_t<R>> to_pol(Vector<RealFunction<R>>* coeffs)
	{
		uint32_t lambda = coeffs->at(0).lambda(), len = coeffs->size();
		RealFunction<_polynomialize_t<R>> ans(lambda);
		Vector<R> v(len);
		for (uint32_t r = 0; r <= lambda; ++r)
		{
			for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs->at(i).at(r));
			ans.at(r) = to_pol(&v);
		}
		return ans;
	}
	template <Ring R>
	std::ostream& operator<<(std::ostream& out, const RealFunction<R>& v)
	{
		auto f = false;
		for (uint32_t k = 0; k <= v.lambda(); ++k)
		{
			if (f) out << " + ";
			out << "(" << v.at(k) << ")";
			if (k > 0) out << " * x";
			if (k > 1) out << " ^ " << k;
			f = true;
		}
		return out;
	}
	// f(x) = (a + b * x) ^ p
	RealFunction<mp::real> power_function(const mp::real& a, const mp::real& b, const mp::real& p, uint32_t lambda);

	// convert a function f(x) of x to a function g(y) of y
	// Let x = x(y), then g(y) = f(x(y))
	// x(0) must be 0
	class RealConverter
	{
		uint32_t lambda_;
		Matrix<mp::real> mat_;

	public:
		void _reset() &&
		{
			lambda_ = 0;
			std::move(mat_)._reset();
		}
		// x = func(y)
		// this converts a functino of x to a function of y
		explicit RealConverter(const RealFunction<mp::real>& func);
		[[nodiscard]] RealConverter inverse() const;
		// convert a function f of x to a function of y where x = func(y)
		template <Ring R>
		[[nodiscard]] RealFunction<R> convert(const RealFunction<R>& f) const
		{
			assert(lambda_ == f.lambda());
			return RealFunction<R>(dot(f.coeffs_, mat_));
		}
		[[nodiscard]] uint32_t _total_memory() const { return mat_.row() * mat_.column(); }
	};
	// real function multiplied by x ^ {pow} of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->at(k) x ^ {k + pow} + O(x ^ {pow + lambda + 1})
	class RealFunctionWithPower
	{
		RealFunction<mp::real> f_;
		mp::real pow_;

	public:
		void _reset() &&
		{
			std::move(f_)._reset();
			std::move(pow_)._reset();
		}
		RealFunctionWithPower(const RealFunction<mp::real>& f, const mp::real& p) : f_(f.clone()), pow_(p) {}
		RealFunctionWithPower(RealFunction<mp::real>&& f, const mp::real& p) : f_(std::move(f)), pow_(p) {}
		[[nodiscard]] const RealFunction<mp::real>& func() const { return f_; }
		[[nodiscard]] const mp::real& power() const { return pow_; }
		// take derivative
		void derivate() &
		{
			for (uint32_t k = 0; k <= f_.lambda(); ++k) f_.at(k) *= pow_ + k;
			pow_ -= 1;
		}
		void swap(RealFunctionWithPower& f) &
		{
			f_.swap(f.f_);
			pow_.swap(f.pow_);
		}
		[[nodiscard]] RealFunctionWithPower clone() const { return RealFunctionWithPower(f_, pow_); }
		// evaluate at x = x
		template <class R>
		[[nodiscard]] mp::real approximate(const R& x) const
		{
			mp::real s{};
			for (uint32_t i = f_.lambda(); i <= f_.lambda(); --i)
			{
				if constexpr (std::same_as<R, mp::real>)
					mp::fma(s, s, x, f_.at(i));
				else
				{
					s *= x;
					s += f_.at(i);
				}
			}
			s *= mp::pow(x, pow_);
			return s;
		}
	};
}  // namespace qboot::algebra

#endif  // QBOOT_ALGEBRA_REAL_FUNCTION_HPP_
