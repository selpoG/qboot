#ifndef QBOOT_REAL_FUNCTION_HPP_
#define QBOOT_REAL_FUNCTION_HPP_

#include <cassert>  // for assert
#include <cstddef>  // for uint32_t
#include <ostream>  // for ostream
#include <utility>  // for swap, move

#include "matrix.hpp"      // for Vector, Matrix
#include "polynomial.hpp"  // for polynomialize_t, to_pol
#include "real.hpp"        // for real, pow

namespace algebra
{
	// real function of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->at(k) x ^ k + O(x ^ {lambda + 1})
	template <class Ring>
	class RealFunction
	{
		template <class Ring2>
		friend class RealFunction;
		uint32_t lambda_;
		Vector<Ring> coeffs_;
		explicit RealFunction(Vector<Ring>&& vec) : lambda_(vec.size() - 1), coeffs_(std::move(vec)) {}

	public:
		template <class Ring2>
		friend class RealConverter;
		RealFunction() : RealFunction(0) {}
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
		[[nodiscard]] Ring& at(uint32_t k) & { return coeffs_[k]; }
		[[nodiscard]] const Ring& at(uint32_t k) const& { return coeffs_[k]; }

		[[nodiscard]] auto abs() const { return coeffs_.abs(); }
		[[nodiscard]] auto norm() const { return coeffs_.norm(); }
		// multiply x ^ p
		void shift(uint32_t p) &
		{
			for (uint32_t i = lambda_; i >= p; i--) at(i) = at(i - p);
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
		template <class T>
		RealFunction& operator*=(const T& r) &
		{
			coeffs_ *= r;
			return *this;
		}
		template <class T>
		RealFunction& operator/=(const T& r) &
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
		friend RealFunction operator+(RealFunction&& x, RealFunction&& y) { return std::move(x += y); }
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
		friend RealFunction operator-(RealFunction&& x, RealFunction&& y) { return std::move(x -= y); }
		friend RealFunction mul(const RealFunction& x, const RealFunction& y)
		{
			assert(x.lambda_ == y.lambda_);
			RealFunction z(x.lambda_);
			for (uint32_t k1 = 0; k1 <= x.lambda_; ++k1)
				for (uint32_t k2 = 0; k1 + k2 <= x.lambda_; ++k2) z.at(k1 + k2) += mul(x.at(k1), y.at(k2));
			return z;
		}
		template <class R>
		friend RealFunction mul_scalar(const R& r, const RealFunction& x)
		{
			return RealFunction(mul_scalar(r, x.coeffs_));
		}
		template <class R>
		friend RealFunction mul_scalar(const R& r, RealFunction&& x)
		{
			return std::move(x *= r);
		}
		template <class R>
		friend RealFunction operator/(const RealFunction& x, const R& r)
		{
			return RealFunction(x.coeffs_ / r);
		}
		template <class R>
		friend RealFunction operator/(RealFunction&& x, const R& r)
		{
			return std::move(x /= r);
		}
		friend bool operator==(const RealFunction& x, const RealFunction& y)
		{
			return x.lambda_ == y.lambda_ && x.coeffs_ == y.coeffs_;
		}
		friend bool operator!=(const RealFunction& x, const RealFunction& y) { return !(x == y); }
		template <class Real>
		[[nodiscard]] RealFunction<evaluated_t<Ring>> eval(const Real& x) const
		{
			return RealFunction<evaluated_t<Ring>>(coeffs_.eval(x), lambda_);
		}
	};
	template <class Ring>
	RealFunction<polynomialize_t<Ring>> to_pol(Vector<RealFunction<Ring>>& coeffs)
	{
		uint32_t lambda = coeffs[0].lambda(), len = coeffs.size();
		RealFunction<polynomialize_t<Ring>> ans(lambda);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r <= lambda; ++r)
		{
			for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs[i].at(r));
			ans.at(r) = to_pol(v);
		}
		return ans;
	}
	template <class R>
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
	template <class Ring>
	RealFunction<Ring> power_function(const Ring& a, const Ring& b, const Ring& p, uint32_t lambda)
	{
		assert(a != 0);
		RealFunction<Ring> f(lambda);
		f.at(0) = mpfr::pow(a, p);
		Ring tmp = b / a;
		for (uint32_t k = 1; k <= lambda; ++k) { f.at(k) = f.at(k - 1) * tmp * (p - (k - 1)) / k; }
		return f;
	}

	// convert a function f(x) of x to a function g(y) of y
	// Let x = x(y), then g(y) = f(x(y))
	// x(0) must be 0
	template <class Ring>
	class RealConverter
	{
		uint32_t lambda_;
		Matrix<Ring> mat_;

	public:
		// x = func(y)
		// this converts a functino of x to a function of y
		explicit RealConverter(const RealFunction<Ring>& func)
		    : lambda_(func.lambda()), mat_(func.lambda() + 1, func.lambda() + 1)
		{
			// assert(func.at(0).iszero());
			mat_.at(0, 0) = 1;
			auto pf = func.clone();
			for (uint32_t n = 1; n <= lambda_; ++n)
			{
				for (uint32_t i = n; i <= lambda_; ++i) mat_.at(n, i) = pf.at(i);
				pf = mul(pf, func);
			}
		}
		[[nodiscard]] RealConverter inverse() const
		{
			auto invmat = mat_.clone();
			invmat.transpose();
			invmat = invmat.lower_triangular_inverse();
			RealFunction<Ring> f(lambda_);
			for (uint32_t k = 1; k <= lambda_; ++k) f.at(k) = invmat.at(k, 1);
			return RealConverter(f);
		}
		// convert a function f of x to a function of y where x = func(y)
		template <class R>
		[[nodiscard]] RealFunction<R> convert(const RealFunction<R>& f) const
		{
			assert(lambda_ == f.lambda());
			return RealFunction<R>(dot(f.coeffs_, mat_));
		}
	};
	// real function multiplied by x ^ {pow} of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->at(k) x ^ {k + pow} + O(x ^ {pow + lambda + 1})
	template <class Ring>
	class RealFunctionWithPower
	{
		RealFunction<Ring> f_;
		Ring pow_;

	public:
		RealFunctionWithPower(const RealFunction<Ring>& f, const Ring& p) : f_(f.clone()), pow_(p) {}
		RealFunctionWithPower(RealFunction<Ring>&& f, const Ring& p) : f_(std::move(f)), pow_(p) {}
		[[nodiscard]] const RealFunction<Ring>& func() const { return f_; }
		[[nodiscard]] const Ring& power() const { return pow_; }
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
		[[nodiscard]] Ring approximate(const R& x) const
		{
			Ring s{};
			for (uint32_t i = f_.lambda(); i <= f_.lambda(); --i)
			{
				s *= x;
				s += f_.at(i);
			}
			s *= mpfr::pow(x, pow_);
			return s;
		}
	};
}  // namespace algebra

#endif  // QBOOT_REAL_FUNCTION_HPP_
