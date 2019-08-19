#ifndef REAL_FUNCTION_HPP_
#define REAL_FUNCTION_HPP_

#include <cstddef>      // for uint32_t
#include <ostream>      // for ostream
#include <type_traits>  // for true_type, false_type, enable_if, enable_if_t, is_same_v
#include <utility>      // for swap, move

#include "matrix.hpp"  // for Vector, Matrix, base_ring, is_intermediate_v, is_mpfr_real_v, is_iaddable_v, is_imultipliable_v, is_idividable_v, is_addable_v, is_subtractable_v, is_multipliable_v, union_ring_t
#include "real.hpp"    // for real, pow

namespace qboot
{
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class RealFunction;
	template <class T>
	struct is_realfunc;
	template <class T>
	inline constexpr bool is_realfunc_v = is_realfunc<T>::value;
	template <class R>
	struct is_realfunc<RealFunction<R>> : std::true_type
	{
	};
	template <class R>
	struct is_realfunc : std::false_type
	{
	};
	// real function of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->get(k) x ^ k + O(x ^ {lambda + 1})
	template <class Ring>
	class RealFunction
	{
		uint32_t lambda_;
		algebra::Vector<Ring> coeffs_;
		using base = typename algebra::base_ring<Ring>::type;
		using ring = Ring;
		using type = RealFunction;

	public:
		explicit RealFunction(uint32_t lambda) : lambda_(lambda), coeffs_(lambda + 1) {}
		template <class T, class = std::enable_if_t<algebra::is_intermediate_v<Ring, T> ||
		                                            (std::is_same_v<T, base> && !std::is_same_v<T, Ring>)>>
		explicit RealFunction(const RealFunction<T>& v) : lambda_(v.lambda_), coeffs_(v.coeffs_)
		{
		}
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		[[nodiscard]] uint32_t size() const { return coeffs_.size(); }
		void swap(RealFunction& other)
		{
			std::swap(lambda_, other.lambda_);
			coeffs_.swap(other.coeffs_);
		}
		void negate() { coeffs_.negate(); }
		[[nodiscard]] bool iszero() const { return coeffs_.iszero(); }
		[[nodiscard]] RealFunction clone() const
		{
			RealFunction f(lambda_);
			f.coeffs_ = coeffs_.clone();
			return f;
		}
		// get coefficient of the term x ^ k
		// 0 <= k <= lambda
		[[nodiscard]] Ring& get(uint32_t k) { return coeffs_[k]; }
		[[nodiscard]] const Ring& get(uint32_t k) const { return coeffs_[k]; }

		[[nodiscard]] algebra::Vector<Ring> as_vector() && { return std::move(coeffs_); }
		template <class = std::enable_if<algebra::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const { return coeffs_.norm(); }
		// multiply x ^ p
		void shift(uint32_t p)
		{
			for (uint32_t i = lambda_; i >= p; i--) get(i) = get(i - p);
			for (uint32_t i = 0; i < p; i++) get(i) = {};
		}

		template <class T, class = std::enable_if_t<algebra::is_iaddable_v<Ring, T>>>
		RealFunction& operator+=(const RealFunction<T>& v)
		{
			coeffs_ += v.coeffs_;
			return *this;
		}
		template <class T, class = std::enable_if_t<algebra::is_imultipliable_v<Ring, T>>>
		RealFunction& operator*=(const T& r)
		{
			coeffs_ *= r;
			return *this;
		}
		template <class T, class = std::enable_if_t<algebra::is_idividable_v<Ring, T>>>
		RealFunction& operator/=(const T& r)
		{
			coeffs_ /= r;
			return *this;
		}
		template <class R, class = std::enable_if_t<algebra::is_addable_v<Ring, R>>>
		friend RealFunction<algebra::union_ring_t<Ring, R>> operator+(const RealFunction& x, const RealFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			RealFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_);
			z.coeffs_ = x.coeffs_ + y.coeffs_;
			return z;
		}
		template <class R, class = std::enable_if_t<algebra::is_subtractable_v<Ring, R>>>
		friend RealFunction<algebra::union_ring_t<Ring, R>> operator-(const RealFunction& x, const RealFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			RealFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_);
			z.coeffs_ = x.coeffs_ - y.coeffs_;
			return z;
		}
		template <class R, class = std::enable_if_t<algebra::is_multipliable_v<Ring, R>>>
		friend RealFunction<algebra::union_ring_t<Ring, R>> operator*(const RealFunction& x, const RealFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			RealFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_);
			for (uint32_t k1 = 0; k1 <= x.lambda_; ++k1)
				for (uint32_t k2 = 0; k1 + k2 <= x.lambda_; ++k2) z.get(k1 + k2) += x.get(k1) * y.get(k2);
			return z;
		}
		template <class R, class = std::enable_if_t<!is_realfunc_v<R> && algebra::is_multipliable_v<Ring, R>>>
		friend RealFunction<algebra::union_ring_t<Ring, R>> operator*(const RealFunction& x, const R& r)
		{
			RealFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_);
			z.coeffs_ = x.coeffs_ * r;
			return z;
		}
		template <class R, class = std::enable_if_t<!is_realfunc_v<R> && algebra::is_multipliable_v<Ring, R>>>
		friend RealFunction<algebra::union_ring_t<Ring, R>> operator*(const R& r, const RealFunction& x)
		{
			return x * r;
		}
		friend std::ostream& operator<<(std::ostream& out, const RealFunction& v)
		{
			auto f = false;
			for (uint32_t k = 0; k <= v.lambda_; ++k)
			{
				if (f) out << " + ";
				out << "(" << v.get(k) << ")";
				if (k > 0) out << " * x";
				if (k > 1) out << " ^ " << k;
				f = true;
			}
			return out;
		}
		template <class Ring2>
		friend class RealFunction;
		template <class Ring2>
		friend class RealConverter;
		template <class Ring2>
		friend class RealFunctionWithPower;
	};
	// f(x) = (a + b * x) ^ p
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	RealFunction<Ring> power_function(const Ring& a, const Ring& b, const Ring& p, uint32_t lambda)
	{
		assert(a != 0);
		RealFunction<Ring> f(lambda);
		f.get(0) = mpfr::pow(a, p);
		Ring tmp = b / a;
		for (uint32_t k = 1; k <= lambda; ++k) { f.get(k) = f.get(k - 1) * tmp * (p - (k - 1)) / k; }
		return f;
	}

	// convert a function f(x) of x to a function g(y) of y
	// Let x = x(y), then g(y) = f(x(y))
	// x(0) must be 0
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class RealConverter
	{
		uint32_t lambda_;
		algebra::Matrix<Ring> mat_;

	public:
		// x = func(y)
		// this converts a functino of x to a function of y
		RealConverter(const RealFunction<Ring>& func)
		    : lambda_(func.lambda()), mat_(func.lambda() + 1, func.lambda() + 1)
		{
			// assert(func.get(0).iszero());
			mat_.get(0, 0) = 1;
			auto pf = func.clone();
			for (uint32_t n = 1; n <= lambda_; ++n)
			{
				for (uint32_t i = n; i <= lambda_; ++i) mat_.get(n, i) = pf.get(i);
				pf = pf * func;
			}
		}
		RealConverter inverse() const
		{
			auto invmat = mat_.clone();
			invmat.transpose();
			invmat = invmat.lower_triangular_inverse();
			RealFunction<Ring> f(lambda_);
			for (uint32_t k = 1; k <= lambda_; ++k) f.get(k) = invmat.get(k, 1);
			return RealConverter(f);
		}
		// convert a function f of x to a function of y where x = func(y)
		template <class R>
		RealFunction<R> convert(const RealFunction<R>& f) const
		{
			assert(lambda_ == f.lambda_);
			RealFunction<R> g(lambda_);
			g.coeffs_ = f.coeffs_ * mat_;
			return g;
		}
		const algebra::Matrix<Ring>& matrix() const { return mat_; }
	};
	// real function multiplied by x ^ {pow} of x at x = 0
	// take derivatives (der x) ^ k upto k <= lambda
	// namely, a function is represented as
	//   \sum_{k = 0}^{lambda} this->get(k) x ^ {k + pow} + O(x ^ {pow + lambda + 1})
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class RealFunctionWithPower
	{
		RealFunction<Ring> f_;
		Ring pow_;

	public:
		RealFunctionWithPower(const RealFunction<Ring>& f, const Ring& p) : f_(f.clone()), pow_(p) {}
		const RealFunction<Ring>& func() const { return f_; }
		const Ring& power() const { return pow_; }
		// take derivative
		void derivate()
		{
			for (uint32_t k = 0; k <= f_.lambda_; ++k) f_.get(k) *= pow_ + k;
			pow_ -= 1;
		}
		// evaluate at x = x
		template <class R, class = std::enable_if_t<algebra::is_imultipliable_v<Ring, R>>>
		[[nodiscard]] Ring eval(const R& x) const
		{
			Ring s{};
			for (uint32_t i = f_.lambda_; i <= f_.lambda_; --i)
			{
				s *= x;
				s += f_.get(i);
			}
			s *= mpfr::pow(x, pow_);
			return s;
		}
	};
}  // namespace qboot

#endif  // REAL_FUNCTION_HPP_
