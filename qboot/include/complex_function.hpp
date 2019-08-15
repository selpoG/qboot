#ifndef COMPLEX_FUNCTION_HPP_
#define COMPLEX_FUNCTION_HPP_

#include <cstddef>  // for uint32_t

#include "matrix.hpp"
#include "real.hpp"

namespace qboot
{
	// holomorphic function of z = x + sqrt(y) (z^* = x - sqrt(y))
	// take derivatives (der x) ^ m (der y) ^ n upto m + 2 * n <= lambda
	// namely, a function is represented as
	//   \sum_{n = 0}^{lambda / 2}\sum_{m = 0}^{lambda - 2 * n} this->get(n, m) (x - x0)^m (y - y0)^n + ...
	// we take (x0, y0) = (1 / 2, 0)
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class ComplexFunction
	{
		uint32_t lambda_;
		algebra::Vector<Ring> coeffs_;

	public:
		explicit ComplexFunction(uint32_t lambda) : lambda_(lambda), coeffs_((lambda + 2) * (lambda + 2) / 4) {}
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		void swap(ComplexFunction& other)
		{
			std::swap(lambda_, other.lambda_);
			coeffs_.swap(other.coeffs_);
		}
		void negate() { coeffs_.negate(); }
		[[nodiscard]] bool iszero() const { return coeffs_.iszero(); }
		[[nodiscard]] ComplexFunction clone() const
		{
			ComplexFunction f(lambda_);
			f.coeffs_ = coeffs_.clone();
			return f;
		}
		// get coefficient of the term (x - x0) ^ {dx} (y - y0) ^ {dy}
		[[nodiscard]] Ring& get(int32_t dx, int32_t dy) { return get(uint32_t(dx), uint32_t(dy)); }
		[[nodiscard]] const Ring& get(int32_t dx, int32_t dy) const { return get(uint32_t(dx), uint32_t(dy)); }
		[[nodiscard]] Ring& get(uint32_t dx, uint32_t dy) { return coeffs_[(lambda_ + 2 - dy) * dy + dx]; }
		[[nodiscard]] const Ring& get(uint32_t dx, uint32_t dy) const { return coeffs_[(lambda_ + 2 - dy) * dy + dx]; }
		[[nodiscard]] algebra::Vector<Ring> as_vector() && { return std::move(coeffs_); }
		template <class = std::enable_if<algebra::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const { return coeffs_.norm(); }
		template <class R, class = std::enable_if_t<algebra::is_subtractable_v<Ring, R>>>
		friend ComplexFunction<algebra::union_ring_t<Ring, R>> operator-(const ComplexFunction& x,
		                                                                 const ComplexFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_);
			z.coeffs_ = x.coeffs_ - y.coeffs_;
			return z;
		}
		friend std::ostream& operator<<(std::ostream& out, const ComplexFunction& v)
		{
			out << "[";
			auto f = false;
			for (uint32_t dy = 0; dy <= v.lambda_ / 2; ++dy)
			{
				if (f) out << ", ";
				auto g = false;
				out << "[";
				for (uint32_t dx = 0; dx + 2 * dy <= v.lambda_; ++dx)
				{
					if (g) out << ", ";
					out << v.get(dx, dy);
					g = true;
				}
				out << "]";
				f = true;
			}
			return out << "]";
		}
	};
}  // namespace qboot

#endif  // COMPLEX_FUNCTION_HPP_
