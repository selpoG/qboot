#ifndef COMPLEX_FUNCTION_HPP_
#define COMPLEX_FUNCTION_HPP_

#include <cstddef>  // for uint32_t

#include "matrix.hpp"
#include "real.hpp"

namespace qboot
{
	// a complex function is even, odd or mixed under z -> 1 - z
	// Mixed means the function may contain both Even and Odd parts.
	enum class FunctionSymmetry : uint32_t
	{
		Even = 0,
		Odd = 1,
		Mixed
	};
	inline constexpr bool matches(FunctionSymmetry s, uint32_t dx) noexcept
	{
		return s == FunctionSymmetry::Mixed || (dx % 2) == uint32_t(s);
	}
	// times(Mixed, _) = times(_, Mixed) = Mixed
	// times(s, s) = Even, times(s, !s) = Odd
	inline constexpr FunctionSymmetry times(FunctionSymmetry s, FunctionSymmetry t) noexcept
	{
		return s == FunctionSymmetry::Mixed || t == FunctionSymmetry::Mixed
		           ? FunctionSymmetry::Mixed
		           : s == t ? FunctionSymmetry::Even : FunctionSymmetry::Odd;
	}
	// plus(s, s) = s
	// plus(Mixed, *) = plus(*, Mixed) = plus(s, !s) = Mixed
	inline constexpr FunctionSymmetry plus(FunctionSymmetry s, FunctionSymmetry t) noexcept
	{
		return s == t ? s : FunctionSymmetry::Mixed;
	}
	inline constexpr uint32_t _triangle_num(uint32_t n) noexcept { return n * (n + 1) / 2; }
	// #{(m, n) | 0 <= m, n and m + 2 n <= lambda and m mod 2 == sym}
	inline constexpr uint32_t function_dimension(uint32_t lambda,
	                                             FunctionSymmetry sym = FunctionSymmetry::Mixed) noexcept
	{
		switch (sym)
		{
		case FunctionSymmetry::Even: return _triangle_num((lambda + 2) / 2);
		case FunctionSymmetry::Odd: return _triangle_num((lambda + 1) / 2);
		case FunctionSymmetry::Mixed:
		default: return (lambda + 2) * (lambda + 2) / 4;  // equals to dim(lambda, Even) + dim(lambda, Odd)
		}
	}
	// complex function of z = x + sqrt(y) (z^* = x - sqrt(y))
	// take derivatives (der x) ^ m (der y) ^ n upto m + 2 * n <= lambda
	// namely, a function is represented as
	//   \sum_{n = 0}^{lambda / 2}\sum_{m = 0}^{lambda - 2 * n} this->get(n, m) (x - x0) ^ m (y - y0) ^ n + ...
	// we take (x0, y0) = (1 / 2, 0)
	// if a nontrivial symmetry even (resp. odd) is given, m runs over even (resp. odd) number only.
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class ComplexFunction
	{
		FunctionSymmetry sym_ = FunctionSymmetry::Mixed;
		uint32_t lambda_;
		algebra::Vector<Ring> coeffs_;
		// indices are aligned by lexicographical order of (dy, dx)
		// i.e., (dx, dy) is the index(dx, dy)-th element in
		//   [(dx, dy) for dy in range(lambda // 2) for dx in range(lambda - 2 * dy) if sym is Mixed or dx % 2 == sym]
		uint32_t index(uint32_t dx, uint32_t dy) const noexcept
		{
			switch (sym_)
			{
			case FunctionSymmetry::Even: return (lambda_ / 2 * 2 + 3 - dy) * dy / 2 + dx / 2;
			case FunctionSymmetry::Odd: return ((lambda_ - 1) / 2 * 2 + 3 - dy) * dy / 2 + (dx - 1) / 2;
			case FunctionSymmetry::Mixed:
			default: return (lambda_ + 2 - dy) * dy + dx;
			}
		}

	public:
		explicit ComplexFunction(uint32_t lambda, FunctionSymmetry sym = FunctionSymmetry::Mixed)
		    : sym_(sym), lambda_(lambda), coeffs_(function_dimension(lambda, sym))
		{
		}
		[[nodiscard]] FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		[[nodiscard]] uint32_t size() const { return coeffs_.size(); }
		void swap(ComplexFunction& other)
		{
			std::swap(sym_, other.sym_);
			std::swap(lambda_, other.lambda_);
			coeffs_.swap(other.coeffs_);
		}
		void negate() { coeffs_.negate(); }
		[[nodiscard]] bool iszero() const { return coeffs_.iszero(); }
		[[nodiscard]] ComplexFunction clone() const
		{
			ComplexFunction f(lambda_, sym_);
			f.coeffs_ = coeffs_.clone();
			return f;
		}
		// get coefficient of the term (x - x0) ^ {dx} (y - y0) ^ {dy}
		// 0 <= dx, dy and dx + 2 dy <= lambda
		// if symmetry is even or odd, the parity of dx must equals symmetry
		[[nodiscard]] Ring& get(int32_t dx, int32_t dy) { return get(uint32_t(dx), uint32_t(dy)); }
		[[nodiscard]] const Ring& get(int32_t dx, int32_t dy) const { return get(uint32_t(dx), uint32_t(dy)); }
		[[nodiscard]] Ring& get(uint32_t dx, uint32_t dy) { return coeffs_[index(dx, dy)]; }
		[[nodiscard]] const Ring& get(uint32_t dx, uint32_t dy) const { return coeffs_[index(dx, dy)]; }
		[[nodiscard]] Ring move_get(uint32_t dx, uint32_t dy) && { return std::move(coeffs_[index(dx, dy)]); }

		[[nodiscard]] algebra::Vector<Ring> as_vector() && { return std::move(coeffs_); }
		// z -> 1 - z
		void flip()
		{
			switch (sym_)
			{
			case qboot::FunctionSymmetry::Even: break;
			case qboot::FunctionSymmetry::Odd: coeffs_.negate(); break;
			case qboot::FunctionSymmetry::Mixed:
			default:
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 1; dx + 2 * dy <= lambda_; dx += 2) get(dx, dy).negate();
				break;
			}
		}
		template <class = std::enable_if<algebra::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const { return coeffs_.norm(); }

		[[nodiscard]] ComplexFunction proj(FunctionSymmetry sym) &&
		{
			if (sym_ == sym) return std::move(*this);
			ComplexFunction z(lambda_, sym);
			if (sym_ == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym, dx)) z.get(dx, dy) = std::move(*this).move_get(dx, dy);
			}
			if (sym == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym_, dx)) z.get(dx, dy) = std::move(*this).move_get(dx, dy);
			}
			return z;
		}
		// project this function to sym-symmetric part
		[[nodiscard]] ComplexFunction proj(FunctionSymmetry sym) const&
		{
			ComplexFunction z(lambda_, sym);
			if (sym_ == sym) z.coeffs_ = coeffs_.clone();
			if (sym_ == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym, dx)) z.get(dx, dy) = get(dx, dy);
			}
			if (sym == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym_, dx)) z.get(dx, dy) = get(dx, dy);
			}
			// otherwise (even to odd or odd to even), proj is vanishing
			return z;
		}

		template <class R, class = std::enable_if_t<algebra::is_addable_v<Ring, R>>>
		friend ComplexFunction<algebra::union_ring_t<Ring, R>> operator+(const ComplexFunction& x,
		                                                                 const ComplexFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_, plus(x.sym_, y.sym_));
			if (x.sym_ == y.sym_)
				z.coeffs_ = x.coeffs_ + y.coeffs_;
			else
				for (uint32_t dy = 0; dy <= x.lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= x.lambda_; ++dx)
					{
						if (matches(x.sym_, dx)) z.get(dx, dy) += x.get(dx, dy);
						if (matches(y.sym_, dx)) z.get(dx, dy) += y.get(dx, dy);
					}
			return z;
		}
		template <class R, class = std::enable_if_t<algebra::is_subtractable_v<Ring, R>>>
		friend ComplexFunction<algebra::union_ring_t<Ring, R>> operator-(const ComplexFunction& x,
		                                                                 const ComplexFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_, plus(x.sym_, y.sym_));
			if (x.sym_ == y.sym_)
				z.coeffs_ = x.coeffs_ - y.coeffs_;
			else
				for (uint32_t dy = 0; dy <= x.lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= x.lambda_; ++dx)
					{
						if (matches(x.sym_, dx)) z.get(dx, dy) += x.get(dx, dy);
						if (matches(y.sym_, dx)) z.get(dx, dy) -= y.get(dx, dy);
					}
			return z;
		}
		template <class R, class = std::enable_if_t<algebra::is_multipliable_v<Ring, R>>>
		friend ComplexFunction<algebra::union_ring_t<Ring, R>> operator*(const ComplexFunction& x,
		                                                                 const ComplexFunction<R>& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction<algebra::union_ring_t<Ring, R>> z(x.lambda_, times(x.sym_, y.sym_));
			for (uint32_t dy1 = 0; dy1 <= x.lambda_ / 2; ++dy1)
				for (uint32_t dx1 = 0; dx1 + 2 * dy1 <= x.lambda_; ++dx1)
				{
					if (!matches(x.sym_, dx1)) continue;
					for (uint32_t dy2 = 0; dy1 + dy2 <= x.lambda_ / 2; ++dy2)
						for (uint32_t dx2 = 0; dx1 + dx2 + 2 * dy1 + 2 * dy2 <= x.lambda_; ++dx2)
						{
							if (!matches(y.sym_, dx2)) continue;
							z.get(dx1 + dx2, dy1 + dy2) += x.get(dx1, dy1) * y.get(dx2, dy2);
						}
				}
			return z;
		}
		friend std::ostream& operator<<(std::ostream& out, const ComplexFunction& v)
		{
			out << "[";
			auto f = false;
			for (uint32_t dy = 0; dy <= v.lambda_ / 2; ++dy)
			{
				// do not show an empty array []
				if (v.sym_ == FunctionSymmetry::Odd && 2 * dy == v.lambda_) continue;
				if (f) out << ", ";
				auto g = false;
				out << "[";
				for (uint32_t dx = 0; dx + 2 * dy <= v.lambda_; ++dx)
				{
					if (!matches(v.sym_, dx)) continue;
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
	// v ^ d = ((1 - z) (1 - z_bar)) ^ d as a function of x, y
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	ComplexFunction<Real> v_to_d(const Real& d, uint32_t lambda)
	{
		// f.get(m, n) = (der x) ^ m (der y) ^ n ((x - 1) ^ 2 - y) ^ d / (n! m!)
		//             = (-1) ^ {n + m} 2 ^ {2 n + m} lf(d, n) lf(2 (d - n), m) / (n! m! 4 ^ d)
		// where lf(x, n) = x * (x - 1) * ... * (x - (n - 1)) (falling factorial)
		ComplexFunction f(lambda);
		f.get(0, 0) = mpfr::pow(Real(0.25), d);
		for (uint32_t n = 0; n <= lambda / 2; ++n)
		{
			if (n > 0) f.get(0u, n) = f.get(0u, n - 1) * 4 * (-d + (n - 1)) / n;
			for (uint32_t m = 1; m + 2 * n <= lambda; ++m)
				f.get(m, n) = f.get(m - 1, n) * (-4 * d + 2 * (m + 2 * n - 1)) / m;
		}
		return f;
	}
}  // namespace qboot

#endif  // COMPLEX_FUNCTION_HPP_
