#ifndef QBOOT_COMPLEX_FUNCTION_HPP_
#define QBOOT_COMPLEX_FUNCTION_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t
#include <ostream>  // for ostream
#include <utility>  // for move, swap

#include "matrix.hpp"      // for Vector
#include "polynomial.hpp"  // for polynomialize_t, to_pol
#include "real.hpp"        // for real, pow

namespace algebra
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
	// take derivatives (der x) ^ m (der y) ^ n upto m + 2 n <= lambda
	// namely, a function is represented as
	//   \sum_{n = 0}^{lambda / 2} \sum_{m = 0}^{lambda - 2 n} this->at(n, m) (x - x0) ^ m (y - y0) ^ n + ...
	// we take (x0, y0) = (1 / 2, 0)
	// if a nontrivial symmetry even (resp. odd) is given, m runs over even (resp. odd) number only.
	template <class Ring>
	class ComplexFunction
	{
		template <class Ring2>
		friend class ComplexFunction;
		FunctionSymmetry sym_ = FunctionSymmetry::Mixed;
		uint32_t lambda_;
		Vector<Ring> coeffs_;
		// indices are aligned by lexicographical order of (dy, dx)
		// i.e., (dx, dy) is the index(dx, dy)-th element in
		//   [(dx, dy) for dy in range(lambda // 2) for dx in range(lambda - 2 * dy) if sym is Mixed or dx % 2 == sym]
		[[nodiscard]] uint32_t index(uint32_t dx, uint32_t dy) const noexcept
		{
			switch (sym_)
			{
			case FunctionSymmetry::Even: return (lambda_ / 2 * 2 + 3 - dy) * dy / 2 + dx / 2;
			case FunctionSymmetry::Odd: return ((lambda_ - 1) / 2 * 2 + 3 - dy) * dy / 2 + (dx - 1) / 2;
			case FunctionSymmetry::Mixed:
			default: return (lambda_ + 2 - dy) * dy + dx;
			}
		}
		ComplexFunction(Vector<Ring>&& v, uint32_t l, FunctionSymmetry sym)
		    : sym_(sym), lambda_(l), coeffs_(std::move(v))
		{
		}

	public:
		ComplexFunction() : ComplexFunction(0) {}
		explicit ComplexFunction(uint32_t lambda, FunctionSymmetry sym = FunctionSymmetry::Mixed)
		    : sym_(sym), lambda_(lambda), coeffs_(function_dimension(lambda, sym))
		{
		}
		[[nodiscard]] FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] uint32_t lambda() const { return lambda_; }
		[[nodiscard]] uint32_t size() const { return coeffs_.size(); }
		void swap(ComplexFunction& other) &
		{
			std::swap(sym_, other.sym_);
			std::swap(lambda_, other.lambda_);
			coeffs_.swap(other.coeffs_);
		}
		void negate() & { coeffs_.negate(); }
		[[nodiscard]] bool iszero() const { return coeffs_.iszero(); }
		[[nodiscard]] ComplexFunction clone() const { return ComplexFunction(coeffs_.clone(), lambda_, sym_); }
		// get coefficient of the term (x - x0) ^ {dx} (y - y0) ^ {dy}
		// 0 <= dx, dy and dx + 2 dy <= lambda
		// if symmetry is even or odd, the parity of dx must equals symmetry
		[[nodiscard]] Ring& at(uint32_t dx, uint32_t dy) & { return coeffs_[index(dx, dy)]; }
		[[nodiscard]] const Ring& at(uint32_t dx, uint32_t dy) const& { return coeffs_[index(dx, dy)]; }
		[[nodiscard]] Vector<Ring> flatten() && { return std::move(coeffs_); }

		// z -> 1 - z
		void flip() &
		{
			switch (sym_)
			{
			case FunctionSymmetry::Even: break;
			case FunctionSymmetry::Odd: coeffs_.negate(); break;
			case FunctionSymmetry::Mixed:
			default:
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 1; dx + 2 * dy <= lambda_; dx += 2) at(dx, dy).negate();
				break;
			}
		}
		[[nodiscard]] auto abs() const { return coeffs_.abs(); }
		[[nodiscard]] auto norm() const { return coeffs_.norm(); }

		// project this function to sym-symmetric part
		[[nodiscard]] ComplexFunction proj(FunctionSymmetry sym) &&
		{
			if (sym_ == sym) return std::move(*this);
			ComplexFunction z(lambda_, sym);
			if (sym_ == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym, dx)) z.at(dx, dy).swap(at(dx, dy));
			}
			if (sym == FunctionSymmetry::Mixed)
			{
				for (uint32_t dy = 0; dy <= lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= lambda_; ++dx)
						if (matches(sym_, dx)) z.at(dx, dy).swap(at(dx, dy));
			}
			// otherwise (even to odd or odd to even), proj is vanishing
			return z;
		}
		// project this function to sym-symmetric part
		[[nodiscard]] ComplexFunction proj(FunctionSymmetry sym) const& { return clone().proj(sym); }
		ComplexFunction operator+() const& { return clone(); }
		ComplexFunction operator+() && { return std::move(*this); }
		ComplexFunction operator-() const& { return ComplexFunction(-coeffs_, lambda_, sym_); }
		ComplexFunction operator-() &&
		{
			negate();
			return std::move(*this);
		}
		ComplexFunction& operator+=(const ComplexFunction& f) &
		{
			// TODO(selpo): remove assert
			assert(sym_ == f.sym_);
			coeffs_ += f.coeffs_;
			return *this;
		}
		ComplexFunction& operator-=(const ComplexFunction& f) &
		{
			assert(sym_ == f.sym_);
			coeffs_ -= f.coeffs_;
			return *this;
		}
		template <class T>
		ComplexFunction& operator*=(const T& r) &
		{
			coeffs_ *= r;
			return *this;
		}
		template <class T>
		ComplexFunction& operator/=(const T& r) &
		{
			coeffs_ /= r;
			return *this;
		}
		// TODO(selpo): implement rvalue versions
		friend ComplexFunction operator+(const ComplexFunction& x, const ComplexFunction& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction z(x.lambda_, plus(x.sym_, y.sym_));
			if (x.sym_ == y.sym_)
				z.coeffs_ = x.coeffs_ + y.coeffs_;
			else
				for (uint32_t dy = 0; dy <= x.lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= x.lambda_; ++dx)
					{
						if (matches(x.sym_, dx)) z.at(dx, dy) += x.at(dx, dy);
						if (matches(y.sym_, dx)) z.at(dx, dy) += y.at(dx, dy);
					}
			return z;
		}
		friend ComplexFunction operator-(const ComplexFunction& x, const ComplexFunction& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction z(x.lambda_, plus(x.sym_, y.sym_));
			if (x.sym_ == y.sym_)
				z.coeffs_ = x.coeffs_ - y.coeffs_;
			else
				for (uint32_t dy = 0; dy <= x.lambda_ / 2; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= x.lambda_; ++dx)
					{
						if (matches(x.sym_, dx)) z.at(dx, dy) -= x.at(dx, dy);
						if (matches(y.sym_, dx)) z.at(dx, dy) -= y.at(dx, dy);
					}
			return z;
		}
		friend ComplexFunction mul(const ComplexFunction& x, const ComplexFunction& y)
		{
			assert(x.lambda_ == y.lambda_);
			ComplexFunction z(x.lambda_, times(x.sym_, y.sym_));
			for (uint32_t dy1 = 0; dy1 <= x.lambda_ / 2; ++dy1)
				for (uint32_t dx1 = 0; dx1 + 2 * dy1 <= x.lambda_; ++dx1)
				{
					if (!matches(x.sym_, dx1)) continue;
					for (uint32_t dy2 = 0; dy1 + dy2 <= x.lambda_ / 2; ++dy2)
						for (uint32_t dx2 = 0; dx1 + dx2 + 2 * dy1 + 2 * dy2 <= x.lambda_; ++dx2)
						{
							if (!matches(y.sym_, dx2)) continue;
							z.at(dx1 + dx2, dy1 + dy2) += mul(x.at(dx1, dy1), y.at(dx2, dy2));
						}
				}
			return z;
		}
		template <class R>
		friend ComplexFunction mul_scalar(const R& r, const ComplexFunction& x)
		{
			return ComplexFunction(mul_scalar(r, x.coeffs_), x.lambda_, x.sym_);
		}
		template <class R>
		friend ComplexFunction mul_scalar(const R& r, ComplexFunction&& x)
		{
			return std::move(x *= r);
		}
		template <class R>
		friend ComplexFunction operator/(const ComplexFunction& x, const R& r)
		{
			return ComplexFunction(x.coeffs_ / r, x.lambda_, x.sym_);
		}
		template <class R>
		friend ComplexFunction operator/(ComplexFunction&& x, const R& r)
		{
			return std::move(x /= r);
		}
		friend bool operator==(const ComplexFunction& x, const ComplexFunction& y)
		{
			return x.lambda_ == y.lambda_ && x.sym_ == y.sym_ && x.coeffs_ == y.coeffs_;
		}
		friend bool operator!=(const ComplexFunction& x, const ComplexFunction& y) { return !(x == y); }
		template <class Real>
		[[nodiscard]] ComplexFunction<evaluated_t<Ring>> eval(const Real& x) const
		{
			return ComplexFunction<evaluated_t<Ring>>(coeffs_.eval(x), lambda_, sym_);
		}
	};
	template <class Ring>
	ComplexFunction<polynomialize_t<Ring>> to_pol(Vector<ComplexFunction<Ring>>& coeffs)
	{
		uint32_t lambda = coeffs[0].lambda(), len = coeffs.size();
		auto sym = coeffs[0].symmetry();
		ComplexFunction<polynomialize_t<Ring>> ans(lambda, sym);
		Vector<Ring> v(len);
		for (uint32_t dy = 0; dy <= lambda / 2; ++dy)
			for (uint32_t dx = 0; dx + 2 * dy <= lambda; ++dx)
				if (matches(sym, dx))
				{
					for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs[i].at(dx, dy));
					ans.at(dx, dy) = to_pol(v);
				}
		return ans;
	}
	template <class R>
	std::ostream& operator<<(std::ostream& out, const ComplexFunction<R>& v)
	{
		out << "[";
		auto f = false;
		for (uint32_t dy = 0; dy <= v.lambda() / 2; ++dy)
		{
			// do not show an empty array []
			if (v.symmetry() == FunctionSymmetry::Odd && 2 * dy == v.lambda()) continue;
			if (f) out << ", ";
			auto g = false;
			out << "[";
			for (uint32_t dx = 0; dx + 2 * dy <= v.lambda(); ++dx)
			{
				if (!matches(v.symmetry(), dx)) continue;
				if (g) out << ", ";
				out << v.at(dx, dy);
				g = true;
			}
			out << "]";
			f = true;
		}
		return out << "]";
	}
	// v ^ d = ((1 - z) (1 - z_bar)) ^ d as a function of x, y
	inline ComplexFunction<mpfr::real> v_to_d(const mpfr::real& d, uint32_t lambda)
	{
		// f.at(m, n) = (der x) ^ m (der y) ^ n ((x - 1) ^ 2 - y) ^ d / (n! m!)
		//            = (-1) ^ {n + m} 2 ^ {2 n + m} lf(d, n) lf(2 (d - n), m) / (n! m! 4 ^ d)
		// where lf(x, n) = x (x - 1) ... (x - (n - 1)) (falling factorial)
		ComplexFunction<mpfr::real> f(lambda);
		f.at(0, 0) = mpfr::pow(4, -d);
		for (uint32_t n = 0; n <= lambda / 2; ++n)
		{
			if (n > 0) f.at(0u, n) = f.at(0u, n - 1) * 4 * (-d + (n - 1)) / n;
			for (uint32_t m = 1; m + 2 * n <= lambda; ++m)
				f.at(m, n) = f.at(m - 1, n) * (-4 * d + 2 * (m + 2 * n - 1)) / m;
		}
		return f;
	}
}  // namespace algebra

#endif  // QBOOT_COMPLEX_FUNCTION_HPP_
