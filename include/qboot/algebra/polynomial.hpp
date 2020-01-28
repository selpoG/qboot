#ifndef QBOOT_ALGEBRA_POLYNOMIAL_HPP_
#define QBOOT_ALGEBRA_POLYNOMIAL_HPP_

#include <cassert>           // for assert
#include <cstdint>           // for uint32_t, int32_t
#include <initializer_list>  // for initializer_list
#include <ostream>           // for ostream
#include <utility>           // for move

#include "qboot/algebra/matrix.hpp"  // for Vector, Matrix
#include "qboot/mp/real.hpp"         // for real

namespace qboot::algebra
{
	class Polynomial;
	Polynomial mul(const Polynomial& p, const Polynomial& q);
	Polynomial operator+(const Polynomial& p, const Polynomial& q);
	Polynomial operator-(const Polynomial& p, const Polynomial& q);
	std::ostream& operator<<(std::ostream& out, const Polynomial& v);
	// \sum_{i} coeff_[i] x ^ i
	// the last value of coeff_ must be non-zero
	// zero polynomial is represented by empty coeff_ (coeff_.size() = 0)
	class Polynomial
	{
		Vector<mp::real> coeff_;

	public:
		void _reset() && { std::move(coeff_)._reset(); }
		Polynomial() : coeff_(0) {}
		Polynomial(Polynomial&&) noexcept = default;
		Polynomial& operator=(Polynomial&&) & noexcept = default;
		Polynomial(const Polynomial&) = delete;
		Polynomial& operator=(const Polynomial&) = delete;
		~Polynomial() = default;
		// pow(x, d)
		explicit Polynomial(uint32_t degree) : coeff_(degree + 1) { coeff_[degree] = mp::real(1); }
		// constant polynomial c
		explicit Polynomial(const mp::real& c) : coeff_(1)
		{
			if (c.iszero())
				coeff_ = {};
			else
				coeff_[0] = c.clone();
		}
		// c x ^ d (c != 0)
		Polynomial(const mp::real& c, uint32_t degree) : coeff_(degree + 1)
		{
			assert(!c.iszero());
			coeff_[degree] = c.clone();
		}
		explicit Polynomial(Vector<mp::real>&& coeffs);
		Polynomial(std::initializer_list<mp::real> coeffs) : coeff_(uint32_t(coeffs.size()))
		{
			uint32_t i = 0;
			for (auto& v : coeffs) coeff_[i++] = v.clone();
			assert(coeff_.size() == 0 || !coeff_[coeff_.size() - 1].iszero());
		}
		// coefficient of x ^ p (p <= deg)
		[[nodiscard]] const mp::real& at(uint32_t p) const { return coeff_.at(p); }
		[[nodiscard]] const mp::real* begin() const& noexcept { return coeff_.begin(); }
		[[nodiscard]] const mp::real* end() const& noexcept { return coeff_.end(); }
		[[nodiscard]] mp::real abs() const { return mp::sqrt(norm()); }
		[[nodiscard]] mp::real norm() const { return coeff_.norm(); }
		[[nodiscard]] bool iszero() const noexcept { return coeff_.size() == 0; }
		[[nodiscard]] int32_t degree() const noexcept { return int32_t(coeff_.size()) - 1; }
		template <class R>
		[[nodiscard]] mp::real eval(const R& x) const
		{
			if (iszero()) return mp::real{};
			auto d = uint32_t(degree());
			auto s = coeff_[d];
			for (uint32_t i = d - 1; i <= d; --i)
			{
				if constexpr (std::is_same_v<R, mp::real>)
					mp::fma(s, s, x, coeff_[i]);
				else
				{
					s *= x;
					s += coeff_[i];
				}
			}
			return s;
		}
		[[nodiscard]] Polynomial clone() const
		{
			Polynomial p{};
			p.coeff_ = +coeff_;
			return p;
		}
		void swap(Polynomial& other) & { coeff_.swap(other.coeff_); }
		void negate() & { coeff_.negate(); }
		void derivate() &;
		Polynomial derivative()
		{
			auto t = clone();
			t.derivate();
			return t;
		}
		Polynomial& operator+=(const Polynomial& p) &;
		Polynomial& operator-=(const Polynomial& p) &;
		template <class R>
		Polynomial& operator*=(const R& c) &
		{
			if (iszero()) return *this;
			if (algebra::iszero(c)) return *this = {};
			coeff_ *= c;
			return *this;
		}
		template <class R>
		Polynomial& operator/=(const R& c) &
		{
			coeff_ /= c;
			return *this;
		}
		// coefficient of x ^ p
		const mp::real& operator[](uint32_t p) const { return coeff_[p]; }
		Polynomial operator+() const& { return clone(); }
		Polynomial operator+() && { return std::move(*this); }
		Polynomial operator-() const&
		{
			if (iszero()) return Polynomial();
			auto q = clone();
			q.negate();
			return q;
		}
		Polynomial operator-() &&
		{
			negate();
			return std::move(*this);
		}
		// this *= a + x
		void _mul_linear(const mp::real& a) &;
		friend Polynomial mul(const Polynomial& p, const Polynomial& q);
		template <class R>
		friend Polynomial mul_scalar(const R& c, const Polynomial& p)
		{
			if (p.iszero() || algebra::iszero(c)) return Polynomial{};
			auto d = uint32_t(p.degree());
			Polynomial r(d);
			r.coeff_ *= c;
			return r;
		}
		friend Polynomial operator+(const Polynomial& p, const Polynomial& q);
		friend Polynomial operator-(const Polynomial& p, const Polynomial& q);
		template <class R>
		friend Polynomial operator/(const Polynomial& p, const R& c)
		{
			if (p.iszero()) return Polynomial{};
			auto d = uint32_t(p.degree());
			Polynomial r(d);
			r.coeff_ = p.coeff_ / c;
			return r;
		}
		friend bool operator==(const Polynomial& p, const Polynomial& q) { return p.coeff_ == q.coeff_; }
		friend bool operator!=(const Polynomial& p, const Polynomial& q) { return !(p == q); }
		friend std::ostream& operator<<(std::ostream& out, const Polynomial& v);
	};

	template <>
	struct _evaluated<Polynomial>
	{
		using type = mp::real;
	};
	template <class T>
	using _polynomialize_t = _substitute_t<T, Polynomial>;
	// schematically, to_pol(Vector<Ring>{a, b, c, ...}) = a + b x + c x ^ 2 + ...
	inline auto to_pol(Vector<mp::real>* coeffs) { return Polynomial(coeffs->clone()); }
	template <class Ring>
	Matrix<_polynomialize_t<Ring>> to_pol(Vector<Matrix<Ring>>* coeffs)
	{
		uint32_t row = coeffs->at(0).row(), column = coeffs->at(0).column(), len = coeffs->size();
		Matrix<_polynomialize_t<Ring>> ans(row, column);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r < row; ++r)
			for (uint32_t c = 0; c < column; ++c)
			{
				for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs->at(i).at(r, c));
				ans.at(r, c) = to_pol(&v);
			}
		return ans;
	}
	template <class Ring>
	Vector<_polynomialize_t<Ring>> to_pol(Vector<Vector<Ring>>* coeffs)
	{
		uint32_t sz = coeffs->at(0).size(), len = coeffs->size();
		Vector<_polynomialize_t<Ring>> ans(sz);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r < sz; ++r)
		{
			for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs->at(i).at(r));
			ans.at(r) = to_pol(&v);
		}
		return ans;
	}
	Matrix<mp::real> interpolation_matrix(const Vector<mp::real>& points);
	template <class Ring>
	_polynomialize_t<Ring> polynomial_interpolate(const Vector<Ring>& vals,
	                                              const Matrix<mp::real>& interpolation_matrix)
	{
		assert(vals.size() == interpolation_matrix.row() && vals.size() > 0 && interpolation_matrix.is_square());
		auto coeffs = dot(interpolation_matrix, vals);
		return to_pol(&coeffs);
	}
	// calculate coefficients c of polynomial f(x) s.t. for each i, f(points[i]) = vals[i]
	// vals[i] = c[0] + c[1] points[i] + c[2] points[i] ^ 2 + ... + c[deg] points[i] ^ {deg}
	// evals(polynomial_interpolate(vals, points), points) == vals (up to rounding errors)
	template <class Ring>
	_polynomialize_t<Ring> polynomial_interpolate(const Vector<Ring>& vals, const Vector<mp::real>& points)
	{
		assert(vals.size() == points.size() && points.size() > 0);
		return polynomial_interpolate(vals, interpolation_matrix(points));
	}
	template <class Ring>
	Vector<_evaluated_t<Ring>> evals(const Ring& v, const Vector<mp::real>& xs)
	{
		Vector<_evaluated_t<Ring>> ans(xs.size());
		for (uint32_t i = 0; i < xs.size(); ++i) ans[i] = v.eval(xs[i]);
		return ans;
	}
	Polynomial determinant(Matrix<Polynomial>&& mat);
	inline Polynomial determinant(const Matrix<Polynomial>& mat) { return determinant(mat.clone()); }
}  // namespace qboot::algebra

#endif  // QBOOT_ALGEBRA_POLYNOMIAL_HPP_
