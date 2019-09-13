#ifndef QBOOT_POLYNOMIAL_HPP_
#define QBOOT_POLYNOMIAL_HPP_

#include <cassert>           // for assert
#include <cstdint>           // for uint32_t, int32_t
#include <initializer_list>  // for initializer_list
#include <ostream>           // for ostream
#include <utility>           // for move

#include "matrix.hpp"  // for Vector, Matrix
#include "real.hpp"    // for real

namespace algebra
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
		Vector<mpfr::real> coeff_;

	public:
		Polynomial() : coeff_(0) {}
		Polynomial(Polynomial&&) noexcept = default;
		Polynomial& operator=(Polynomial&&) & noexcept = default;
		Polynomial(const Polynomial&) = delete;
		Polynomial& operator=(const Polynomial&) = delete;
		~Polynomial() = default;
		// pow(x, d)
		explicit Polynomial(uint32_t degree) : coeff_(degree + 1) { coeff_[degree] = mpfr::real(1); }
		// constant polynomial c
		explicit Polynomial(const mpfr::real& c) : coeff_(1)
		{
			if (c.iszero())
				coeff_ = {};
			else
				coeff_[0] = c.clone();
		}
		// c x ^ d (c != 0)
		Polynomial(const mpfr::real& c, uint32_t degree) : coeff_(degree + 1)
		{
			assert(!c.iszero());
			coeff_[degree] = c.clone();
		}
		explicit Polynomial(Vector<mpfr::real>&& coeffs);
		Polynomial(std::initializer_list<mpfr::real> coeffs) : coeff_(uint32_t(coeffs.size()))
		{
			uint32_t i = 0;
			for (auto& v : coeffs) coeff_[i++] = v.clone();
			assert(coeff_.size() == 0 || !coeff_[coeff_.size() - 1].iszero());
		}
		[[nodiscard]] const mpfr::real* begin() const& noexcept { return coeff_.begin(); }
		[[nodiscard]] const mpfr::real* end() const& noexcept { return coeff_.end(); }
		[[nodiscard]] mpfr::real abs() const { return mpfr::sqrt(norm()); }
		[[nodiscard]] mpfr::real norm() const { return coeff_.norm(); }
		[[nodiscard]] bool iszero() const noexcept { return coeff_.size() == 0; }
		[[nodiscard]] int32_t degree() const noexcept { return int32_t(coeff_.size()) - 1; }
		template <class R>
		[[nodiscard]] mpfr::real eval(const R& x) const
		{
			if (iszero()) return mpfr::real{};
			auto d = uint32_t(degree());
			auto s = coeff_[d];
			for (uint32_t i = d - 1; i <= d; --i)
			{
				s *= x;
				s += coeff_[i];
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
		const mpfr::real& operator[](uint32_t p) const { return coeff_[p]; }
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
	struct evaluated<Polynomial>
	{
		using type = mpfr::real;
	};
	template <class T>
	using polynomialize_t = substitute_t<T, Polynomial>;
	// schematically, to_pol(Vector<Ring>{a, b, c, ...}) = a + b x + c x ^ 2 + ...
	inline auto to_pol(Vector<mpfr::real>* coeffs) { return Polynomial(coeffs->clone()); }
	template <class Ring>
	Matrix<polynomialize_t<Ring>> to_pol(Vector<Matrix<Ring>>* coeffs)
	{
		uint32_t row = coeffs->at(0).row(), column = coeffs->at(0).column(), len = coeffs->size();
		Matrix<polynomialize_t<Ring>> ans(row, column);
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
	Vector<polynomialize_t<Ring>> to_pol(Vector<Vector<Ring>>* coeffs)
	{
		uint32_t sz = coeffs->at(0).size(), len = coeffs->size();
		Vector<polynomialize_t<Ring>> ans(sz);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r < sz; ++r)
		{
			for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs->at(i).at(r));
			ans.at(r) = to_pol(&v);
		}
		return ans;
	}
	// calculate coefficients c of polynomial f(x) s.t. for each i, f(points[i]) = vals[i]
	// vals[i] = c[0] + c[1] points[i] + c[2] points[i] ^ 2 + ... + c[deg] points[i] ^ {deg}
	// evals(polynomial_interpolate(vals, points), points) == vals (up to rounding errors)
	template <class Ring>
	polynomialize_t<Ring> polynomial_interpolate(const Vector<Ring>& vals, const Vector<mpfr::real>& points)
	{
		assert(vals.size() == points.size() && points.size() > 0);
		auto deg = points.size() - 1;
		Matrix<mpfr::real> mat(deg + 1, deg + 1);
		for (uint32_t i = 0; i <= deg; ++i)
		{
			mat.at(i, 0) = 1;
			for (uint32_t j = 1; j <= deg; ++j) mat.at(i, j) = points[i] * mat.at(i, j - 1);
		}
		auto coeffs = dot(inverse(mat), vals);
		return to_pol(&coeffs);
	}
	template <class Ring>
	Vector<evaluated_t<Ring>> evals(const Ring& v, const Vector<mpfr::real>& xs)
	{
		Vector<evaluated_t<Ring>> ans(xs.size());
		for (uint32_t i = 0; i < xs.size(); ++i) ans[i] = v.eval(xs[i]);
		return ans;
	}
}  // namespace algebra

#endif  // QBOOT_POLYNOMIAL_HPP_
