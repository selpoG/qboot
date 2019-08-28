#ifndef QBOOT_POLYNOMIAL_HPP_
#define QBOOT_POLYNOMIAL_HPP_

#include <cassert>           // for assert
#include <cstdint>           // for uint32_t, int32_t
#include <initializer_list>  // for initializer_list
#include <ostream>           // for ostream
#include <utility>           // for move

#include "matrix.hpp"   // for Vector, Matrix
#include "real.hpp"     // for real, pow, sgn
#include "real_io.hpp"  // for operator<<

namespace algebra
{
	// \sum_{i} coeff_[i] x ^ i
	// the last value of coeff_ must be non-zero
	// zero polynomial is represented by empty coeff_ (coeff_.size() = 0)
	template <class Ring = mpfr::real<1000, MPFR_RNDN>, class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
	class Polynomial
	{
		Vector<Ring> coeff_;

	public:
		Polynomial() : coeff_(0) {}
		Polynomial(Polynomial&&) noexcept = default;
		Polynomial& operator=(Polynomial&&) noexcept = default;
		Polynomial(const Polynomial&) = delete;
		Polynomial& operator=(const Polynomial&) = delete;
		~Polynomial() = default;
		// pow(x, d)
		explicit Polynomial(uint32_t degree) : coeff_(degree + 1) { coeff_[degree] = Ring(1); }
		// constant polynomial c
		explicit Polynomial(const Ring& c) : coeff_(1)
		{
			if (c.iszero())
				coeff_ = {};
			else
				coeff_[0] = c.clone();
		}
		// c x ^ d (c != 0)
		Polynomial(const Ring& c, uint32_t degree) : coeff_(degree + 1)
		{
			assert(!c.iszero());
			coeff_[degree] = c.clone();
		}
		explicit Polynomial(Vector<Ring>&& coeffs) : coeff_(0)
		{
			if (coeffs.size() > 0 && !coeffs[coeffs.size() - 1].iszero())
			{
				coeff_ = std::move(coeffs);
				return;
			}
			int32_t deg = int32_t(coeffs.size()) - 1;
			while (deg >= 0 && coeffs[uint32_t(deg)].iszero()) deg--;
			if (deg < 0) return;
			coeff_ = Vector<Ring>{uint32_t(deg + 1)};
			for (uint32_t i = 0; i <= uint32_t(deg); ++i) coeff_[i].swap(coeffs[i]);
		}
		Polynomial(std::initializer_list<Ring> coeffs) : coeff_(uint32_t(coeffs.size()))
		{
			uint32_t i = 0;
			for (auto& v : coeffs) coeff_[i++] = v.clone();
			assert(coeff_.size() == 0 || !coeff_[coeff_.size() - 1].iszero());
		}
		[[nodiscard]] auto abs() const { return norm().sqrt(); }
		[[nodiscard]] auto norm() const { return coeff_.norm(); }
		[[nodiscard]] bool iszero() const noexcept { return coeff_.size() == 0; }
		[[nodiscard]] int32_t degree() const noexcept { return int32_t(coeff_.size()) - 1; }
		template <class R>
		[[nodiscard]] Ring eval(const R& x) const
		{
			if (iszero()) return Ring{};
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
		void swap(Polynomial& other) { coeff_.swap(other.coeff_); }
		void negate() { coeff_.negate(); }
		Polynomial& operator+=(const Polynomial& p)
		{
			if (p.iszero()) return *this;
			if (iszero())
			{
				coeff_ = p.coeff_.clone();
				return *this;
			}
			auto pd = uint32_t(p.degree());
			auto d = uint32_t(degree());
			if (pd < d || (pd == d && !(coeff_[d] + p.coeff_[d]).iszero()))
			{
				for (uint32_t i = 0; i <= pd; ++i) coeff_[i] += p.coeff_[i];
				return *this;
			}
			return *this = *this + p;
		}
		Polynomial& operator-=(const Polynomial& p)
		{
			if (p.iszero()) return *this;
			if (iszero())
			{
				coeff_ = p.coeff_.clone();
				negate();
				return *this;
			}
			auto pd = uint32_t(p.degree());
			auto d = uint32_t(degree());
			if (pd < d || (pd == d && coeff_[d] != p.coeff_[d]))
			{
				for (uint32_t i = 0; i <= pd; ++i) coeff_[i] -= p.coeff_[i];
				return *this;
			}
			*this = *this - p;
			return *this;
		}
		template <class R>
		Polynomial& operator*=(const R& c)
		{
			if (iszero()) return *this;
			if (algebra::iszero(c)) return *this = {};
			coeff_ *= c;
			return *this;
		}
		template <class R>
		Polynomial& operator/=(const R& c)
		{
			coeff_ /= c;
			return *this;
		}
		// coefficient of x ^ p
		const Ring& operator[](uint32_t p) const { return coeff_[p]; }
		Polynomial operator+() const { return clone(); }
		Polynomial operator-() const
		{
			if (iszero()) return Polynomial();
			auto q = clone();
			q.negate();
			return q;
		}
		friend Polynomial mul(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero() || q.iszero()) return Polynomial{};
			auto dp = uint32_t(p.degree());
			auto dq = uint32_t(q.degree());
			Polynomial r(dp + dq);
			r.coeff_[dp + dq] = {};
			for (uint32_t i = 0; i <= dp; ++i)
				for (uint32_t j = 0; j <= dq; ++j) r.coeff_[i + j] += mul(p[i], q[j]);
			return r;
		}
		template <class R>
		friend Polynomial mul_scalar(const R& c, const Polynomial& p)
		{
			if (p.iszero() || algebra::iszero(c)) return Polynomial{};
			auto d = uint32_t(p.degree());
			Polynomial r(d);
			r.coeff_ *= c;
			return r;
		}
		friend Polynomial operator+(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero()) return +q;
			if (q.iszero()) return +p;
			auto dp = uint32_t(p.degree());
			auto dq = uint32_t(q.degree());
			if (dp > dq)
			{
				Polynomial r = +p;
				for (uint32_t i = 0; i <= dq; ++i) r.coeff_[i] += q.coeff_[i];
				return r;
			}
			if (dp < dq)
			{
				Polynomial r = +q;
				for (uint32_t i = 0; i <= dp; ++i) r.coeff_[i] += p.coeff_[i];
				return r;
			}
			while (dp <= dq && (p.coeff_[dp] + q.coeff_[dp]).iszero()) --dp;
			if (dp > dq) return Polynomial{};
			Polynomial r(dp);
			for (uint32_t i = 0; i <= dp; ++i) r.coeff_[i] = p.coeff_[i] + q.coeff_[i];
			return r;
		}
		friend Polynomial operator-(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero()) return -q;
			if (q.iszero()) return +p;
			auto dp = uint32_t(p.degree());
			auto dq = uint32_t(q.degree());
			if (dp > dq)
			{
				Polynomial r = +p;
				for (uint32_t i = 0; i <= dq; ++i) r.coeff_[i] -= q.coeff_[i];
				return r;
			}
			if (dp < dq)
			{
				Polynomial r = -q;
				for (uint32_t i = 0; i <= dp; ++i) r.coeff_[i] += p.coeff_[i];
				return r;
			}
			while (dp <= dq && (p.coeff_[dp] - q.coeff_[dp]).iszero()) --dp;
			if (dp > dq) return Polynomial{};
			Polynomial r(dp);
			for (uint32_t i = 0; i <= dp; ++i) r.coeff_[i] = p.coeff_[i] - q.coeff_[i];
			return r;
		}
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
	};

	template <class Ring>
	std::ostream& operator<<(std::ostream& out, const Polynomial<Ring>& v)
	{
		if (v.iszero()) return out << "0";
		auto f = false;
		{
			auto& val = v[0];
			if (!val.iszero())
			{
				out << val;
				f = true;
			}
		}
		auto d = uint32_t(v.degree());
		for (uint32_t i = 1; i <= d; ++i)
		{
			auto& val = v[i];
			auto sgn = mpfr::sgn(val);
			if (sgn == 0) continue;
			if (f)
			{
				if (val == 1)
					out << " + ";
				else if (val == -1)
					out << " - ";
				else if (sgn > 0)
					out << " + " << val << " * ";
				else
					out << " - " << -val << " * ";
			}
			else if (val == -1)
				out << "-";
			else if (val != -1)
				out << val << " * ";
			out << "x";
			if (i > 1) out << "^" << i;
			f = true;
		}
		return out;
	}
	template <class Ring>
	struct evaluated<Polynomial<Ring>>
	{
		using type = Ring;
	};
	template <class T>
	using Polynomial_ = Polynomial<T>;
	template <class T>
	using polynomialize_t = substitute_t<T, Polynomial_>;
	// do not allow nested polynomial
	template <class R, template <class> class F>
	struct substitute<Polynomial<R>, F>
	{
	};
	// schematically, to_pol(Vector<Ring>{a, b, c, ...}) = a + b x + c x ^ 2 + ...
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	auto to_pol(Vector<mpfr::real<prec, rnd>>& coeffs)
	{
		return Polynomial<mpfr::real<prec, rnd>>(coeffs.clone());
	}
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	Matrix<polynomialize_t<Ring>> to_pol(Vector<Matrix<Ring>>& coeffs)
	{
		uint32_t row = coeffs[0].row(), column = coeffs[0].column(), len = coeffs.size();
		Matrix<polynomialize_t<Ring>> ans(row, column);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r < row; ++r)
			for (uint32_t c = 0; c < column; ++c)
			{
				for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs[i].at(r, c));
				ans.at(r, c) = to_pol(v);
			}
		return ans;
	}
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	Vector<polynomialize_t<Ring>> to_pol(Vector<Vector<Ring>>& coeffs)
	{
		uint32_t sz = coeffs[0].size(), len = coeffs.size();
		Vector<polynomialize_t<Ring>> ans(sz);
		Vector<Ring> v(len);
		for (uint32_t r = 0; r < sz; ++r)
		{
			for (uint32_t i = 0; i < len; ++i) v[i].swap(coeffs[i].at(r));
			ans.at(r) = to_pol(v);
		}
		return ans;
	}
	// calculate coefficients c of polynomial f(x) s.t. for each i, f(points[i]) = vals[i]
	// vals[i] = c[0] + c[1] points[i] + c[2] points[i] ^ 2 + ... + c[deg] points[i] ^ {deg}
	// evals(polynomial_interpolate(vals, points), points) == vals (up to rounding errors)
	template <class Ring, class Real = mpfr::real<1000, MPFR_RNDN>>
	polynomialize_t<Ring> polynomial_interpolate(const Vector<Ring>& vals, const Vector<Real>& points)
	{
		assert(vals.size() == points.size() && points.size() > 0);
		auto deg = points.size() - 1;
		Matrix<Real> mat(deg + 1, deg + 1);
		for (uint32_t i = 0; i <= deg; ++i)
		{
			mat.at(i, 0) = 1;
			for (uint32_t j = 1; j <= deg; ++j) mat.at(i, j) = points[i] * mat.at(i, j - 1);
		}
		auto coeffs = dot(mat.inverse(), vals);
		return to_pol(coeffs);
	}
	template <class Ring, class Real>
	Vector<evaluated_t<Ring>> evals(const Ring& v, const Vector<Real>& xs)
	{
		Vector<evaluated_t<Ring>> ans(xs.size());
		for (uint32_t i = 0; i < xs.size(); ++i) ans[i] = v.eval(xs[i]);
		return ans;
	}
}  // namespace algebra

#endif  // QBOOT_POLYNOMIAL_HPP_
