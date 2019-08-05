#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include "matrix.hpp"

namespace algebra
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Polynomial
	{
		// sum(coeff_[i] * pow(x, i), i)
		// the last value of coeff_ must be non-zero
		// zero polynomial is represented by empty coeff_ (coeff_.size() == 0)
		Vector<Real> coeff_;
		// Real& operator[](size_t p) { return coeff_[p]; }

	public:
		Polynomial() : coeff_(0) {}
		// pow(x, d)
		Polynomial(size_t degree) : coeff_(degree + 1) { coeff_[degree] = 1; }
		// constant polynomial c (c != 0)
		Polynomial(const Real& c) : coeff_(1)
		{
			assert(!c.iszero());
			coeff_[0] = c.clone();
		}
		// c * pow(x, d) (c != 0)
		Polynomial(const Real& c, size_t degree) : coeff_(degree + 1)
		{
			assert(!c.iszero());
			coeff_[degree] = c.clone();
		}
		Polynomial(std::initializer_list<Real> coeffs) : coeff_(coeffs.size())
		{
			size_t i = 0;
			for (auto& v : coeffs) coeff_[i++] = v.clone();
			assert(coeff_.size() == 0 || !coeff_[coeff_.size() - 1].iszero());
		}
		bool iszero() const noexcept { return coeff_.size() == 0; }
		int64_t degree() const noexcept { return int64_t(coeff_.size()) - 1; }
		Real eval(const Real& x) const
		{
			if (iszero()) return Real(0);
			Real s(0);
			for (size_t i = degree(); i <= degree(); i--)
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
		Polynomial& operator+=(const Polynomial& p)
		{
			if (p.iszero()) return *this;
			if (iszero())
			{
				coeff_ = p.coeff_.clone();
				return *this;
			}
			int64_t pd = p.degree(), d = degree();
			if (pd < d || (pd == d && coeff_[d] + p.coeff_[d] != 0))
			{
				for (size_t i = 0; i <= pd; i++) coeff_[i] += p.coeff_[i];
				return *this;
			}
			*this = const_cast<const Polynomial&>(*this) + p;
			return *this;
		}
		Polynomial& operator-=(const Polynomial& p)
		{
			if (p.iszero()) return *this;
			if (iszero())
			{
				coeff_ = p.coeff_.clone();
				coeff_.negate();
				return *this;
			}
			int64_t pd = p.degree(), d = degree();
			if (pd < d || (pd == d && coeff_[d] != p.coeff_[d]))
			{
				for (size_t i = 0; i <= pd; i++) coeff_[i] -= p.coeff_[i];
				return *this;
			}
			*this = const_cast<const Polynomial&>(*this) - p;
			return *this;
		}
		Polynomial& operator*=(const Polynomial& p)
		{
			*this = const_cast<const Polynomial&>(*this) * p;
			return *this;
		}
		// coefficient of pow(x, p)
		const Real& operator[](size_t p) const { return coeff_[p]; }
		friend Polynomial operator+(const Polynomial& p) { return p.clone(); }
		friend Polynomial operator-(const Polynomial& p)
		{
			if (p.iszero()) return Polynomial();
			auto q = p.clone();
			q.coeff_.negate();
			return q;
		}
		friend Polynomial operator*(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero()) return Polynomial();
			if (q.iszero()) return Polynomial();
			Polynomial r(p.degree() + q.degree());
			r.coeff_[p.degree() + q.degree()] = 0;
			for (size_t i = 0; i <= p.degree(); i++)
				for (size_t j = 0; j <= q.degree(); j++) r.coeff_[i + j] += p[i] * q[j];
			return r;
		}
		friend Polynomial operator+(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero()) return +q;
			if (q.iszero()) return +p;
			size_t pd = p.degree(), qd = q.degree(), deg = std::max(pd, qd);
			if (pd == qd)
			{
				while (deg >= 0 && (p[deg] + q[deg]).iszero()) deg--;
				if (deg < 0) return Polynomial();
			}
			Polynomial r(deg);
			for (size_t i = deg; i <= deg; i--) r.coeff_[i] = p[i] + q[i];
			return r;
		}
		friend Polynomial operator-(const Polynomial& p, const Polynomial& q)
		{
			if (p.iszero()) return -q;
			if (q.iszero()) return +p;
			size_t pd = p.degree(), qd = q.degree(), deg = std::max(pd, qd);
			if (pd == qd)
			{
				while (deg >= 0 && p[deg] == q[deg]) deg--;
				if (deg < 0) return Polynomial();
			}
			Polynomial r(deg);
			for (size_t i = deg; i <= deg; i--) r.coeff_[i] = p[i] - q[i];
			return r;
		}
		friend Polynomial operator*(const Polynomial& p, const Real& c)
		{
			if (p.iszero()) return Polynomial();
			assert(!c.iszero());
			Polynomial r(p.degree());
			for (size_t i = 0; i <= p.degree(); i++) r.coeff_[i] = p[i] * c;
			return r;
		}
		friend Polynomial operator*(const Real& c, const Polynomial& p) { return p * c; }
		template <class T = Polynomial>
		friend std::enable_if_t<is_mpfr_real_v<Real>, T> operator/(const Polynomial& p, const Real& c)
		{
			if (p.iszero()) return Polynomial();
			Polynomial r(p.degree());
			for (size_t i = 0; i <= p.degree(); i++) r.coeff_[i] = p[i] / c;
			return r;
		}
		friend bool operator==(const Polynomial& p, const Polynomial& q) { return p.coeff_ == q.coeff_; }
		friend bool operator!=(const Polynomial& p, const Polynomial& q) { return !(p == q); }
		static const Polynomial& one()
		{
			static Polynomial val(0);
			return val;
		}
		static const Polynomial& var()
		{
			static Polynomial val(1);
			return val;
		}
		// pow(a * x + b, d)
		template <class T = Polynomial>
		static std::enable_if_t<is_mpfr_real_v<Real>, T> linear_power(const Real& a, const Real& b, size_t degree)
		{
			if (degree == 0) return Polynomial(0);
			if (a.iszero()) return Polynomial(mpfr::pow(b, static_cast<unsigned long>(degree)));
			if (b.iszero()) return Polynomial(mpfr::pow(a, static_cast<unsigned long>(degree)), degree);
			Polynomial p(degree);
			// p = sum(binom(d, i) * pow(a, i) * pow(b, d - i) * pow(x, i), i)
			p.coeff_[0] = p.coeff_[degree] = 1;
			for (size_t i = 1; 2 * i <= degree; i++)
				p.coeff_[i] = p.coeff_[degree - i] = Real(degree - i + 1) * p[i - 1] / Real(i);
			Real pow(1);
			for (size_t i = 1; i <= degree; i++)
			{
				pow *= a;
				p.coeff_[i] *= pow;
			}
			pow = 1;
			for (size_t i = degree - 1; i <= degree; i--)
			{
				pow *= b;
				p.coeff_[i] *= pow;
			}
			return p;
		}
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
		for (size_t i = 1; i <= v.degree(); i++)
		{
			auto& val = v[i];
			if (val.iszero()) continue;
			if (f)
				out << " + (" << val << ") * ";
			else
				out << "(" << val << ") * ";
			if (i == 1)
				out << "x";
			else
				out << "x^" << i;
			f = true;
		}
		return out;
	}
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	std::ostream& operator<<(std::ostream& out, const Polynomial<mpfr::real<prec, rnd>>& v)
	{
		if (v.iszero()) return out << "0";
		auto f = false;
		{
			auto& val = v[0];
			if (mpfr::sgn(val) != 0)
			{
				out << val;
				f = true;
			}
		}
		for (size_t i = 1; i <= v.degree(); i++)
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
			if (i == 1)
				out << "x";
			else
				out << "x^" << i;
			f = true;
		}
		return out;
	}
}  // namespace algebra

#endif  // !POLYNOMIAL_HPP_
