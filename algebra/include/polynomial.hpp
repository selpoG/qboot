#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include <cstdint>           // for uint32_t, int32_t
#include <initializer_list>  // for initializer_list
#include <ostream>           // for ostream
#include <type_traits>       // for enable_if, enable_if_t, is_same_v
#include <utility>           // for move

#include "matrix.hpp"   // for Vector, is_intermediate_v, base_ring, ring_dimension_v
#include "real.hpp"     // for real, pow, sgn
#include "real_io.hpp"  // for operator<<

namespace algebra
{
	template <class Ring>
	class Polynomial
	{
		// \sum_{i} coeff_[i] x ^ i
		// the last value of coeff_ must be non-zero
		// zero polynomial is represented by empty coeff_ (coeff_.size() = 0)
		Vector<Ring> coeff_;

	public:
		using base = typename base_ring<Ring>::type;
		using ring = Ring;
		using type = Polynomial<Ring>;
		Polynomial() : coeff_(0) {}
		Polynomial(Polynomial&&) = default;
		Polynomial& operator=(Polynomial&&) = default;
		Polynomial(const Polynomial&) = delete;
		Polynomial& operator=(const Polynomial&) = delete;
		~Polynomial() = default;
		// pow(x, d)
		explicit Polynomial(uint32_t degree) : coeff_(degree + 1) { coeff_[degree] = Ring(base(1)); }
		template <class = std::enable_if<!std::is_same_v<Ring, base>>>
		explicit Polynomial(const base& c) : coeff_(1)
		{
			assert(!c.iszero());
			coeff_[0] = Ring(c);
		}
		template <class T, class = std::enable_if_t<is_intermediate_v<Polynomial, T>>>
		explicit Polynomial(const T& c) : coeff_(1)
		{
			assert(!c.iszero());
			coeff_[0] = Ring(c);
		}
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
			for (uint32_t i = 0; i <= uint32_t(deg); i++) coeff_[i] = coeffs[i].clone();
		}
		Polynomial(std::initializer_list<Ring> coeffs) : coeff_(uint32_t(coeffs.size()))
		{
			uint32_t i = 0;
			for (auto& v : coeffs) coeff_[i++] = v.clone();
			assert(coeff_.size() == 0 || !coeff_[coeff_.size() - 1].iszero());
		}
		[[nodiscard]] bool iszero() const noexcept { return coeff_.size() == 0; }
		[[nodiscard]] int32_t degree() const noexcept { return int32_t(coeff_.size()) - 1; }
		template <class R, class = std::enable_if_t<is_imultipliable_v<Ring, R>>>
		[[nodiscard]] Ring eval(const R& x) const
		{
			if (iszero()) return Ring{};
			Ring s{};
			auto d = uint32_t(degree());
			for (uint32_t i = d; i <= d; --i)
			{
				s *= x;
				s += coeff_[i];
			}
			return s;
		}
		[[nodiscard]] Polynomial deriv() const
		{
			if (degree() <= 0) return Polynomial{};
			Polynomial p(uint32_t(degree() - 1));
			for (uint32_t i = 1; i < coeff_.size(); ++i) p.coeff_[i - 1] = i * coeff_[i];
			return p;
		}
		[[nodiscard]] Polynomial clone() const
		{
			Polynomial p{};
			p.coeff_ = +coeff_;
			return p;
		}
		void swap(Polynomial& other) { coeff_.swap(other.coeff_); }
		void negate() { coeff_.negate(); }
		template <class R, class = std::enable_if_t<is_iaddable_v<Ring, R>>>
		Polynomial& operator+=(const Polynomial<R>& p)
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
			*this = *this + p;
			return *this;
		}
		template <class R, class = std::enable_if_t<!is_polynomial_v<R> && is_iaddable_v<Ring, R>>>
		Polynomial& operator+=(const R& c)
		{
			if (c == 0) return *this;
			if (iszero()) return *this = Polynomial(Ring(c));
			auto t = coeff_[0] + c;
			if (degree() == 0 && t.iszero()) return *this = {};
			coeff_[0] = std::move(t);
			return *this;
		}
		template <class R, class = std::enable_if_t<is_isubtractable_v<Ring, R>>>
		Polynomial& operator-=(const Polynomial<R>& p)
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
		template <class R, class = std::enable_if_t<!is_polynomial_v<R> && is_isubtractable_v<Ring, R>>>
		Polynomial& operator-=(const R& c)
		{
			if (c.iszero()) return *this;
			if (iszero()) return *this = Polynomial(Ring(-c));
			auto t = coeff_[0] - c;
			if (degree() == 0 && t.iszero()) return *this = {};
			coeff_[0] = std::move(t);
			return *this;
		}
		template <class R, class = std::enable_if_t<is_imultipliable_v<Ring, R>>>
		Polynomial& operator*=(const Polynomial<R>& p)
		{
			*this = *this * p;
			return *this;
		}
		template <class R, class = std::enable_if_t<!is_polynomial_v<R> && is_imultipliable_v<Ring, R>>>
		Polynomial& operator*=(const R& c)
		{
			if (iszero()) return *this;
			if constexpr (is_ring_v<R>)
			{
				if (c.iszero()) return *this = {};
			}
			else
			{
				if (c == 0) return *this = {};
			}
			for (uint32_t i = 0; i < coeff_.size(); ++i) coeff_[i] *= c;
			return *this;
		}
		template <class R, class = std::enable_if_t<is_idividable_v<Ring, R>>>
		Polynomial& operator/=(const R& c)
		{
			if (iszero()) return *this;
			for (uint32_t i = 0; i < coeff_.size(); ++i) coeff_[i] /= c;
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
		template <class R, class = std::enable_if_t<is_multipliable_v<Ring, R>>>
		friend union_ring_t<Polynomial, Polynomial<R>> operator*(const Polynomial& p, const Polynomial<R>& q)
		{
			if (p.iszero() || q.iszero()) return union_ring_t<Polynomial, Polynomial<R>>{};
			auto dp = uint32_t(p.degree());
			auto dq = uint32_t(q.degree());
			union_ring_t<Polynomial, Polynomial<R>> r(dp + dq);
			r.coeff_[dp + dq] = {};
			for (uint32_t i = 0; i <= dp; ++i)
				for (uint32_t j = 0; j <= dq; ++j) r.coeff_[i + j] += p[i] * q[j];
			return r;
		}
		template <class R,
		          class = std::enable_if_t<!is_polynomial_v<R> && !is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Polynomial operator*(const Polynomial& p, const R& c)
		{
			if (p.iszero() || c == 0) return Polynomial{};
			auto d = uint32_t(p.degree());
			Polynomial r(d);
			for (uint32_t i = 0; i <= d; ++i) r.coeff_[i] = p[i] * c;
			return r;
		}
		template <class R,
		          class = std::enable_if_t<!is_polynomial_v<R> && !is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Polynomial operator*(const R& c, const Polynomial& p)
		{
			return p * c;
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
		template <class R, class = std::enable_if_t<is_iaddable_v<Ring, R>>>
		friend Polynomial operator+(const Polynomial& p, const R& c)
		{
			if (p.iszero()) return Polynomial(c);
			auto t = p.coeff_[0] + c;
			if (p.degree() == 0 && t.iszero()) return Polynomial{};
			Polynomial q = +p;
			q.coeff_[0] = std::move(t);
			return q;
		}
		template <class R, class = std::enable_if_t<is_iaddable_v<Ring, R>>>
		friend Polynomial operator+(const R& c, const Polynomial& p)
		{
			return p + c;
		}
		template <class R, class = std::enable_if_t<is_isubtractable_v<Ring, R>>>
		friend Polynomial operator-(const Polynomial& p, const R& c)
		{
			if (p.iszero()) return Polynomial(Ring(-c));
			auto t = p.coeff_[0] - c;
			if (p.degree() == 0 && t.iszero()) return Polynomial{};
			Polynomial q = +p;
			q.coeff_[0] = std::move(t);
			return q;
		}
		template <class R, class = std::enable_if_t<is_isubtractable_v<Ring, R>>>
		friend Polynomial operator-(const R& c, const Polynomial& p)
		{
			if (p.iszero()) return Polynomial(Ring(c));
			auto t = p.coeff_[0] - c;
			t.negate();
			if (p.degree() == 0 && t.iszero()) return Polynomial{};
			Polynomial q = -p;
			q.coeff_[0] = std::move(t);
			return q;
		}
		template <class R, class = std::enable_if_t<is_idividable_v<Ring, R>>>
		friend Polynomial operator/(const Polynomial& p, const R& c)
		{
			if (p.iszero()) return Polynomial{};
			auto d = uint32_t(p.degree());
			Polynomial r(d);
			for (uint32_t i = 0; i <= d; ++i) r.coeff_[i] = p[i] / c;
			return r;
		}
		friend bool operator==(const Polynomial& p, const Polynomial& q) { return p.coeff_ == q.coeff_; }
		friend bool operator!=(const Polynomial& p, const Polynomial& q) { return !(p == q); }
		static const Polynomial& one()
		{
			static Polynomial val{0u};
			return val;
		}
		static const Polynomial& var()
		{
			static Polynomial val{1u};
			return val;
		}
		// (a * x + b) ^ d
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		static Polynomial linear_power(const Ring& a, const Ring& b, uint32_t degree)
		{
			if (degree == 0) return Polynomial(0);
			if constexpr (is_mpfr_real_v<Ring>)
			{
				if (a.iszero()) return Polynomial(mpfr::pow(b, degree));
				if (b.iszero()) return Polynomial(mpfr::pow(a, degree), degree);
				Polynomial p(degree);
				// p = \sum_{i = 0}^{degree} binom(d, i) a ^ i b ^ {d - i} x ^ i
				p.coeff_[0] = p.coeff_[degree] = 1;
				for (uint32_t i = 1; 2 * i <= degree; ++i)
					p.coeff_[i] = p.coeff_[degree - i] = base(degree - i + 1) * p[i - 1] / base(i);
				Ring pow(base(1));
				for (uint32_t i = 1; i <= degree; ++i)
				{
					pow *= a;
					p.coeff_[i] *= pow;
				}
				pow = base(1);
				for (uint32_t i = degree - 1; i <= degree; --i)
				{
					pow *= b;
					p.coeff_[i] *= pow;
				}
				return p;
			}
			else
			{
				Polynomial p{0}, q{+b, +a};
				while (degree > 0)
				{
					if (degree % 2 == 1) p = p * q;
					q = q * q;
					degree /= 2;
				}
				return p;
			}
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
		auto d = uint32_t(v.degree());
		for (uint32_t i = 1; i <= d; ++i)
		{
			auto& val = v[i];
			if constexpr (is_mpfr_real_v<Ring>)
			{
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
			}
			else
			{
				if (val.iszero()) continue;
				if (f)
					out << " + (" << val << ") * ";
				else
					out << "(" << val << ") * ";
				out << "x" << ring_dimension_v<Ring>;
			}
			if (i > 1) out << "^" << i;
			f = true;
		}
		return out;
	}
}  // namespace algebra

#endif  // !POLYNOMIAL_HPP_
