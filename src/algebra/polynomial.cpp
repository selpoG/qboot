#include "qboot/algebra/polynomial.hpp"

#include "qboot/mp/real.hpp"  // for real, pow, sgn

using qboot::mp::real;
using std::move, std::ostream;

namespace qboot::algebra
{
	Polynomial::Polynomial(Vector<real>&& coeffs) : coeff_(0)
	{
		if (coeffs.size() > 0 && !coeffs[coeffs.size() - 1].iszero())
		{
			coeff_ = move(coeffs);
			return;
		}
		int32_t deg = int32_t(coeffs.size()) - 1;
		while (deg >= 0 && coeffs[uint32_t(deg)].iszero()) deg--;
		if (deg < 0)
		{
			move(coeffs)._reset();
			return;
		}
		coeff_ = Vector<real>{uint32_t(deg + 1)};
		for (uint32_t i = 0; i <= uint32_t(deg); ++i) coeff_[i].swap(coeffs[i]);
		move(coeffs)._reset();
	}
	Polynomial& Polynomial::operator+=(const Polynomial& p) &
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
	Polynomial& Polynomial::operator-=(const Polynomial& p) &
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
	void Polynomial::_mul_linear(const real& a) &
	{
		Vector<real> v(coeff_.size() + 1);
		v[0] = coeff_[0] * a;
		for (uint32_t i = 1; i < coeff_.size(); ++i) mp::fma(v[i], coeff_[i], a, coeff_[i - 1]);
		v[coeff_.size()] = coeff_[coeff_.size() - 1];
		coeff_ = move(v);
	}
	Polynomial mul(const Polynomial& p, const Polynomial& q)
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
	Polynomial operator+(const Polynomial& p, const Polynomial& q)
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
	Polynomial operator-(const Polynomial& p, const Polynomial& q)
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
	ostream& operator<<(ostream& out, const Polynomial& v)
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
			auto sgn = mp::sgn(val);
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
			else if (val != 1)
				out << val << " * ";
			out << "x";
			if (i > 1) out << "^" << i;
			f = true;
		}
		return out;
	}
}  // namespace qboot::algebra
