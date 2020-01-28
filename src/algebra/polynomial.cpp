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
	void Polynomial::derivate() &
	{
		if (iszero()) return;
		if (degree() == 0)
		{
			*this = {};
			return;
		}
		Vector<real> new_coeff(coeff_.size() - 1);
		for (uint32_t i = 1; i < coeff_.size(); ++i) new_coeff[i - 1] = i * move(coeff_[i]);
		coeff_ = move(new_coeff);
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
	Matrix<real> interpolation_matrix(const Vector<real>& points)
	{
		assert(points.size() > 0);
		auto deg = points.size() - 1;
		Matrix<real> mat(deg + 1, deg + 1);
		for (uint32_t i = 0; i <= deg; ++i)
		{
			mat.at(i, 0) = 1;
			for (uint32_t j = 1; j <= deg; ++j) mat.at(i, j) = points[i] * mat.at(i, j - 1);
		}
		return inverse(mat);
	}
	static int32_t _sign(uint32_t N, const std::vector<uint32_t>& perm);
	int32_t _sign(uint32_t N, const std::vector<uint32_t>& perm)
	{
		int32_t sign = 1;
		std::vector<bool> done(N, false);
		for (uint32_t i = 0; i < N; ++i)
		{
			if (done[i]) continue;
			uint32_t len = 0, j = i;
			while (!done[j])
			{
				done[j] = true;
				j = perm[j];
				++len;
			}
			if (len % 2 == 0) sign *= -1;
		}
		return sign;
	}
	Polynomial determinant(Matrix<Polynomial>&& mat)
	{
		assert(mat.is_square());
		uint32_t N = mat.row();
		if (N == 0) return Polynomial(0u);
		if (N == 1) return move(mat.at(0, 0));
		if (N == 2) return mul(move(mat.at(0, 0)), mat.at(1, 1)) - mul(move(mat.at(0, 1)), mat.at(1, 0));
		std::vector<uint32_t> perm(N);
		for (uint32_t i = 0; i < N; ++i) perm[i] = i;
		Polynomial ans;
		do
		{
			Polynomial prod(real(_sign(N, perm)), 0u);
			for (uint32_t i = 0; i < N; ++i) prod = mul(prod, mat.at(i, perm[i]));
			ans += prod;
		} while (std::next_permutation(perm.begin(), perm.end()));
		return ans;
	}
}  // namespace qboot::algebra
