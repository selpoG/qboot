#include "matrix.hpp"

#include <algorithm>  // for max, min, sort
#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <iomanip>    // for operator<<, setprecision
#include <iostream>   // for cout
#include <set>        // for set
#include <sstream>    // for basic_stringbuf<>::int_type, basic_stringbuf<>:...

#include "complex_io.hpp"

using std::map, std::set, std::vector, std::string, std::ostringstream, std::cout, std::endl, std::initializer_list;

namespace algebra
{
	const Real epsilon = mpfr::pow2<_prec, _rnd>(1, -_chop);
	void SparseRow::add(const SparseRow& v, const Real& r)
	{
		if (mpfr::iszero(r) != 0) return;
		for (auto&& it : v)
		{
			auto it2 = elements_.find(it.first);
			if (it2 != elements_.end())
			{
				auto x = it2->second + it.second * r;
				if (iszero(x))
					elements_.erase(it2);
				else
					it2->second = x;
			}
			else
				elements_[it.first] = it.second * r;
		}
	}
	void SparseRow::add(const SparseRow& v, const Complex& r)
	{
		if (iszero(r)) return;
		for (auto&& it : v)
		{
			auto it2 = elements_.find(it.first);
			if (it2 != elements_.end())
			{
				auto x = it2->second + it.second * r;
				if (iszero(x))
					elements_.erase(it2);
				else
					it2->second = x;
			}
			else
				elements_[it.first] = it.second * r;
		}
	}

	SparseRow& SparseRow::operator*=(const Real& v)
	{
		if (mpfr::iszero(v) != 0)
			elements_.clear();
		else
			for (auto&& it : elements_) it.second *= v;
		return *this;
	}
	SparseRow& SparseRow::operator*=(const Complex& v)
	{
		if (iszero(v))
			elements_.clear();
		else
			for (auto&& it : elements_) it.second *= v;
		return *this;
	}

	SparseRow& SparseRow::operator+=(const SparseRow& v)
	{
		add(v, Real(1));
		return *this;
	}

	SparseRow& SparseRow::operator-=(const SparseRow& v)
	{
		add(v, Real(-1));
		return *this;
	}

	const Complex& SparseRow::operator[](size_t c) const { return elements_.find(c)->second; }

	string SparseRow::str() const
	{
		ostringstream oss;
		// oss << std::fixed;
		// oss << std::setprecision(6);
		auto f = false;
		for (auto& e : elements_)
		{
			if (f) oss << ", ";
			oss << e.first << ": " << e.second;
			f = true;
		}
		return oss.str();
	}
	string SparseRow::str(size_t C) const
	{
		ostringstream oss;
		// oss << std::fixed;
		// oss << std::setprecision(6);
		for (size_t c = 0; c < C; c++)
		{
			if (c > 0) oss << " ";
			oss << get(c);
		}
		return oss.str();
	}

	Real chop(const Real& x)
	{
		if (mpfr::abs(x) < epsilon) return Real();
		return x;
	}

	Complex chop(const Complex& z) { return Complex(chop(z.real()), chop(z.imag())); }

	SparseRow conj(const SparseRow& v)
	{
		auto w = v;
		return w.conjugate();
	}

	constexpr double threshold_diagonal = 1e-12;
	constexpr size_t debug_interval = 1000;
	Matrix Matrix::null_space()
	{
		for (auto& r : a)
		{
			auto m = Real(-1);
			for (auto& elm : r) m = std::max(m, abs(elm.second));
			if (m > 0) r *= 1 / m;
		}
		cout << "density: " << double(nonzero_) / double(row_ * column_) << endl;
		map<size_t, size_t> cnt;
		for (auto& r : a) cnt[r.size()]++;
		for (auto& t : cnt) cout << t.first << " -> " << t.second << endl;
		set<size_t> zero;
		for (auto& r : a)
			if (r.size() == 1)
				for (auto& e : r) zero.insert(e.first);
		cout << "zero: " << zero.size() << endl;
		for (size_t r = 0; r < row_; r++)
		{
			if (a[r].size() < 2) continue;
			vector<size_t> rem;
			for (auto& it : a[r])
				if (zero.count(it.first) != 0) rem.push_back(it.first);
			for (auto& c : rem) set_at(r, c, Complex(0));
		}
		std::sort(a.begin(), a.end(), [](const SparseRow& _a, const SparseRow& _b) { return _a.size() < _b.size(); });
		size_t pos = 0, ma = 0;
		vector<size_t> iperm(row_), diag, sing;
		cout << "density: " << double(nonzero_) / double(row_ * column_) << endl;
		for (size_t j = 0; j < column_; j++)
		{
			size_t p = pos;
			// Complex max = 0;
			while (p < row_ && abs(at(p, j)) < threshold_diagonal) p++;
			// for (size_t i = p; i < row_; i++)
			// 	if (abs(at(i, j)) > max)
			// 	{
			// 		max = abs(at(i, j));
			// 		p = i;
			// 	}
			if (p >= row_)
			{
				sing.push_back(j);
				continue;
			}
			diag.push_back(j);
			// cout << "j = " << j << ", pos = " << pos << ", p = " << p << endl;
			if (p != pos)
			{
				swap_row(p, pos);
				auto t = iperm[p];
				iperm[p] = iperm[pos];
				iperm[pos] = t;
			}
			multiply_row(pos, 1 / at(pos, j));
			for (size_t i = pos + 1; i < row_; i++) add_row(pos, i, -at(i, j));
			pos++;
			ma = std::max(ma, nonzero_);
			if (pos % debug_interval == 0)
				cout << "(" << pos << ") density: " << double(nonzero_) / double(row_ * column_) << endl;
		}
		cout << "max density: " << double(ma) / double(row_ * column_) << endl;
		cout << "density: " << double(nonzero_) / double(row_ * column_) << endl;
		if (pos == 0) return identity(column_);
		while (--pos > 0)
			for (size_t r = 0; r < pos; r++) add_row(pos, r, chop(-at(r, diag[pos])));
		Matrix ans(sing.size(), column_);
		for (size_t i = 0; i < sing.size(); i++) ans.set_at(i, sing[i], Complex(1));
		for (size_t i = 0; i < sing.size(); i++)
			for (size_t r = 0; r < diag.size(); r++) ans.set_at(i, diag[r], chop(-at(r, sing[i])));
		cout << "diag = " << diag.size() << endl;
		for (auto x : sing) cout << x << endl;
		return ans;
	}

	Matrix Matrix::orthogonalize() const
	{
		Matrix m(row_, column_);
		for (size_t r = 0; r < row_; r++)
		{
			m.a[r] = a[r];
			for (size_t k = 0; k < r; k++) m.a[r].add(m.a[k], -inner_product(m.a[k], m.a[r]));
			m.a[r] *= 1 / sqrt(inner_product(m.a[r], m.a[r]).real());
		}
		return m;
	}

	size_t Matrix::rank() const
	{
		Matrix m = *this;
		size_t pos = 0;
		for (size_t j = 0; j < column_; j++)
		{
			size_t p = pos;
			Real max{};
			for (size_t i = p; i < row_; i++)
				if (abs(m.at(i, j)) > max)
				{
					max = abs(m.at(i, j));
					p = i;
				}
			if (max < threshold_diagonal) continue;
			m.swap_row(p, pos);
			m.multiply_row(pos, 1 / m.at(pos, j));
			for (size_t i = pos + 1; i < row_; i++) m.add_row(pos, i, -m.at(i, j));
			pos++;
		}
		return pos;
	}

	void Matrix::row_reduce()
	{
		size_t pos = 0;
		vector<size_t> diag;
		for (size_t j = 0; j < column_; j++)
		{
			size_t p = pos;
			Real max{};
			for (size_t i = p; i < row_; i++)
				if (abs(at(i, j)) > max)
				{
					max = abs(at(i, j));
					p = i;
				}
			if (max < threshold_diagonal) continue;
			diag.push_back(j);
			swap_row(p, pos);
			multiply_row(pos, 1 / at(pos, j));
			for (size_t i = pos + 1; i < row_; i++) add_row(pos, i, -at(i, j));
			pos++;
		}
		if (pos == 0)
		{
			a.clear();
			nonzero_ = 0;
			row_ = 0;
			return;
		}
		while (--pos > 0)
			for (size_t r = pos - 1; r < pos; r--) add_row(pos, r, -at(r, diag[pos]));
	}

	Matrix Matrix::constant(const Complex& v, size_t len)
	{
		Matrix m(len, len);
		for (size_t r = 0; r < len; r++) m.a[r].set(r, v);
		return m;
	}
	Matrix Matrix::identity(size_t size) { return constant(Complex(1), size); }

	SparseRow& Matrix::operator[](size_t r) { return a[r]; }

	const SparseRow& Matrix::operator[](size_t r) const { return a[r]; }

	const Complex& Matrix::operator[](initializer_list<size_t> rc) const
	{
		assert(rc.size() == 2);
		auto it = rc.begin();
		auto r = *it;
		auto c = *(++it);
		return a[r][c];
	}

	Matrix Matrix::inverse() const
	{
		Matrix m = *this, inv = constant(Complex(1), row_);
		for (size_t j = 0; j < column_; j++)
		{
			size_t p = j;
			auto max = Real(0);
			for (size_t i = p; i < row_; i++)
				if (abs(m.at(i, j)) > max)
				{
					max = abs(m.at(i, j));
					p = i;
				}
			inv.swap_row(p, j);
			m.swap_row(p, j);
			inv.multiply_row(j, 1 / m.at(j, j));
			m.multiply_row(j, 1 / m.at(j, j));
			for (size_t i = j + 1; i < row_; i++)
			{
				inv.add_row(j, i, -m.at(i, j));
				m.add_row(j, i, -m.at(i, j));
			}
		}
		for (size_t j = column_ - 1; j < column_; j--)
			for (size_t r = j - 1; r < j; r--)
			{
				inv.add_row(j, r, -m.at(r, j));
				m.add_row(j, r, -m.at(r, j));
			}
		return inv;
	}

	string Matrix::str(size_t R1, size_t C) const
	{
		ostringstream oss;
		for (size_t r = 0; r < std::min(R1, row_); r++)
		{
			if (r > 0) oss << "\n";
			oss << a[r].str(std::min(C, column_));
		}
		return oss.str();
	}
	constexpr std::streamsize default_prec = 300;
	string Matrix::str2() const
	{
		ostringstream oss;
		oss << std::setprecision(default_prec);
		size_t c = 0;
		for (size_t r = 0; r < row_; r++)
			for (auto& it : a[r])
				if (abs(it.second) > threshold_diagonal)
				{
					oss << "(" << r << ", " << it.first << "): " << it.second << "\n";
					c++;
				}
		oss << c << " elements";
		return oss.str();
	}

	Matrix Matrix::kronecker_product(const Matrix& x, const Matrix& y)
	{
		Matrix ans(x.row_ * y.row_, x.column_ * y.column_);
		for (size_t rx = 0; rx < x.row_; rx++)
			for (size_t ry = 0; ry < y.row_; ry++)
				for (auto& tx : x[rx])
					for (auto& ty : y[ry])
						ans.set_at(rx * y.row_ + ry, tx.first * y.column_ + ty.first, tx.second * ty.second);
		return ans;
	}
	Vector Vector::kronecker_product(const Vector& x, const Vector& y)
	{
		Vector ans(x.dim_ * y.dim_);
		for (auto& tx : x)
			for (auto& ty : y) ans.set(tx.first * y.dim_ + ty.first, tx.second * ty.second);
		return ans;
	}
}  // namespace algebra
