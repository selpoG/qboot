#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <array>             // for array
#include <cassert>           // for assert
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <iostream>          // for <<
#include <memory>            // for unique_ptr
#include <ostream>           // for basic_ostream
#include <string>            // for string
#include <vector>            // for vector

#include "real.hpp"  // for iszero, complex, real_prec_t, real_rnd_t

namespace algebra
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Vector
	{
		std::unique_ptr<Real[]> arr_;
		size_t sz_;
		Vector(Real* ptr, size_t len) : arr_(ptr), sz_(len) {}

	public:
		explicit Vector(size_t len) : arr_(std::make_unique<Real[]>(len)), sz_(len) {}
		Vector(size_t len, const Real& val) : arr_(std::make_unique<Real[]>(len)), sz_(len)
		{
			for (size_t i = 0; i < len; i++) arr_[i] = val;
		}
		Vector(Vector&& v) noexcept = default;
		Vector& operator=(Vector&& v) noexcept = default;
		~Vector() = default;
		Vector(const Vector& v) = delete;
		Vector& operator=(const Vector& v) = default;
		[[nodiscard]] Real* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Real* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Real& get(std::size_t i) { return arr_[i]; }
		[[nodiscard]] const Real& get(std::size_t i) const { return arr_[i]; }
		[[nodiscard]] Real& operator[](std::size_t i) { return get(i); }
		[[nodiscard]] const Real& operator[](std::size_t i) const { return get(i); }
		[[nodiscard]] const size_t& size() const noexcept { return sz_; }
		[[nodiscard]] Vector clone() const
		{
			Vector v(sz_);
			for (size_t i = 0; i < sz_; i++) v.arr_[i] = arr_[i];
			return v;
		}
		void negate()
		{
			for (size_t i = 0; i < sz_; i++) arr_[i] = -arr_[i];
		}
		[[nodiscard]] Real abs() const { return mpfr::sqrt(norm()); }
		[[nodiscard]] Real norm() const
		{
			Real s(0);
			for (size_t i = 0; i < sz_; i++) s += arr_[i] * arr_[i];
			return s;
		}
		Vector& operator+=(const Vector& v)
		{
			assert(v.sz_ == sz_);
			for (size_t i = 0; i < sz_; i++) arr_[i] += v.arr_[i];
			return *this;
		}
		Vector& operator-=(const Vector& v)
		{
			assert(v.sz_ == sz_);
			for (size_t i = 0; i < sz_; i++) arr_[i] -= v.arr_[i];
			return *this;
		}
		Vector& operator*=(const Real& r)
		{
			for (size_t i = 0; i < sz_; i++) arr_[i] *= r;
			return *this;
		}
		Vector& operator/=(const Real& r)
		{
			for (size_t i = 0; i < sz_; i++) arr_[i] /= r;
			return *this;
		}
		[[nodiscard]] Vector convolve(const Vector& other) const
		{
			assert(sz_ == other.sz_);
			Vector z(sz_);
			Real s;
			for (size_t i = 0; i < sz_; i++)
			{
				s = 0;
				for (size_t j = 0; j <= i; j++) s += arr_[j] * other.arr_[i - j];
				z.arr_[i] = s;
			}
			return z;
		}
		friend Vector operator+(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Vector z(x.sz_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		friend Vector operator-(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Vector z(x.sz_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		friend Vector operator*(const Vector& x, const Real& r)
		{
			Vector z(x.sz_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] * r;
			return z;
		}
		friend Vector operator*(const Real& r, const Vector& x) { return x * r; }
		friend Vector operator/(const Vector& x, const Real& r)
		{
			Vector z(x.sz_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		friend Real operator*(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Real s(0);
			for (size_t i = 0; i < x.sz_; i++) s += x.arr_[i] * y.arr_[i];
			return s;
		}
		friend bool operator==(const Vector& x, const Vector& y)
		{
			if (x.sz_ != y.sz_) return false;
			for (size_t i = 0; i < x.sz_; i++)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Vector& x, const Vector& y) { return !(x == y); }
		friend Vector operator+(const Vector& x) { return x.clone(); }
		friend Vector operator-(const Vector& x)
		{
			auto y = x.clone();
			y.negate();
			return y;
		}
		template <class Real2>
		friend class Matrix;
		template <class Real2>
		friend class Tensor;
	};

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Matrix
	{
		std::unique_ptr<Real[]> arr_;
		size_t row_, col_, sz_;

	public:
		Matrix(size_t r, size_t c) : arr_(std::make_unique<Real[]>(r * c)), row_(r), col_(c), sz_(r * c) {}
		Matrix(Matrix&& v) noexcept = default;
		Matrix& operator=(Matrix&& v) noexcept = default;
		~Matrix() = default;
		Matrix(const Matrix& v) = delete;
		Matrix& operator=(const Matrix& v) = default;
		[[nodiscard]] Real* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Real* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Real& get(size_t r, size_t c) { return arr_[r * col_ + c]; }
		[[nodiscard]] const Real& get(size_t r, size_t c) const { return arr_[r * col_ + c]; }
		[[nodiscard]] Real& operator[](std::array<std::size_t, 2> i) { return get(i[0], i[1]); }
		[[nodiscard]] const Real& operator[](std::array<std::size_t, 2> i) const { return get(i[0], i[1]); }
		[[nodiscard]] const size_t& size() const noexcept { return sz_; }
		[[nodiscard]] const size_t& row() const noexcept { return row_; }
		[[nodiscard]] const size_t& column() const noexcept { return col_; }
		[[nodiscard]] bool is_square() const noexcept { return row_ == col_; }
		[[nodiscard]] Real abs() const { return mpfr::sqrt(norm()); }
		[[nodiscard]] Real norm() const
		{
			Real s(0);
			for (size_t i = 0; i < sz_; i++) s += arr_[i] * arr_[i];
			return s;
		}
		[[nodiscard]] Matrix clone() const
		{
			Matrix v(row_, col_);
			for (size_t i = 0; i < sz_; i++) v.arr_[i] = arr_[i];
			return v;
		}
		void negate()
		{
			for (size_t i = 0; i < sz_; i++) arr_[i] = -arr_[i];
		}
		Vector<Real> flatten()
		{
			Vector<Real> v(arr_.release(), sz_);
			row_ = col_ = sz_ = 0;
			return v;
		}
		void transpose()
		{
			if (row_ == col_)
			{
				Real t;
				for (size_t i = 0; i < row_; i++)
					for (size_t j = 0; j < i; j++) get(i, j).swap(get(j, i));
			}
			else
			{
				auto t = clone();
				for (size_t i = 0; i < row_; i++)
					for (size_t j = 0; j < col_; j++) get(i, j) = std::move(t.get(j, i));
			}
		}
		friend Matrix operator+(const Matrix& x) { return x.clone(); }
		friend Matrix operator-(const Matrix& x)
		{
			auto z = x.clone();
			z.negate();
			return z;
		}
		friend Matrix operator+(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			Matrix z(x.row_, x.col_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		friend Matrix operator-(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			Matrix z(x.row_, x.col_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		friend Matrix operator*(const Matrix& x, const Real& r)
		{
			Matrix z(x.row_, x.col_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] * r;
			return z;
		}
		friend Matrix operator*(const Real& r, const Matrix& x) { return x * r; }
		friend Matrix operator/(const Matrix& x, const Real& r)
		{
			Matrix z(x.row_, x.col_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		friend Vector<Real> operator*(const Matrix& x, const Vector<Real>& y)
		{
			assert(x.col_ == y.size());
			Vector<Real> z(x.row_);
			Real s;
			for (size_t i = 0; i < x.row_; i++)
			{
				s = 0;
				for (size_t j = 0, p = i * x.col_; j < x.col_; j++, p++) s += x.arr_[p] * y[j];
				z[i] = s;
			}
			return z;
		}
		friend Vector<Real> operator*(const Vector<Real>& x, const Matrix& y)
		{
			assert(x.size() == y.row_);
			Vector<Real> z(y.col_);
			Real s;
			for (size_t i = 0; i < y.col_; i++)
			{
				s = 0;
				for (size_t j = 0; j < y.row_; j++) s += x[j] * y.get(j, i);
				z[i] = s;
			}
			return z;
		}
		friend Matrix operator*(const Matrix& x, const Matrix& y)
		{
			assert(x.col_ == y.row_);
			Matrix<Real> z(x.row_, y.col_);
			Real s;
			for (size_t i = 0; i < x.row_; i++)
			{
				for (size_t j = 0; j < y.col_; j++)
				{
					s = 0;
					for (size_t k = 0; k < x.col_; k++) s += x.get(i, k) * y.get(k, j);
					z.get(i, j) = s;
				}
			}
			return z;
		}
		friend bool operator==(const Matrix& x, const Matrix& y)
		{
			if (x.row_ != y.row_ || x.col_ != y.col_) return false;
			for (size_t i = 0; i < x.sz_; i++)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Matrix& x, const Matrix& y) { return !(x == y); }
		// calculate the cholesky decomposition L of positive definite matrix, by Cholesky–Banachiewicz algorithm.
		// this == L L^t and L is lower triangular.
		[[nodiscard]] Matrix cholesky_decomposition() const
		{
			assert(is_square());
			Matrix L(row_, row_);
			Real s;
			for (size_t i = 0; i < row_; i++)
			{
				for (size_t j = 0; j <= i; j++)
				{
					s = 0;
					for (size_t k = 0; k < j; k++) s += L.get(i, k) * L.get(j, k);
					s = get(i, j) - s;
					L.get(i, j) = i == j ? mpfr::sqrt(s) : s / L.get(j, j);
				}
			}
			return L;
		}
		// calculate the inverse matrix of lower triangular matrix
		[[nodiscard]] Matrix lower_triangular_inverse() const
		{
			// this must be lower triangular, i.e., get(i, j) == 0 for i < j
			assert(is_square());
			Matrix res(row_, row_);
			Real s;
			for (size_t i = 0; i < row_; i++)
			{
				res.get(i, i) = 1 / get(i, i);
				for (size_t j = 0; j < i; j++)
				{
					s = 0;
					for (size_t k = j; k < i; k++) s += get(i, k) * res.get(k, j);
					res.get(i, j) = -s / get(i, i);
				}
			}
			return res;
		}
	};

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class Tensor
	{
		std::unique_ptr<Real[]> arr_;
		size_t row_, col_, len_, sz_;

	public:
		Tensor(size_t r, size_t c, size_t l)
		    : arr_(std::make_unique<Real[]>(r * c * l)), row_(r), col_(c), len_(l), sz_(r * c * l)
		{
		}
		Tensor(Tensor&& v) noexcept = default;
		Tensor& operator=(Tensor&& v) noexcept = default;
		~Tensor() = default;
		Tensor(const Tensor& v) = delete;
		Tensor& operator=(const Tensor& v) = default;
		[[nodiscard]] Real* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Real* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Real& get(size_t r, size_t c, size_t i) { return arr_[(r * col_ + c) * len_ + i]; }
		[[nodiscard]] const Real& get(size_t r, size_t c, size_t i) const { return arr_[(r * col_ + c) * len_ + i]; }
		[[nodiscard]] Real& operator[](std::array<std::size_t, 3> i) { return get(i[0], i[1], i[2]); }
		[[nodiscard]] const Real& operator[](std::array<std::size_t, 3> i) const { return get(i[0], i[1], i[2]); }
		[[nodiscard]] const size_t& size() const noexcept { return sz_; }
		[[nodiscard]] const size_t& row() const noexcept { return row_; }
		[[nodiscard]] const size_t& column() const noexcept { return col_; }
		[[nodiscard]] const size_t& length() const noexcept { return len_; }
		[[nodiscard]] Real abs() const { return mpfr::sqrt(norm()); }
		[[nodiscard]] Real norm() const
		{
			Real s(0);
			for (size_t i = 0; i < sz_; i++) s += arr_[i] * arr_[i];
			return s;
		}
		[[nodiscard]] Tensor clone() const
		{
			Tensor v(row_, col_, len_);
			for (size_t i = 0; i < sz_; i++) v.arr_[i] = arr_[i];
			return v;
		}
		void negate()
		{
			for (size_t i = 0; i < sz_; i++) arr_[i] = -arr_[i];
		}
		Vector<Real> flatten()
		{
			Vector<Real> v(arr_.release(), sz_);
			row_ = col_ = len_ = sz_ = 0;
			return v;
		}
		friend Tensor operator+(const Tensor& x) { return x.clone(); }
		friend Tensor operator-(const Tensor& x)
		{
			auto z = x.clone();
			z.negate();
			return z;
		}
		friend Tensor operator+(const Tensor& x, const Tensor& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_ && x.len_ == y.len_);
			Tensor z(x.row_, x.col_, x.len_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		friend Tensor operator-(const Tensor& x, const Tensor& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_ && x.len_ == y.len_);
			Tensor z(x.row_, x.col_, x.len_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		friend Tensor operator*(const Tensor& x, const Real& r)
		{
			Tensor z(x.row_, x.col_, x.len_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] * r;
			return z;
		}
		friend Tensor operator*(const Real& r, const Tensor& x) { return x * r; }
		friend Tensor operator/(const Tensor& x, const Real& r)
		{
			Tensor z(x.row_, x.col_, x.len_);
			for (size_t i = 0; i < x.sz_; i++) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		friend bool operator==(const Tensor& x, const Tensor& y)
		{
			if (x.row_ != y.row_ || x.col_ != y.col_ || x.len_ != y.len_) return false;
			for (size_t i = 0; i < x.sz_; i++)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Tensor& x, const Tensor& y) { return !(x == y); }
	};

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	std::ostream& operator<<(std::ostream& out, const algebra::Vector<Real>& v)
	{
		out << "[";
		auto f = false;
		for (size_t i = 0; i < v.size(); i++)
		{
			if (f) out << ", ";
			out << v[i];
			f = true;
		}
		return out << "]";
	}

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	std::ostream& operator<<(std::ostream& out, const algebra::Matrix<Real>& v)
	{
		out << "[";
		auto f = false;
		for (size_t i = 0; i < v.row(); i++)
		{
			if (f) out << ", ";
			auto g = false;
			out << "[";
			for (size_t j = 0; j < v.column(); j++)
			{
				if (g) out << ", ";
				out << v.get(i, j);
				g = true;
			}
			out << "]";
			f = true;
		}
		return out << "]";
	}

}  // namespace algebra

#endif  // MATRIX_HPP_
