#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cassert>           // for assert
#include <cstdint>           // for uint32_t
#include <initializer_list>  // for initializer_list
#include <memory>            // for unique_ptr, make_unique
#include <ostream>           // for ostream
#include <type_traits>       // for true_type, false_type, is_same_v, enable_if, void_t
#include <utility>           // for move, swap

#include "real.hpp"     // for real, abs, mpfr_prec_t, mpfr_rnd_t
#include "real_io.hpp"  // for operator<<

namespace algebra
{
	template <class, template <class> class, class = std::void_t<>>
	struct detect : std::false_type
	{
	};
	template <class T, template <class> class Check>
	struct detect<T, Check, std::void_t<Check<T>>> : std::true_type
	{
	};
	template <class T>
	using has_iszero_checker = decltype(std::declval<const T&>().iszero());
	template <class T>
	constexpr bool has_iszero = detect<T, has_iszero_checker>::value;
	template <class R>
	bool iszero(const R& v)
	{
		if constexpr (has_iszero<R>)
			return v.iszero();
		else
			return v == 0;
	}
	template <class T>
	struct evaluated;
	template <class T>
	using evaluated_t = typename evaluated<T>::type;
	template <class R, template <class> class Vec>
	struct evaluated<Vec<R>>
	{
		using type = Vec<evaluated_t<R>>;
	};
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct evaluated<mpfr::real<prec, rnd>>
	{
		using type = mpfr::real<prec, rnd>;
	};
	template <class Ring, template <class> class F>
	struct substitute;
	// substitute the most inner template argument R by F<R>,
	// i.e. substitute_t<A<B<...<C<R>>...>>> = A<B<...<C<F<R>>>...>>
	template <class T, template <class> class F>
	using substitute_t = typename substitute<T, F>::type;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd, template <class> class F>
	struct substitute<mpfr::real<prec, rnd>, F>
	{
		using type = F<mpfr::real<prec, rnd>>;
	};
	template <class R, template <class> class Vec, template <class> class F>
	struct substitute<Vec<R>, F>
	{
		using type = Vec<substitute_t<R, F>>;
	};
	template <mpfr_prec_t prec, mpfr_rnd_t rnd, class Real>
	mpfr::real<prec, rnd> eval(const mpfr::real<prec, rnd>& v, [[maybe_unused]] const Real& x)
	{
		return v;
	}
	// interface Swappable<T> {
	//	void swap(T&);
	// };
	// interface Cloneable<T> {
	//	T clone() const;
	// };
	// interface ZeroCheckable<T> {
	//	bool iszero() const;
	// }
	// assert(T{}.iszero());  // where T: DefaultConstructible, ZeroCheckable
	// Ring: Swappable, Clonable, ZeroCheckable, DefaultConstructible
	// Polynomial: Swappable, Clonable, ZeroCheckable, DefaultConstructible
	// Vector<Ring>: Swappable, Clonable
	// Matrix<Ring>: Swappable, Clonable
	// definitions
	template <class Ring>
	class Vector
	{
		std::unique_ptr<Ring[]> arr_;
		uint32_t sz_;

	public:
		Vector() : Vector(0) {}
		explicit Vector(uint32_t len) : arr_(std::make_unique<Ring[]>(len)), sz_(len) {}
		Vector(uint32_t len, const Ring& val) : Vector(len)
		{
			for (uint32_t i = 0; i < len; ++i) arr_[i] = val.clone();
		}
		Vector(Vector&& v) noexcept = default;
		Vector& operator=(Vector&& v) noexcept = default;
		~Vector() = default;
		Vector(const Vector& v) = delete;
		Vector& operator=(const Vector& v) = default;
		Vector(std::initializer_list<Ring> v) : Vector(uint32_t(v.size()))
		{
			uint32_t i = 0;
			for (auto& t : v) arr_[i++] = t.clone();
		}
		[[nodiscard]] Ring& at(uint32_t i) { return arr_[i]; }
		[[nodiscard]] const Ring& at(uint32_t i) const { return arr_[i]; }
		[[nodiscard]] Ring& operator[](uint32_t i) { return at(i); }
		[[nodiscard]] const Ring& operator[](uint32_t i) const { return at(i); }
		[[nodiscard]] const uint32_t& size() const noexcept { return sz_; }
		[[nodiscard]] const Ring* begin() const noexcept { return arr_.get(); }
		[[nodiscard]] const Ring* end() const noexcept { return arr_.get() + sz_; }
		[[nodiscard]] Ring* begin() noexcept { return arr_.get(); }
		[[nodiscard]] Ring* end() noexcept { return arr_.get() + sz_; }
		[[nodiscard]] Vector clone() const
		{
			Vector v(sz_);
			for (uint32_t i = 0; i < sz_; ++i) v.arr_[i] = arr_[i].clone();
			return v;
		}
		void swap(Vector& other)
		{
			arr_.swap(other.arr_);
			std::swap(sz_, other.sz_);
		}
		void negate()
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = -arr_[i];
		}
		[[nodiscard]] bool iszero() const noexcept
		{
			for (uint32_t i = 0; i < sz_; ++i)
				if (!arr_[i].iszero()) return false;
			return true;
		}
		[[nodiscard]] auto abs() const { return norm().sqrt(); }
		[[nodiscard]] auto norm() const
		{
			auto s = arr_[0].norm();
			for (uint32_t i = 1; i < sz_; ++i) s += arr_[i].norm();
			return s;
		}
		Vector& operator+=(const Vector& v)
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] += v.arr_[i];
			return *this;
		}
		Vector& operator-=(const Vector& v)
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] -= v.arr_[i];
			return *this;
		}
		template <class T>
		Vector& operator*=(const T& v)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] *= v;
			return *this;
		}
		template <class T>
		Vector& operator/=(const T& v)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] /= v;
			return *this;
		}
		friend Vector operator+(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		friend Vector operator-(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		template <class R>
		friend Vector mul_scalar(const R& r, const Vector& x)
		{
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = mpfr::mul_scalar(r, x.arr_[i]);
			return z;
		}
		template <class R>
		friend Vector operator/(const Vector& x, const R& r)
		{
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		friend bool operator==(const Vector& x, const Vector& y)
		{
			if (x.sz_ != y.sz_) return false;
			for (uint32_t i = 0; i < x.sz_; ++i)
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
		template <class Ring2 = Ring, class = std::enable_if<std::is_same_v<Ring, Ring2>>>
		friend Ring dot(const Vector& x, const Vector<Ring2>& y)
		{
			assert(x.sz_ == y.sz_);
			if (x.sz_ == 0) return {};
			auto s = mul(x.arr_[0], y.arr_[0]);
			for (uint32_t i = 1; i < x.sz_; ++i) s += mul(x.arr_[i], y.arr_[i]);
			return s;
		}
		template <class Real>
		[[nodiscard]] Vector<evaluated_t<Ring>> eval(const Real& x) const
		{
			Vector<evaluated_t<Ring>> ans(sz_);
			for (uint32_t i = 0; i < sz_; i++) ans[i] = at(i).eval(x);
			return ans;
		}
	};

	template <class Ring>
	class Matrix
	{
		template <class Ring2>
		friend class Matrix;
		Vector<Ring> arr_;
		uint32_t row_, col_;
		Matrix(Vector<Ring>&& v, uint32_t r, uint32_t c) : arr_(std::move(v)), row_(r), col_(c) {}

	public:
		Matrix() : Matrix(0, 0) {}
		Matrix(uint32_t r, uint32_t c) : arr_(r * c), row_(r), col_(c) {}
		Matrix(Matrix&& v) noexcept = default;
		Matrix& operator=(Matrix&& v) noexcept = default;
		~Matrix() = default;
		Matrix(const Matrix& v) = delete;
		Matrix& operator=(const Matrix& v) = default;
		[[nodiscard]] Ring& at(uint32_t r, uint32_t c) { return arr_[r * col_ + c]; }
		[[nodiscard]] const Ring& at(uint32_t r, uint32_t c) const { return arr_[r * col_ + c]; }
		[[nodiscard]] const uint32_t& row() const noexcept { return row_; }
		[[nodiscard]] const uint32_t& column() const noexcept { return col_; }
		[[nodiscard]] bool is_square() const noexcept { return row_ == col_; }
		[[nodiscard]] auto abs() const { return norm().sqrt(); }
		[[nodiscard]] auto norm() const { return arr_.norm(); }
		[[nodiscard]] Matrix clone() const { return Matrix(arr_.clone(), row_, col_); }
		void swap(Matrix& other)
		{
			arr_.swap(other.arr_);
			std::swap(row_, other.row_);
			std::swap(col_, other.col_);
		}
		[[nodiscard]] bool iszero() const noexcept { return arr_.iszero(); }
		void negate() { arr_.negate(); }
		Matrix& operator+=(const Matrix& v)
		{
			assert(v.row_ == row_ && v.col_ == col_);
			arr_ += v.arr_;
			return *this;
		}
		Matrix& operator-=(const Matrix& v)
		{
			assert(v.row_ == row_ && v.col_ == col_);
			arr_ -= v.arr_;
			return *this;
		}
		template <class T>
		Matrix& operator*=(const T& v)
		{
			arr_ *= v;
			return *this;
		}
		template <class T>
		Matrix& operator/=(const T& v)
		{
			arr_ /= v;
			return *this;
		}
		Vector<Ring> flatten() && { return std::move(arr_); }
		void transpose()
		{
			assert(is_square());
			for (uint32_t i = 0; i < row_; ++i)
				for (uint32_t j = 0; j < i; ++j) at(i, j).swap(at(j, i));
		}
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix inverse() const&
		{
			return clone().inverse();
		}
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix inverse() &&
		{
			assert(row_ == col_);
			Matrix inv = Matrix::constant(Ring(1), row_);
			for (uint32_t j = 0; j < row_; ++j)
			{
				uint32_t p = j;
				auto max = Ring{};
				for (uint32_t i = p; i < row_; ++i)
				{
					auto abs = mpfr::abs(at(i, j));
					if (abs > max)
					{
						max = abs;
						p = i;
					}
				}
				inv.swap_row(p, j);
				swap_row(p, j, j);
				auto t = 1 / at(j, j);
				inv.multiply_row(j, t);
				multiply_row(j, t, j);
				for (uint32_t i = j + 1; i < row_; ++i)
				{
					inv.add_row(j, i, -at(i, j));
					add_row(j, i, -at(i, j), j);
				}
			}
			for (uint32_t j = row_ - 1; j < row_; --j)
				for (uint32_t r = j - 1; r < j; --r) inv.add_row(j, r, -at(r, j));
			return inv;
		}
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		void add_row(uint32_t f, uint32_t t, const Ring& x, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) at(t, c) += mul_scalar(x, at(f, c));
		}
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		void multiply_row(uint32_t r, const Ring& x, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) at(r, c) *= x;
		}
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		void swap_row(uint32_t r1, uint32_t r2, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) at(r1, c).swap(at(r2, c));
		}
		static Matrix constant(const Ring& c, uint32_t n)
		{
			Matrix m(n, n);
			for (uint32_t i = 0; i < n; ++i) m.at(i, i) = c.clone();
			return m;
		}
		// calculate the cholesky decomposition L of positive definite matrix, by Cholesky–Banachiewicz algorithm.
		// this = L L^t and L is lower triangular.
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix cholesky_decomposition() const
		{
			assert(is_square());
			Matrix L(row_, row_);
			Ring s;
			for (uint32_t i = 0; i < row_; ++i)
			{
				for (uint32_t j = 0; j <= i; ++j)
				{
					s = {};
					for (uint32_t k = 0; k < j; ++k) s += L.at(i, k) * L.at(j, k);
					s = at(i, j) - s;
					L.at(i, j) = i == j ? s.sqrt() : s / L.at(j, j);
				}
			}
			return L;
		}
		// calculate the inverse matrix of lower triangular matrix
		template <class = std::enable_if<mpfr::is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix lower_triangular_inverse() const
		{
			// this must be lower triangular, i.e., at(i, j) = 0 for i < j
			assert(is_square());
			Matrix res(row_, row_);
			Ring s;
			for (uint32_t i = 0; i < row_; ++i)
			{
				res.at(i, i) = 1 / at(i, i);
				for (uint32_t j = 0; j < i; ++j)
				{
					s = {};
					for (uint32_t k = j; k < i; ++k) s += at(i, k) * res.at(k, j);
					res.at(i, j) = -s / at(i, i);
				}
			}
			return res;
		}
		friend Matrix operator+(const Matrix& x) { return x.clone(); }
		friend Matrix operator-(const Matrix& x) { return Matrix(-x.arr_, x.row_, x.col_); }
		friend Matrix operator+(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			return Matrix(x.arr_ + y.arr_, x.row_, x.col_);
		}
		friend Matrix operator-(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			return Matrix(x.arr_ - y.arr_, x.row_, x.col_);
		}
		template <class R>
		friend Matrix mul_scalar(const R& r, const Matrix& x)
		{
			return Matrix(mul_scalar(r, x.arr_), x.row_, x.col_);
		}
		template <class R>
		friend Matrix operator/(const Matrix& x, const R& r)
		{
			return Matrix(x.arr_ / r, x.row_, x.col_);
		}
		friend bool operator==(const Matrix& x, const Matrix& y)
		{
			return x.row_ == y.row_ && x.col_ == y.col_ && x.arr_ == y.arr_;
		}
		friend bool operator!=(const Matrix& x, const Matrix& y) { return !(x == y); }
		friend Matrix dot(const Matrix& x, const Matrix& y)
		{
			assert(x.col_ == y.row_);
			Matrix z(x.row_, y.col_);
			if (x.col_ > 0)
				for (uint32_t i = 0; i < x.row_; ++i)
					for (uint32_t j = 0; j < y.col_; ++j)
					{
						z.at(i, j) = mul(x.at(i, 0), y.at(0, j));
						for (uint32_t k = 1; k < x.col_; ++k) z.at(i, j) += mul(x.at(i, k), y.at(k, j));
					}
			return z;
		}
		friend Matrix mul(const Matrix& x, const Matrix& y) { return dot(x, y); }
		template <class R>
		friend Vector<R> dot(const Matrix& x, const Vector<R>& y)
		{
			assert(x.col_ == y.size());
			Vector<R> z(x.row_);
			if (x.col_ > 0)
				for (uint32_t i = 0; i < x.row_; ++i)
				{
					z[i] = mul_scalar(x.arr_[i * x.col_], y[0]);
					for (uint32_t j = 1, p = i * x.col_ + 1; j < x.col_; ++j, ++p) z[i] += mul_scalar(x.arr_[p], y[j]);
				}
			return z;
		}
		template <class R>
		friend Vector<R> dot(const Vector<R>& x, const Matrix& y)
		{
			assert(x.size() == y.row_);
			Vector<R> z(y.col_);
			if (y.row_ > 0)
				for (uint32_t i = 0; i < y.col_; ++i)
				{
					z[i] = mul(x[0], y.at(0, i));
					for (uint32_t j = 1; j < y.row_; ++j) z[i] += mul(x[j], y.at(j, i));
				}
			return z;
		}
		template <class Real>
		Matrix<evaluated_t<Ring>> eval(const Real& x) const
		{
			return Matrix<evaluated_t<Ring>>(arr_.eval(x), row_, col_);
		}
	};

	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	std::ostream& operator<<(std::ostream& out, const Vector<Ring>& v)
	{
		out << "[";
		auto f = false;
		for (uint32_t i = 0; i < v.size(); ++i)
		{
			if (f) out << ", ";
			out << v[i];
			f = true;
		}
		return out << "]";
	}

	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	std::ostream& operator<<(std::ostream& out, const Matrix<Ring>& v)
	{
		out << "[";
		auto f = false;
		for (uint32_t i = 0; i < v.row(); ++i)
		{
			if (f) out << ", ";
			auto g = false;
			out << "[";
			for (uint32_t j = 0; j < v.column(); ++j)
			{
				if (g) out << ", ";
				out << v.at(i, j);
				g = true;
			}
			out << "]";
			f = true;
		}
		return out << "]";
	}
}  // namespace algebra

#endif  // MATRIX_HPP_
