#ifndef QBOOT_ALGEBRA_MATRIX_HPP_
#define QBOOT_ALGEBRA_MATRIX_HPP_

#include <cassert>           // for assert
#include <concepts>          // for is_constructible_v, is_nothrow_destructible_v, same_as
#include <cstdint>           // for uint32_t
#include <initializer_list>  // for initializer_list
#include <memory>            // for unique_ptr, make_unique
#include <ostream>           // for ostream
#include <type_traits>       // for true_type, false_type, is_same_v, enable_if, void_t
#include <utility>           // for move, swap

#include "qboot/mp/real.hpp"  // for real

namespace qboot::algebra
{
	template <class, template <class> class, class = std::void_t<>>
	struct _detect : std::false_type
	{
	};
	template <class T, template <class> class Check>
	struct _detect<T, Check, std::void_t<Check<T>>> : std::true_type
	{
	};
	template <class T>
	using _has_iszero_checker = decltype(std::declval<const T&>().iszero());
	template <class T>
	inline constexpr bool _has_iszero = _detect<T, _has_iszero_checker>::value;
	template <class R>
	bool iszero(const R& v)
	{
		if constexpr (_has_iszero<R>)
			return v.iszero();
		else
			return v == 0;
	}
	template <class T>
	struct _evaluated;
	template <class T>
	using _evaluated_t = typename _evaluated<T>::type;
	template <class R, template <class> class Vec>
	struct _evaluated<Vec<R>>
	{
		using type = Vec<_evaluated_t<R>>;
	};
	template <>
	struct _evaluated<mp::real>
	{
		using type = mp::real;
	};
	template <class R, class S>
	struct _substitute;
	// substitute the most inner template argument R by S,
	// i.e. _substitute_t<A<B<...<C<R>>...>>> = A<B<...<C<S>>...>>
	template <class T, class S>
	using _substitute_t = typename _substitute<T, S>::type;
	template <class S>
	struct _substitute<mp::real, S>
	{
		using type = S;
	};
	template <class R, template <class> class Vec, class S>
	struct _substitute<Vec<R>, S>
	{
		using type = Vec<_substitute_t<R, S>>;
	};
	inline mp::real eval(const mp::real& v, [[maybe_unused]] const mp::real& x) { return v; }

	// clang does not implement std::default_initializable
	template <class T>
	concept _destructible = std::is_nothrow_destructible_v<T>;
	template <class T, class... Args>
	concept _constructible_from = _destructible<T>&& std::is_constructible_v<T, Args...>;
	template <class T>
	concept _default_initializable = _constructible_from<T>&& requires
	{
		T{};
	}
	&&requires { ::new (static_cast<void*>(nullptr)) T; };

	template <class T>
	concept _ring_base = requires(T x, T y, const T c)
	{
		{
			x.swap(y)
		}
		->std::same_as<void>;
		{
			c.clone()
		}
		->std::same_as<T>;
		// assert(T{}.iszero());
		{
			c.iszero()
		}
		->std::same_as<bool>;
	}
	&&_default_initializable<T>;
	template <class T>
	concept _ring_ops = requires(const T& x, const T& y, T& r)
	{
		r += y;
		r -= y;
		r.negate();
		{
			x + y
		}
		->std::same_as<T>;
		{
			+x
		}
		->std::same_as<T>;
		{
			-x
		}
		->std::same_as<T>;
		{
			x == y
		}
		->std::same_as<bool>;
		x.norm();
	};
	template <class T>
	concept Ring = _ring_base<T>&& _ring_ops<T>;
	template <class T>
	concept Algebra = requires(const T x, const T y)
	{
		{
			mul(x, y)
		}
		->std::same_as<T>;
	};
	template <Ring R>
	class Vector
	{
		std::unique_ptr<R[]> arr_;
		uint32_t sz_;

	public:
		void _reset() &&
		{
			arr_.reset();
			sz_ = 0;
		}
		Vector() : arr_{}, sz_(0) {}
		explicit Vector(uint32_t len) : arr_(std::make_unique<R[]>(len)), sz_(len) {}
		Vector(uint32_t len, const R& val) : Vector(len)
		{
			for (uint32_t i = 0; i < len; ++i) arr_[i] = val.clone();
		}
		Vector(Vector&& v) noexcept = default;
		Vector& operator=(Vector&& v) & noexcept = default;
		~Vector() = default;
		Vector(const Vector& v) = delete;
		Vector& operator=(const Vector& v) = delete;
		Vector(std::initializer_list<R> v) : Vector(uint32_t(v.size()))
		{
			uint32_t i = 0;
			for (auto& t : v) arr_[i++] = t.clone();
		}
		[[nodiscard]] R& at(uint32_t i) & { return arr_[i]; }
		[[nodiscard]] const R& at(uint32_t i) const& { return arr_[i]; }
		[[nodiscard]] R& operator[](uint32_t i) & { return at(i); }
		[[nodiscard]] const R& operator[](uint32_t i) const& { return at(i); }
		[[nodiscard]] const uint32_t& size() const noexcept { return sz_; }
		[[nodiscard]] const R* begin() const& noexcept { return arr_.get(); }
		[[nodiscard]] const R* end() const& noexcept { return arr_.get() + sz_; }
		[[nodiscard]] R* begin() & noexcept { return arr_.get(); }
		[[nodiscard]] R* end() & noexcept { return arr_.get() + sz_; }
		[[nodiscard]] Vector clone() const
		{
			Vector v(sz_);
			for (uint32_t i = 0; i < sz_; ++i) v.arr_[i] = arr_[i].clone();
			return v;
		}
		void swap(Vector& other) &
		{
			arr_.swap(other.arr_);
			std::swap(sz_, other.sz_);
		}
		void negate() &
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = -arr_[i];
		}
		[[nodiscard]] bool iszero() const noexcept
		{
			for (uint32_t i = 0; i < sz_; ++i)
				if (!arr_[i].iszero()) return false;
			return true;
		}
		[[nodiscard]] auto abs() const { return mp::sqrt(norm()); }
		[[nodiscard]] auto norm() const
		{
			auto s = arr_[0].norm();
			for (uint32_t i = 1; i < sz_; ++i) s += arr_[i].norm();
			return s;
		}
		Vector& operator+=(const Vector& v) &
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] += v.arr_[i];
			return *this;
		}
		Vector& operator-=(const Vector& v) &
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] -= v.arr_[i];
			return *this;
		}
		template <class S>
		Vector& operator*=(const S& v) &
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] *= v;
			return *this;
		}
		template <class S>
		Vector& operator/=(const S& v) &
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
		friend Vector operator+(Vector&& x, const Vector& y) { return std::move(x += y); }
		friend Vector operator+(const Vector& x, Vector&& y) { return std::move(y += x); }
		friend Vector operator+(Vector&& x, Vector&& y)
		{
			x += y;
			std::move(y)._reset();
			return std::move(x);
		}
		friend Vector operator-(const Vector& x, const Vector& y)
		{
			assert(x.sz_ == y.sz_);
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		friend Vector operator-(Vector&& x, const Vector& y) { return std::move(x -= y); }
		friend Vector operator-(const Vector& x, Vector&& y)
		{
			assert(x.sz_ == y.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) y.arr_[i] = x.arr_[i] - y.arr_[i];
			return std::move(y);
		}
		friend Vector operator-(Vector&& x, Vector&& y)
		{
			x -= y;
			std::move(y)._reset();
			return std::move(x);
		}
		template <class S>
		friend Vector mul_scalar(const S& r, const Vector& x)
		{
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = mul_scalar(r, x.arr_[i]);
			return z;
		}
		template <class S>
		friend Vector mul_scalar(const S& r, Vector&& x)
		{
			return std::move(x *= r);
		}
		template <class S>
		friend Vector operator/(const Vector& x, const S& r)
		{
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		template <class S>
		friend Vector operator/(Vector&& x, const S& r)
		{
			return std::move(x /= r);
		}
		friend bool operator==(const Vector& x, const Vector& y)
		{
			if (x.sz_ != y.sz_) return false;
			for (uint32_t i = 0; i < x.sz_; ++i)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Vector& x, const Vector& y) { return !(x == y); }
		Vector operator+() const& { return clone(); }
		Vector operator+() && { return std::move(*this); }
		Vector operator-() const&
		{
			auto y = clone();
			y.negate();
			return y;
		}
		Vector operator-() &&
		{
			negate();
			return std::move(*this);
		}
		friend R dot(const Vector& x, const Vector& y) requires Algebra<R>
		{
			assert(x.sz_ == y.sz_);
			if (x.sz_ == 0) return {};
			R s = mul(x.arr_[0], y.arr_[0]);
			for (uint32_t i = 1; i < x.sz_; ++i) s += mul(x.arr_[i], y.arr_[i]);
			return s;
		}
		[[nodiscard]] Vector<_evaluated_t<R>> eval(const mp::real& x) const
		{
			Vector<_evaluated_t<R>> ans(sz_);
			for (uint32_t i = 0; i < sz_; ++i) ans[i] = at(i).eval(x);
			return ans;
		}
	};

	template <Ring R>
	class Matrix
	{
		template <Ring R2>
		friend class Matrix;
		Vector<R> arr_;
		uint32_t row_, col_;
		Matrix(Vector<R>&& v, uint32_t r, uint32_t c) : arr_(std::move(v)), row_(r), col_(c) {}

	public:
		void _reset() &&
		{
			std::move(arr_)._reset();
			row_ = col_ = 0;
		}
		Matrix() : Matrix(0, 0) {}
		Matrix(uint32_t r, uint32_t c) : arr_(r * c), row_(r), col_(c) {}
		Matrix(Matrix&& v) noexcept = default;
		Matrix& operator=(Matrix&& v) & noexcept = default;
		~Matrix() = default;
		Matrix(const Matrix& v) = delete;
		Matrix& operator=(const Matrix& v) = delete;
		[[nodiscard]] R& at(uint32_t r, uint32_t c) & { return arr_[r * col_ + c]; }
		[[nodiscard]] const R& at(uint32_t r, uint32_t c) const& { return arr_[r * col_ + c]; }
		[[nodiscard]] const uint32_t& row() const noexcept { return row_; }
		[[nodiscard]] const uint32_t& column() const noexcept { return col_; }
		[[nodiscard]] bool is_square() const noexcept { return row_ == col_; }
		[[nodiscard]] auto abs() const { return mp::sqrt(norm()); }
		[[nodiscard]] auto norm() const { return arr_.norm(); }
		[[nodiscard]] Matrix clone() const { return Matrix(arr_.clone(), row_, col_); }
		void swap(Matrix& other) &
		{
			arr_.swap(other.arr_);
			std::swap(row_, other.row_);
			std::swap(col_, other.col_);
		}
		[[nodiscard]] bool iszero() const noexcept { return arr_.iszero(); }
		void negate() & { arr_.negate(); }
		Matrix& operator+=(const Matrix& v) &
		{
			assert(v.row_ == row_ && v.col_ == col_);
			arr_ += v.arr_;
			return *this;
		}
		Matrix& operator-=(const Matrix& v) &
		{
			assert(v.row_ == row_ && v.col_ == col_);
			arr_ -= v.arr_;
			return *this;
		}
		template <class S>
		Matrix& operator*=(const S& v) &
		{
			arr_ *= v;
			return *this;
		}
		template <class S>
		Matrix& operator/=(const S& v) &
		{
			arr_ /= v;
			return *this;
		}
		Vector<R> flatten() && { return std::move(arr_); }
		void transpose() &
		{
			assert(is_square());
			for (uint32_t i = 0; i < row_; ++i)
				for (uint32_t j = 0; j < i; ++j) at(i, j).swap(at(j, i));
		}
		static Matrix constant(const R& c, uint32_t n)
		{
			Matrix m(n, n);
			for (uint32_t i = 0; i < n; ++i) m.at(i, i) = c.clone();
			return m;
		}
		Matrix operator+() const& { return clone(); }
		Matrix operator+() && { return std::move(*this); }
		Matrix operator-() const& { return Matrix(-arr_, row_, col_); }
		Matrix operator-() &&
		{
			negate();
			return std::move(*this);
		}
		// v^t M v
		template <class = std::enable_if<std::is_same_v<R, mp::real>>>
		[[nodiscard]] R inner_product(const Vector<R>& v) const  // requires std::same_as<R, mp::real>
		{
			assert(is_square() && row_ == v.size());
			mp::real s{};
			for (uint32_t r = 0; r < row_; ++r)
				for (uint32_t c = 0; c < row_; ++c) mp::fma(s, mul(v[r], v[c]), at(r, c), s);
			return s;
		}
		friend Matrix operator+(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			return Matrix(x.arr_ + y.arr_, x.row_, x.col_);
		}
		friend Matrix operator+(Matrix&& x, const Matrix& y) { return std::move(x += y); }
		friend Matrix operator+(const Matrix& x, Matrix&& y) { return std::move(y += x); }
		friend Matrix operator+(Matrix&& x, Matrix&& y)
		{
			x += y;
			std::move(y)._reset();
			return std::move(x);
		}
		friend Matrix operator-(const Matrix& x, const Matrix& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			return Matrix(x.arr_ - y.arr_, x.row_, x.col_);
		}
		friend Matrix operator-(Matrix&& x, const Matrix& y) { return std::move(x -= y); }
		friend Matrix operator-(const Matrix& x, Matrix&& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			y.arr_ = x.arr_ - std::move(y.arr_);
			return std::move(y);
		}
		friend Matrix operator-(Matrix&& x, Matrix&& y)
		{
			x -= y;
			std::move(y)._reset();
			return std::move(x);
		}
		template <class S>
		friend Matrix mul_scalar(const S& r, const Matrix& x)
		{
			return Matrix(mul_scalar(r, x.arr_), x.row_, x.col_);
		}
		template <class S>
		friend Matrix mul_scalar(const S& r, Matrix&& x)
		{
			return std::move(x *= r);
		}
		template <class S>
		friend Matrix operator/(const Matrix& x, const S& r)
		{
			return Matrix(x.arr_ / r, x.row_, x.col_);
		}
		template <class S>
		friend Matrix operator/(Matrix&& x, const S& r)
		{
			return std::move(x /= r);
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
		template <class S>
		friend Vector<S> dot(const Matrix& x, const Vector<S>& y)
		{
			assert(x.col_ == y.size());
			Vector<S> z(x.row_);
			if (x.col_ > 0)
				for (uint32_t i = 0; i < x.row_; ++i)
				{
					z[i] = mul_scalar(x.arr_[i * x.col_], y[0]);
					for (uint32_t j = 1, p = i * x.col_ + 1; j < x.col_; ++j, ++p) z[i] += mul_scalar(x.arr_[p], y[j]);
				}
			return z;
		}
		template <class S>
		friend Vector<S> dot(const Vector<S>& x, const Matrix& y)
		{
			assert(x.size() == y.row_);
			Vector<S> z(y.col_);
			if (y.row_ > 0)
				for (uint32_t i = 0; i < y.col_; ++i)
				{
					z[i] = mul(x[0], y.at(0, i));
					for (uint32_t j = 1; j < y.row_; ++j) z[i] += mul(x[j], y.at(j, i));
				}
			return z;
		}
		[[nodiscard]] Matrix<_evaluated_t<R>> eval(const mp::real& x) const
		{
			return Matrix<_evaluated_t<R>>(arr_.eval(x), row_, col_);
		}
	};

	mp::real determinant(Matrix<mp::real>&& mat);
	inline mp::real determinant(const Matrix<mp::real>& mat) { return determinant(mat.clone()); }
	Matrix<mp::real> inverse(Matrix<mp::real>&& mat);
	inline Matrix<mp::real> inverse(const Matrix<mp::real>& mat) { return inverse(mat.clone()); }
	// calculate the cholesky decomposition L of positive definite matrix, by Cholesky–Banachiewicz algorithm.
	// this = L L^t and L is lower triangular.
	[[nodiscard]] Matrix<mp::real> cholesky_decomposition(const Matrix<mp::real>& mat);
	// calculate the inverse matrix of lower triangular matrix
	[[nodiscard]] Matrix<mp::real> lower_triangular_inverse(const Matrix<mp::real>& mat);

	template <Ring R>
	std::ostream& operator<<(std::ostream& out, const Vector<R>& v)
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

	template <Ring R>
	std::ostream& operator<<(std::ostream& out, const Matrix<R>& v)
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
}  // namespace qboot::algebra

#endif  // QBOOT_ALGEBRA_MATRIX_HPP_
