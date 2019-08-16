#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <array>             // for array
#include <cassert>           // for assert
#include <cstdint>           // for uint32_t
#include <initializer_list>  // for initializer_list
#include <memory>            // for unique_ptr, make_unique
#include <ostream>           // for ostream
#include <type_traits>       // for true_type, false_type, is_same_v, enable_if_t, enable_if, integral_constant
#include <utility>           // for move, swap

#include "real.hpp"     // for real, abs, mpfr_prec_t, mpfr_rnd_t
#include "real_io.hpp"  // for operator<<

namespace algebra
{
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class Polynomial;
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class Vector;
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class Matrix;
	template <class Ring = mpfr::real<1000, MPFR_RNDN>>
	class Tensor;
	// meta functions
	template <class, class, template <class, class> class, class = std::void_t<>>
	struct detect2 : std::false_type
	{
	};
	template <class T1, class T2, template <class, class> class Check>
	struct detect2<T1, T2, Check, std::void_t<Check<T1, T2>>> : std::true_type
	{
	};
	template <class T1, class T2>
	using is_addable_checker = decltype(std::declval<const T1&>() + std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_addable_v = detect2<T1, T2, is_addable_checker>::value;
	template <class T1, class T2>
	using is_subtractable_checker = decltype(std::declval<const T1&>() - std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_subtractable_v = detect2<T1, T2, is_subtractable_checker>::value;
	template <class T1, class T2>
	using is_multipliable_checker = decltype(std::declval<const T1&>() * std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_multipliable_v = detect2<T1, T2, is_multipliable_checker>::value;
	template <class T1, class T2>
	using is_dividable_checker = decltype(std::declval<const T1&>() / std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_dividable_v = detect2<T1, T2, is_dividable_checker>::value;
	template <class T1, class T2>
	using is_iaddable_checker = decltype(std::declval<T1&>() += std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_iaddable_v = detect2<T1, T2, is_iaddable_checker>::value;
	template <class T1, class T2>
	using is_isubtractable_checker = decltype(std::declval<T1&>() -= std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_isubtractable_v = detect2<T1, T2, is_isubtractable_checker>::value;
	template <class T1, class T2>
	using is_imultipliable_checker = decltype(std::declval<T1&>() *= std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_imultipliable_v = detect2<T1, T2, is_imultipliable_checker>::value;
	template <class T1, class T2>
	using is_idividable_checker = decltype(std::declval<T1&>() /= std::declval<const T2&>());
	template <class T1, class T2>
	constexpr bool is_idividable_v = detect2<T1, T2, is_idividable_checker>::value;
	// T is a ring or not (ring is mpfr::real or a polynomial ring of some ring)
	template <class T>
	struct is_ring;
	template <class T>
	inline constexpr bool is_ring_v = is_ring<T>::value;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct is_ring<mpfr::real<prec, rnd>> : std::true_type
	{
	};
	template <class R>
	struct is_ring<Polynomial<R>> : is_ring<R>
	{
	};
	template <class R>
	struct is_ring : std::false_type
	{
	};
	template <class T>
	struct is_linear_space;
	template <class T>
	inline constexpr bool is_linear_space_v = is_linear_space<T>::value;
	template <class R>
	struct is_linear_space<Vector<R>> : std::true_type
	{
	};
	template <class R>
	struct is_linear_space<Matrix<R>> : std::true_type
	{
	};
	template <class R>
	struct is_linear_space<Tensor<R>> : std::true_type
	{
	};
	template <class R>
	struct is_linear_space : std::false_type
	{
	};
	template <class T>
	struct is_polynomial;
	template <class T>
	inline constexpr bool is_polynomial_v = is_polynomial<T>::value;
	template <class R>
	struct is_polynomial<Polynomial<R>> : std::true_type
	{
	};
	template <class R>
	struct is_polynomial : std::false_type
	{
	};
	// T is mpfr::real or not
	template <class T>
	struct is_mpfr_real;
	template <class T>
	inline constexpr bool is_mpfr_real_v = is_mpfr_real<T>::value;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct is_mpfr_real<mpfr::real<prec, rnd>> : std::true_type
	{
	};
	template <class T>
	struct is_mpfr_real : std::false_type
	{
	};
	// calculate base ring (= mpfr::real) of T
	template <class T>
	struct base_ring;
	template <class T>
	using base_ring_t = typename base_ring<T>::type;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct base_ring<mpfr::real<prec, rnd>>
	{
		using type = mpfr::real<prec, rnd>;
	};
	template <class Ring>
	struct base_ring<Polynomial<Ring>>
	{
		using type = base_ring_t<Ring>;
	};
	template <class Ring>
	struct base_ring<Vector<Ring>>
	{
		using type = base_ring_t<Ring>;
	};
	template <class Ring>
	struct base_ring<Matrix<Ring>>
	{
		using type = base_ring_t<Ring>;
	};
	template <class Ring>
	struct base_ring<Tensor<Ring>>
	{
		using type = base_ring_t<Ring>;
	};
	// R is a non-trivial (R != E and R != mpfr::real) intermediate ring of E / mpfr::real or not
	template <class E, class R, class = void>
	struct is_intermediate;
	template <class E, class R>
	constexpr bool is_intermediate_v = is_intermediate<E, R>::value;
	template <class E>
	struct is_intermediate<Polynomial<E>, E, std::enable_if_t<!is_mpfr_real_v<E>>> : std::true_type
	{
	};
	template <class E, class R>
	struct is_intermediate<
	    Polynomial<E>, R,
	    std::enable_if_t<!std::is_same_v<Polynomial<E>, R> && !std::is_same_v<E, R> && !is_mpfr_real_v<R>>>
	    : is_intermediate<E, R>
	{
	};
	template <class E, class R, class>
	struct is_intermediate : std::false_type
	{
	};
	template <class R, uint32_t n>
	struct add_variables;
	template <class R, uint32_t n>
	using add_variables_t = typename add_variables<R, n>::type;
	template <class R>
	struct add_variables<R, 0>
	{
		using type = R;
	};
	template <class R, uint32_t n>
	struct add_variables
	{
		using type = Polynomial<add_variables_t<R, n - 1>>;
	};
	template <class R>
	struct ring_dimension;
	template <class R>
	constexpr uint32_t ring_dimension_v = ring_dimension<R>::value;
	template <mpfr_prec_t prec, mpfr_rnd_t rnd>
	struct ring_dimension<mpfr::real<prec, rnd>> : std::integral_constant<uint32_t, 0>
	{
	};
	template <class R>
	struct ring_dimension<Polynomial<R>> : std::integral_constant<uint32_t, ring_dimension_v<R> + 1>
	{
	};
	template <class R>
	struct ring_dimension : std::integral_constant<uint32_t, 0>
	{
	};
	template <class R1, class R2, class = void>
	struct union_ring;
	template <class R1, class R2>
	using union_ring_t = typename union_ring<R1, R2>::type;
	template <class R1, class R2>
	struct union_ring<
	    R1, R2,
	    std::enable_if_t<is_ring_v<R1> && is_ring_v<R2> && std::is_same_v<base_ring_t<R1>, base_ring_t<R2>> &&
	                     ring_dimension_v<R1> >= ring_dimension_v<R2>>>
	{
		using type = R1;
	};
	template <class R1, class R2>
	    struct union_ring < R1,
	    R2,
	    std::enable_if_t<is_ring_v<R1> && is_ring_v<R2> && std::is_same_v<base_ring_t<R1>, base_ring_t<R2>> &&
	                     ring_dimension_v<R1><ring_dimension_v<R2>>>
	{
		using type = R2;
	};
	template <class R1, class R2>
	struct union_ring<R1, R2, std::enable_if_t<is_ring_v<R1> && !is_ring_v<R2>>>
	{
		using type = R1;
	};
	template <class R1, class R2>
	struct union_ring<R1, R2, std::enable_if_t<!is_ring_v<R1> && is_ring_v<R2>>>
	{
		using type = R2;
	};
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
	// Tensor<Ring>: Swappable, Clonable
	// definitions
	template <class Ring>
	class Vector
	{
		std::unique_ptr<Ring[]> arr_;
		uint32_t sz_;
		Vector(Ring* ptr, uint32_t len) : arr_(ptr), sz_(len) {}

	public:
		template <class Ring2>
		friend class Vector;
		template <class Ring2>
		friend class Matrix;
		template <class Ring2>
		friend class Tensor;
		using base = typename base_ring<Ring>::type;
		using ring = Ring;
		using type = Vector<Ring>;
		explicit Vector(uint32_t len) : arr_(std::make_unique<Ring[]>(len)), sz_(len) {}
		Vector(uint32_t len, const Ring& val) : arr_(std::make_unique<Ring[]>(len)), sz_(len)
		{
			for (uint32_t i = 0; i < len; ++i) arr_[i] = val.clone();
		}
		Vector(Vector&& v) noexcept = default;
		Vector& operator=(Vector&& v) noexcept = default;
		~Vector() = default;
		Vector(const Vector& v) = delete;
		Vector& operator=(const Vector& v) = default;
		template <class T, class = std::enable_if_t<is_intermediate_v<Ring, T> ||
		                                            (std::is_same_v<T, base> && !std::is_same_v<T, Ring>)>>
		explicit Vector(const Vector<T>& v) : Vector(v.sz_)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = Ring(v.arr_[i]);
		}
		Vector(std::initializer_list<Ring> v) : Vector(uint32_t(v.size()))
		{
			uint32_t i = 0;
			for (auto& t : v) arr_[i++] = t.clone();
		}
		[[nodiscard]] Ring* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Ring* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Ring& get(uint32_t i) { return arr_[i]; }
		[[nodiscard]] const Ring& get(uint32_t i) const { return arr_[i]; }
		[[nodiscard]] Ring& operator[](uint32_t i) { return get(i); }
		[[nodiscard]] const Ring& operator[](uint32_t i) const { return get(i); }
		[[nodiscard]] const uint32_t& size() const noexcept { return sz_; }
		[[nodiscard]] const Ring* cbegin() const noexcept { return arr_.get(); }
		[[nodiscard]] const Ring* cend() const noexcept { return arr_.get() + sz_; }
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
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const
		{
			Ring s{};
			for (uint32_t i = 0; i < sz_; ++i) s += arr_[i] * arr_[i];
			return s;
		}
		template <class T, class = std::enable_if_t<is_iaddable_v<Ring, T>>>
		Vector& operator+=(const Vector<T>& v)
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] += v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_isubtractable_v<Ring, T>>>
		Vector& operator-=(const Vector<T>& v)
		{
			assert(v.sz_ == sz_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] -= v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_imultipliable_v<Ring, T>>>
		Vector& operator*=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] *= r;
			return *this;
		}
		template <class T, class = std::enable_if_t<is_idividable_v<Ring, T>>>
		Vector& operator/=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] /= r;
			return *this;
		}
		[[nodiscard]] Vector convolve(const Vector& other) const
		{
			assert(sz_ == other.sz_);
			Vector z(sz_);
			Ring s;
			for (uint32_t i = 0; i < sz_; ++i)
			{
				s = {};
				for (uint32_t j = 0; j <= i; ++j) s += arr_[j] * other.arr_[i - j];
				z.arr_[i] = std::move(s);
			}
			return z;
		}
		template <class R, class = std::enable_if_t<is_addable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator+(const Vector& x, const Vector<R>& y)
		{
			assert(x.sz_ == y.sz_);
			Vector<union_ring_t<Ring, R>> z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<is_subtractable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator-(const Vector& x, const Vector<R>& y)
		{
			assert(x.sz_ == y.sz_);
			Vector<union_ring_t<Ring, R>> z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator*(const Vector& x, const R& r)
		{
			Vector<union_ring_t<Ring, R>> z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z[i] = x.arr_[i] * r;
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator*(const R& r, const Vector& x)
		{
			return x * r;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_dividable_v<Ring, R>>>
		friend Vector operator/(const Vector& x, const R& r)
		{
			Vector z(x.sz_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		template <class R, class = std::enable_if_t<is_multipliable_v<Ring, R>>>
		friend union_ring_t<Ring, R> operator*(const Vector& x, const Vector<R>& y)
		{
			assert(x.sz_ == y.sz_);
			union_ring_t<Ring, R> s{};
			for (uint32_t i = 0; i < x.sz_; ++i) s += x.arr_[i] * y.arr_[i];
			return s;
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
	};

	template <class Ring>
	class Matrix
	{
		std::unique_ptr<Ring[]> arr_;
		uint32_t row_, col_, sz_;

	public:
		using base = typename base_ring<Ring>::type;
		using ring = Ring;
		using type = Matrix<Ring>;
		Matrix(uint32_t r, uint32_t c) : arr_(std::make_unique<Ring[]>(r * c)), row_(r), col_(c), sz_(r * c) {}
		Matrix(Matrix&& v) noexcept = default;
		Matrix& operator=(Matrix&& v) noexcept = default;
		~Matrix() = default;
		Matrix(const Matrix& v) = delete;
		template <class T, class = std::enable_if_t<is_intermediate_v<Ring, T> ||
		                                            (std::is_same_v<T, base> && !std::is_same_v<T, Ring>)>>
		explicit Matrix(const Matrix<T>& v) : Matrix(v.row_, v.col_)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = Ring(v.arr_[i]);
		}
		Matrix& operator=(const Matrix& v) = default;
		[[nodiscard]] Ring* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Ring* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Ring& get(uint32_t r, uint32_t c) { return arr_[r * col_ + c]; }
		[[nodiscard]] const Ring& get(uint32_t r, uint32_t c) const { return arr_[r * col_ + c]; }
		[[nodiscard]] Ring& operator[](std::array<std::uint32_t, 2> i) { return get(i[0], i[1]); }
		[[nodiscard]] const Ring& operator[](std::array<std::uint32_t, 2> i) const { return get(i[0], i[1]); }
		[[nodiscard]] const uint32_t& size() const noexcept { return sz_; }
		[[nodiscard]] const uint32_t& row() const noexcept { return row_; }
		[[nodiscard]] const uint32_t& column() const noexcept { return col_; }
		[[nodiscard]] bool is_square() const noexcept { return row_ == col_; }
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const
		{
			Ring s{};
			for (uint32_t i = 0; i < sz_; ++i) s += arr_[i] * arr_[i];
			return s;
		}
		[[nodiscard]] Matrix clone() const
		{
			Matrix v(row_, col_);
			for (uint32_t i = 0; i < sz_; ++i) v.arr_[i] = arr_[i].clone();
			return v;
		}
		void swap(Matrix& other)
		{
			arr_.swap(other.arr_);
			std::swap(row_, other.row_);
			std::swap(col_, other.col_);
			std::swap(sz_, other.sz_);
		}
		[[nodiscard]] bool iszero() const noexcept
		{
			for (uint32_t i = 0; i < sz_; ++i)
				if (!arr_[i].iszero()) return false;
			return true;
		}
		void negate()
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = -arr_[i];
		}
		template <class T, class = std::enable_if_t<is_iaddable_v<Ring, T>>>
		Matrix& operator+=(const Matrix<T>& v)
		{
			assert(v.row_ == row_ && v.col_ == col_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] += v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_isubtractable_v<Ring, T>>>
		Matrix& operator-=(const Matrix<T>& v)
		{
			assert(v.row_ == row_ && v.col_ == col_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] -= v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_imultipliable_v<Ring, T>>>
		Matrix& operator*=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] *= r;
			return *this;
		}
		template <class T, class = std::enable_if_t<is_idividable_v<Ring, T>>>
		Matrix& operator/=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] /= r;
			return *this;
		}
		Vector<Ring> flatten()
		{
			Vector<Ring> v(arr_.release(), sz_);
			row_ = col_ = sz_ = 0;
			return v;
		}
		void transpose()
		{
			if (row_ == col_)
			{
				for (uint32_t i = 0; i < row_; ++i)
					for (uint32_t j = 0; j < i; ++j) get(i, j).swap(get(j, i));
			}
			else
			{
				auto t = clone();
				for (uint32_t i = 0; i < row_; ++i)
					for (uint32_t j = 0; j < col_; ++j) get(i, j) = std::move(t.get(j, i));
			}
		}
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix inverse() const
		{
			assert(row_ == col_);
			Matrix m = clone(), inv = Matrix::constant(base(1), row_);
			for (uint32_t j = 0; j < row_; ++j)
			{
				uint32_t p = j;
				auto max = Ring{};
				for (uint32_t i = p; i < row_; ++i)
				{
					auto abs = mpfr::abs(m.get(i, j));
					if (abs > max)
					{
						max = abs;
						p = i;
					}
				}
				inv.swap_row(p, j);
				m.swap_row(p, j, j);
				inv.multiply_row(j, 1 / m.get(j, j));
				m.multiply_row(j, 1 / m.get(j, j), j);
				for (uint32_t i = j + 1; i < row_; ++i)
				{
					inv.add_row(j, i, -m.get(i, j));
					m.add_row(j, i, -m.get(i, j), j);
				}
			}
			for (uint32_t j = row_ - 1; j < row_; --j)
				for (uint32_t r = j - 1; r < j; --r)
				{
					inv.add_row(j, r, -m.get(r, j));
					m.add_row(j, r, -m.get(r, j), j);
				}
			return inv;
		}
		void add_row(uint32_t f, uint32_t t, const Ring& x, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) get(t, c) += x * get(f, c);
		}
		void multiply_row(uint32_t r, const Ring& x, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) get(r, c) *= x;
		}
		void swap_row(uint32_t r1, uint32_t r2, uint32_t c0 = 0)
		{
			for (uint32_t c = c0; c < col_; ++c) get(r1, c).swap(get(r2, c));
		}
		static Matrix constant(const Ring& c, uint32_t n)
		{
			Matrix m(n, n);
			for (uint32_t i = 0; i < n; ++i) m.get(i, i) = c.clone();
			return m;
		}
		friend Matrix operator+(const Matrix& x) { return x.clone(); }
		friend Matrix operator-(const Matrix& x)
		{
			auto z = x.clone();
			z.negate();
			return z;
		}
		template <class R, class = std::enable_if_t<is_addable_v<Ring, R>>>
		friend Matrix<union_ring_t<Ring, R>> operator+(const Matrix& x, const Matrix<R>& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			Matrix<union_ring_t<Ring, R>> z(x.row_, x.col_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<is_subtractable_v<Ring, R>>>
		friend Matrix<union_ring_t<Ring, R>> operator-(const Matrix& x, const Matrix<R>& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_);
			Matrix<union_ring_t<Ring, R>> z(x.row_, x.col_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Matrix<union_ring_t<Ring, R>> operator*(const Matrix& x, const R& r)
		{
			Matrix<union_ring_t<Ring, R>> z(x.row_, x.col_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] * r;
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Matrix<union_ring_t<Ring, R>> operator*(const R& r, const Matrix& x)
		{
			return x * r;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_dividable_v<Ring, R>>>
		friend Matrix operator/(const Matrix& x, const R& r)
		{
			Matrix z(x.row_, x.col_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		template <class R, class = std::enable_if_t<is_multipliable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator*(const Matrix& x, const Vector<R>& y)
		{
			assert(x.col_ == y.size());
			Vector<union_ring_t<Ring, R>> z(x.row_);
			union_ring_t<Ring, R> s;
			for (uint32_t i = 0; i < x.row_; ++i)
			{
				s = {};
				for (uint32_t j = 0, p = i * x.col_; j < x.col_; ++j, ++p) s += x.arr_[p] * y[j];
				z[i] = std::move(s);
			}
			return z;
		}
		template <class R, class = std::enable_if_t<is_multipliable_v<Ring, R>>>
		friend Vector<union_ring_t<Ring, R>> operator*(const Vector<R>& x, const Matrix& y)
		{
			assert(x.size() == y.row_);
			Vector<union_ring_t<Ring, R>> z(y.col_);
			union_ring_t<Ring, R> s;
			for (uint32_t i = 0; i < y.col_; ++i)
			{
				s = {};
				for (uint32_t j = 0; j < y.row_; ++j) s += x[j] * y.get(j, i);
				z[i] = std::move(s);
			}
			return z;
		}
		template <class R, class = std::enable_if_t<is_multipliable_v<Ring, R>>>
		friend Matrix<union_ring_t<Ring, R>> operator*(const Matrix& x, const Matrix<R>& y)
		{
			assert(x.col_ == y.row_);
			Matrix<union_ring_t<Ring, R>> z(x.row_, y.col_);
			union_ring_t<Ring, R> s;
			for (uint32_t i = 0; i < x.row_; ++i)
			{
				for (uint32_t j = 0; j < y.col_; ++j)
				{
					s = {};
					for (uint32_t k = 0; k < x.col_; ++k) s += x.get(i, k) * y.get(k, j);
					z.get(i, j) = std::move(s);
				}
			}
			return z;
		}
		friend bool operator==(const Matrix& x, const Matrix& y)
		{
			if (x.row_ != y.row_ || x.col_ != y.col_) return false;
			for (uint32_t i = 0; i < x.sz_; ++i)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Matrix& x, const Matrix& y) { return !(x == y); }
		// calculate the cholesky decomposition L of positive definite matrix, by Cholesky–Banachiewicz algorithm.
		// this = L L^t and L is lower triangular.
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
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
					for (uint32_t k = 0; k < j; ++k) s += L.get(i, k) * L.get(j, k);
					s = get(i, j) - s;
					L.get(i, j) = i == j ? s.sqrt() : s / L.get(j, j);
				}
			}
			return L;
		}
		// calculate the inverse matrix of lower triangular matrix
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		[[nodiscard]] Matrix lower_triangular_inverse() const
		{
			// this must be lower triangular, i.e., get(i, j) = 0 for i < j
			assert(is_square());
			Matrix res(row_, row_);
			Ring s;
			for (uint32_t i = 0; i < row_; ++i)
			{
				res.get(i, i) = 1 / get(i, i);
				for (uint32_t j = 0; j < i; ++j)
				{
					s = {};
					for (uint32_t k = j; k < i; ++k) s += get(i, k) * res.get(k, j);
					res.get(i, j) = -s / get(i, i);
				}
			}
			return res;
		}
	};

	template <class Ring>
	class Tensor
	{
		std::unique_ptr<Ring[]> arr_;
		uint32_t row_, col_, len_, sz_;

	public:
		using base = typename base_ring<Ring>::type;
		using ring = Ring;
		using type = Tensor<Ring>;
		Tensor(uint32_t r, uint32_t c, uint32_t l)
		    : arr_(std::make_unique<Ring[]>(r * c * l)), row_(r), col_(c), len_(l), sz_(r * c * l)
		{
		}
		Tensor(Tensor&& v) noexcept = default;
		Tensor& operator=(Tensor&& v) noexcept = default;
		~Tensor() = default;
		Tensor(const Tensor& v) = delete;
		Tensor& operator=(const Tensor& v) = default;
		template <class T, class = std::enable_if_t<is_intermediate_v<Ring, T> ||
		                                            (std::is_same_v<T, base> && !std::is_same_v<T, Ring>)>>
		explicit Tensor(const Tensor<T>& v) : Tensor(v.row_, v.col_, v.len_)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = Ring(v.arr_[i]);
		}
		[[nodiscard]] Ring* get() noexcept { return arr_.get(); }
		[[nodiscard]] const Ring* get() const noexcept { return arr_.get(); }
		[[nodiscard]] Ring& get(uint32_t r, uint32_t c, uint32_t i) { return arr_[(r * col_ + c) * len_ + i]; }
		[[nodiscard]] const Ring& get(uint32_t r, uint32_t c, uint32_t i) const
		{
			return arr_[(r * col_ + c) * len_ + i];
		}
		[[nodiscard]] Ring& operator[](std::array<std::uint32_t, 3> i) { return get(i[0], i[1], i[2]); }
		[[nodiscard]] const Ring& operator[](std::array<std::uint32_t, 3> i) const { return get(i[0], i[1], i[2]); }
		[[nodiscard]] const uint32_t& size() const noexcept { return sz_; }
		[[nodiscard]] const uint32_t& row() const noexcept { return row_; }
		[[nodiscard]] const uint32_t& column() const noexcept { return col_; }
		[[nodiscard]] const uint32_t& length() const noexcept { return len_; }
		template <class = std::enable_if<is_mpfr_real_v<Ring>>>
		[[nodiscard]] Ring abs() const
		{
			return norm().sqrt();
		}
		[[nodiscard]] Ring norm() const
		{
			Ring s{};
			for (uint32_t i = 0; i < sz_; ++i) s += arr_[i] * arr_[i];
			return s;
		}
		[[nodiscard]] Tensor clone() const
		{
			Tensor v(row_, col_, len_);
			for (uint32_t i = 0; i < sz_; ++i) v.arr_[i] = arr_[i].clone();
			return v;
		}
		void swap(Tensor& other)
		{
			arr_.swap(other.arr_);
			std::swap(row_, other.row_);
			std::swap(col_, other.col_);
			std::swap(len_, other.len_);
			std::swap(sz_, other.sz_);
		}
		[[nodiscard]] bool iszero() const noexcept
		{
			for (uint32_t i = 0; i < sz_; ++i)
				if (!arr_[i].iszero()) return false;
			return true;
		}
		void negate()
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] = -arr_[i];
		}
		template <class T, class = std::enable_if_t<is_iaddable_v<Ring, T>>>
		Tensor& operator+=(const Tensor<T>& v)
		{
			assert(v.row_ == row_ && v.col_ == col_ && v.len_ == len_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] += v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_isubtractable_v<Ring, T>>>
		Tensor& operator-=(const Tensor<T>& v)
		{
			assert(v.row_ == row_ && v.col_ == col_ && v.len_ == len_);
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] -= v.arr_[i];
			return *this;
		}
		template <class T, class = std::enable_if_t<is_imultipliable_v<Ring, T>>>
		Tensor& operator*=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] *= r;
			return *this;
		}
		template <class T, class = std::enable_if_t<is_idividable_v<Ring, T>>>
		Tensor& operator/=(const T& r)
		{
			for (uint32_t i = 0; i < sz_; ++i) arr_[i] /= r;
			return *this;
		}
		Vector<Ring> flatten()
		{
			Vector<Ring> v(arr_.release(), sz_);
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
		template <class R, class = std::enable_if_t<is_addable_v<Ring, R>>>
		friend Tensor<union_ring_t<Ring, R>> operator+(const Tensor& x, const Tensor<R>& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_ && x.len_ == y.len_);
			Tensor z(x.row_, x.col_, x.len_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] + y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<is_subtractable_v<Ring, R>>>
		friend Tensor<union_ring_t<Ring, R>> operator-(const Tensor& x, const Tensor<R>& y)
		{
			assert(x.row_ == y.row_ && x.col_ == y.col_ && x.len_ == y.len_);
			Tensor z(x.row_, x.col_, x.len_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] - y.arr_[i];
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Tensor<union_ring_t<Ring, R>> operator*(const Tensor& x, const R& r)
		{
			Tensor z(x.row_, x.col_, x.len_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] * r;
			return z;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_multipliable_v<Ring, R>>>
		friend Tensor<union_ring_t<Ring, R>> operator*(const R& r, const Tensor& x)
		{
			return x * r;
		}
		template <class R, class = std::enable_if_t<!is_linear_space_v<R> && is_dividable_v<Ring, R>>>
		friend Tensor operator/(const Tensor& x, const R& r)
		{
			Tensor z(x.row_, x.col_, x.len_);
			for (uint32_t i = 0; i < x.sz_; ++i) z.arr_[i] = x.arr_[i] / r;
			return z;
		}
		friend bool operator==(const Tensor& x, const Tensor& y)
		{
			if (x.row_ != y.row_ || x.col_ != y.col_ || x.len_ != y.len_) return false;
			for (uint32_t i = 0; i < x.sz_; ++i)
				if (x.arr_[i] != y.arr_[i]) return false;
			return true;
		}
		friend bool operator!=(const Tensor& x, const Tensor& y) { return !(x == y); }
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
