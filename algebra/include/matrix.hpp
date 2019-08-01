#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cassert>           // for assert
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <istream>           // for basic_istream
#include <map>               // for map, _Rb_tree_const_iterator, _Rb_tree...
#include <ostream>           // for basic_ostream
#include <string>            // for string
#include <tuple>             // for tuple
#include <vector>            // for vector

#include "complex.hpp"  // for iszero, complex, real_prec_t, real_rnd_t

namespace algebra
{
	using std::size_t;
	constexpr mpfr::real_prec_t _chop = 1000;
	constexpr mpfr::real_prec_t _prec = _chop + 30;
	constexpr mpfr::real_rnd_t _rnd = MPFR_REAL_CLASS_RND_DFLT;
	using Complex = complex<_prec, _rnd>;
	using Real = Complex::value_type;
	extern const Real epsilon;
	Real chop(const Real& x);
	Complex chop(const Complex& z);
	class SparseRow
	{
	private:
		std::map<size_t, Complex> elements_ = {};

	protected:
		void elements(const std::map<size_t, Complex>& e) { elements_ = e; }
		[[nodiscard]] const std::map<size_t, Complex>& elements() const noexcept { return elements_; }

	public:
		SparseRow() = default;
		~SparseRow() = default;
		SparseRow(const SparseRow&) = default;
		SparseRow(SparseRow&&) = default;
		SparseRow& operator=(const SparseRow&) = default;
		SparseRow& operator=(SparseRow&&) = default;
		[[nodiscard]] size_t size() const noexcept { return elements_.size(); }
		int set(size_t c, const Complex& v)
		{
			auto it = elements_.find(c);
			if (it != elements_.end())
				if (iszero(v))
				{
					elements_.erase(it);
					return -1;
				}
				else
					it->second = v;
			else if (!iszero(v))
			{
				elements_[c] = v;
				return 1;
			}
			return 0;
		}
		[[nodiscard]] Complex get(size_t c) const
		{
			auto it = elements_.find(c);
			return it != elements_.end() ? it->second : Complex(0);
		}
		// auto begin() noexcept { return elements_.begin(); }
		[[nodiscard]] auto begin() const noexcept { return elements_.begin(); }
		// auto end() noexcept { return elements_.end(); }
		[[nodiscard]] auto end() const noexcept { return elements_.end(); }
		SparseRow& conjugate()
		{
			for (auto& t : elements_) t.second = conj(t.second);
			return *this;
		}
		void add(const SparseRow& v, const Real& r);
		void add(const SparseRow& v, const Complex& r);
		SparseRow& operator*=(const Real& v);
		SparseRow& operator*=(const Complex& v);
		SparseRow& operator+=(const SparseRow& v);
		SparseRow& operator-=(const SparseRow& v);
		// Complex& operator[](size_t c);
		const Complex& operator[](size_t c) const;
		[[nodiscard]] std::string str() const;
		[[nodiscard]] std::string str(size_t C) const;

		friend class Matrix;
		friend Complex operator*(const SparseRow& r1, const SparseRow& r2)
		{
			auto sum = Complex();
			for (auto& it : r1)
			{
				auto f = r2.elements_.find(it.first);
				if (f != r2.elements_.end()) sum += it.second * f->second;
			}
			return sum;
		}
		friend Complex inner_product(const SparseRow& r1, const SparseRow& r2)
		{
			auto sum = Complex(0);
			for (auto& it : r1)
			{
				auto f = r2.elements_.find(it.first);
				if (f != r2.elements_.end()) sum += inner_product(it.second, f->second);
			}
			return sum;
		}
	};
	SparseRow conj(const SparseRow& v);

	class Matrix
	{
		std::vector<SparseRow> a;
		size_t row_ = 1, column_ = 1, nonzero_ = 0;
		void swap_row(size_t r1, size_t r2)
		{
			if (r1 == r2) return;
			auto t = a[r1];
			a[r1] = a[r2];
			a[r2] = t;
		}
		void multiply_row(size_t r, const Complex& x)
		{
			auto p = a[r].size();
			a[r] *= x;
			nonzero_ += a[r].size() - p;
		}
		void add_row(size_t f, size_t t, const Complex& x)
		{
			auto p = a[t].size();
			a[t].add(a[f], x);
			nonzero_ += a[t].size() - p;
		}

	public:
		Matrix() : a(1) {}
		Matrix(size_t r, size_t c) : a(r), row_(r), column_(c) {}
		[[nodiscard]] size_t row() const noexcept { return row_; }
		[[nodiscard]] size_t column() const noexcept { return column_; }
		[[nodiscard]] Complex at(size_t r, size_t c) const { return a[r].get(c); }
		void set_at(size_t r, size_t c, const Complex& v) { nonzero_ += size_t(a[r].set(c, v)); }
		void set_row(size_t r, const SparseRow& v)
		{
			nonzero_ -= a[r].size();
			a[r] = v;
			nonzero_ += a[r].size();
		}
		Matrix null_space();
		[[nodiscard]] Matrix orthogonalize() const;
		[[nodiscard]] size_t rank() const;
		void row_reduce();
		template <class InIt>
		static Matrix diagonal(InIt first, InIt last)
		{
			Matrix m(0, 0);
			size_t r = 0;
			while (first != last)
			{
				SparseRow v;
				v.set(r++, *(first++));
				m.a.push_back(v);
			}
			m.row_ = m.column_ = r;
			return m;
		}
		static Matrix constant(const Complex& v, size_t len);
		[[nodiscard]] Matrix inverse() const;
		[[nodiscard]] std::string str(size_t R1, size_t C) const;
		[[nodiscard]] std::string str2() const;
		static Matrix kronecker_product(const Matrix& x, const Matrix& y);
		static Matrix identity(size_t size);
		SparseRow& operator[](size_t r);
		const SparseRow& operator[](size_t r) const;
		const Complex& operator[](std::initializer_list<size_t> rc) const;
	};

	template <class Char, class Traits>
	std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os, const Matrix& g)
	{
		os << "{";
		for (size_t r = 0; r < g.row(); r++)
		{
			if (r > 0) os << ", ";
			os << '{';
			for (size_t c = 0; c < g.column(); c++)
			{
				if (c > 0) os << ", ";
				os << g.at(r, c);
			}
			os << '}';
		}
		return os << "}";
	}
	class Tensor2 : public SparseRow
	{
		size_t dim2_;

	public:
		[[nodiscard]] const size_t& dim2() const noexcept { return dim2_; }
		explicit Tensor2(size_t d2) : dim2_(d2) {}
		Tensor2() noexcept : Tensor2(1) {}
		Tensor2(size_t d2, const SparseRow& v) : SparseRow(v), dim2_(d2) {}
		~Tensor2() = default;
		Tensor2(const Tensor2& v) = default;
		Tensor2& operator=(const Tensor2& r)
		{
			if (this != &r)
			{
				dim2_ = r.dim2_;
				elements(r.elements());
			}
			return *this;
		}
		Tensor2(Tensor2&&) = default;
		Tensor2& operator=(Tensor2&&) = default;
		[[nodiscard]] Complex get(const std::tuple<size_t, size_t>& id) const
		{
			auto [a, b] = id;
			return get(a, b);
		}
		[[nodiscard]] Complex get(size_t a, size_t b) const { return SparseRow::get(dim2_ * a + b); }
		int set(size_t a, size_t b, const Complex& v) { return SparseRow::set(dim2_ * a + b, v); }
		int set(const std::tuple<size_t, size_t>& id, const Complex& v)
		{
			auto [a, b] = id;
			return set(a, b, v);
		}
		[[nodiscard]] std::tuple<size_t, size_t> to_index(size_t x) const
		{
			auto b = x % dim2_;
			x /= dim2_;
			return std::tuple(x, b);
		}
	};
	class Tensor3 : public SparseRow
	{
		size_t dim2_, dim3_;

	public:
		[[nodiscard]] const size_t& dim2() const { return dim2_; }
		[[nodiscard]] const size_t& dim3() const { return dim3_; }
		Tensor3() noexcept : Tensor3(1, 1) {}
		Tensor3(size_t d2, size_t d3) : dim2_(d2), dim3_(d3) {}
		Tensor3(size_t d2, size_t d3, const SparseRow& v) : SparseRow(v), dim2_(d2), dim3_(d3) {}
		~Tensor3() = default;
		Tensor3(const Tensor3& v) = default;
		Tensor3& operator=(const Tensor3& r)
		{
			if (this != &r)
			{
				dim2_ = r.dim2_;
				dim3_ = r.dim3_;
				elements(r.elements());
			}
			return *this;
		}
		Tensor3(Tensor3&&) = default;
		Tensor3& operator=(Tensor3&&) = default;
		[[nodiscard]] Complex get(const std::tuple<size_t, size_t, size_t>& id) const
		{
			auto [a, b, c] = id;
			return get(a, b, c);
		}
		[[nodiscard]] Complex get(size_t a, size_t b, size_t c) const
		{
			return SparseRow::get((dim2_ * a + b) * dim3_ + c);
		}
		int set(const std::tuple<size_t, size_t, size_t>& id, const Complex& v)
		{
			auto [a, b, c] = id;
			return set(a, b, c, v);
		}
		int set(size_t a, size_t b, size_t c, const Complex& v)
		{
			return SparseRow::set((dim2_ * a + b) * dim3_ + c, v);
		}
		[[nodiscard]] std::tuple<size_t, size_t, size_t> to_index(size_t x) const
		{
			auto c = x % dim3_;
			x /= dim3_;
			auto b = x % dim2_;
			x /= dim2_;
			return std::tuple(x, b, c);
		}
	};
	class Tensor4 : public SparseRow
	{
		size_t dim2_, dim3_, dim4_;

	public:
		[[nodiscard]] const size_t& dim2() const { return dim2_; }
		[[nodiscard]] const size_t& dim3() const { return dim3_; }
		[[nodiscard]] const size_t& dim4() const { return dim4_; }
		Tensor4() noexcept : Tensor4(1, 1, 1) {}
		Tensor4(size_t d2, size_t d3, size_t d4) : dim2_(d2), dim3_(d3), dim4_(d4) {}
		Tensor4(size_t d2, size_t d3, size_t d4, const SparseRow& v) : SparseRow(v), dim2_(d2), dim3_(d3), dim4_(d4) {}
		~Tensor4() = default;
		Tensor4(const Tensor4& v) = default;
		Tensor4& operator=(const Tensor4& r)
		{
			if (this != &r)
			{
				dim2_ = r.dim2_;
				dim3_ = r.dim3_;
				dim4_ = r.dim4_;
				elements(r.elements());
			}
			return *this;
		}
		Tensor4(Tensor4&&) = default;
		Tensor4& operator=(Tensor4&&) = default;
		[[nodiscard]] Complex get(const std::tuple<size_t, size_t, size_t, size_t>& id) const
		{
			auto [a, b, c, d] = id;
			return get(a, b, c, d);
		}
		[[nodiscard]] Complex get(size_t a, size_t b, size_t c, size_t d) const
		{
			return SparseRow::get(((dim2_ * a + b) * dim3_ + c) * dim4_ + d);
		}
		int set(const std::tuple<size_t, size_t, size_t, size_t>& id, const Complex& v)
		{
			auto [a, b, c, d] = id;
			return set(a, b, c, d, v);
		}
		int set(size_t a, size_t b, size_t c, size_t d, const Complex& v)
		{
			return SparseRow::set(((dim2_ * a + b) * dim3_ + c) * dim4_ + d, v);
		}
		[[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> to_index(size_t x) const
		{
			auto d = x % dim4_;
			x /= dim4_;
			auto c = x % dim3_;
			x /= dim3_;
			auto b = x % dim2_;
			x /= dim2_;
			return std::tuple(x, b, c, d);
		}
	};
	class Vector : public SparseRow
	{
		size_t dim_;

	public:
		[[nodiscard]] const size_t& dim() const { return dim_; }
		Vector() noexcept : Vector(1) {}
		explicit Vector(size_t d) : dim_(d) {}
		Vector(size_t d, const SparseRow& v) : SparseRow(v), dim_(d) {}
		Vector(const Vector& v) = default;
		~Vector() = default;
		Vector& operator=(const Vector& r)
		{
			if (this != &r)
			{
				dim_ = r.dim_;
				elements(r.elements());
			}
			return *this;
		}
		Vector(Vector&&) = default;
		Vector& operator=(Vector&&) = default;
		static Vector kronecker_product(const Vector& x, const Vector& y);
	};
}  // namespace algebra

#endif  // MATRIX_HPP_
