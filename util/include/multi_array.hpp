#ifndef MULTI_ARRAY_HPP_
#define MULTI_ARRAY_HPP_

#include <array>             // for array
#include <cassert>           // for assert
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <vector>            // for vector

namespace util
{
	template <size_t rank_>
	constexpr size_t _product(const std::array<size_t, rank_>& a) noexcept
	{
		size_t p = 1;
		for (auto& x : a) p *= x;
		return p;
	}
	template <size_t rank_>
	constexpr std::array<size_t, rank_ + 1> _products(const std::array<size_t, rank_>& a) noexcept
	{
		std::array<size_t, rank_ + 1> p{};
		p[rank_] = 1;
		for (size_t i = rank_ - 1; i < rank_; --i) p[i] = p[i + 1] * a[i];
		return p;
	}
	using std::size_t;
	template <class T, size_t rank_, size_t fixed_>
	class multi_array_view;
	template <class T, size_t rank_, size_t fixed_>
	class multi_array_const_view;
	template <class T, size_t rank_>
	class multi_array
	{
		std::vector<T> buf;
		std::array<size_t, rank_> dims;
		std::array<size_t, rank_ + 1> dim_products;
		multi_array(const std::array<size_t, rank_>& dim, const std::array<size_t, rank_ + 1>& pdim)
		    : buf(pdim[0]), dims(dim), dim_products(pdim)
		{
		}
		multi_array(const T& val, const std::array<size_t, rank_>& dim, const std::array<size_t, rank_ + 1>& pdim)
		    : buf(pdim[0], val), dims(dim), dim_products(pdim)
		{
		}

	public:
		explicit multi_array(const std::array<size_t, rank_>& dim) : multi_array(dim, _products(dim)) {}
		multi_array(const T& val, const std::array<size_t, rank_>& dim) : multi_array(val, dim, _products(dim)) {}
		multi_array() : buf(1), dims{}, dim_products{}
		{
			for (size_t i = 0; i < rank_; ++i) dims[i] = 1;
			for (size_t i = 0; i <= rank_; ++i) dim_products[i] = 1;
		}
		void resize(const std::array<size_t, rank_>& dim)
		{
			dims = dim;
			dim_products = _products(dims);
			buf.resize(dim_products[0]);
		}
		void resize(const std::array<size_t, rank_>& dim, const T& val)
		{
			dims = dim;
			dim_products = _products(dims);
			buf.resize(dim_products[0], val);
		}
		multi_array_view<T, rank_, 0> get_view() { return multi_array_view<T, rank_, 0>(this); }
		multi_array_view<T, rank_, 1> operator[](size_t k) { return get_view()[k]; }
		multi_array_const_view<T, rank_, 0> get_view() const { return multi_array_const_view<T, rank_, 0>(*this); }
		multi_array_const_view<T, rank_, 1> operator[](size_t k) const { return get_view()[k]; }
		size_t size() const noexcept { return buf.size(); }
		T* data() noexcept { return buf.data(); }
		const T* data() const noexcept { return buf.data(); }
		typename std::vector<T>::const_iterator begin() const noexcept { return buf.cbegin(); }
		typename std::vector<T>::const_iterator end() const noexcept { return buf.cend(); }
		typename std::vector<T>::const_iterator cbegin() const noexcept { return buf.cbegin(); }
		typename std::vector<T>::const_iterator cend() const noexcept { return buf.cend(); }
		typename std::vector<T>::iterator begin() noexcept { return buf.begin(); }
		typename std::vector<T>::iterator end() noexcept { return buf.end(); }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_view;
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_const_view;
	};
	template <class T, size_t rank_, size_t fixed_>
	class multi_array_view
	{
		multi_array<T, rank_>* view;
		const size_t begin_;
		multi_array_view(const multi_array_view<T, rank_, fixed_ - 1>& a, size_t id)
		    : view(a.view), begin_(a.begin_ + id * a.view->dim_products[fixed_])
		{
			static_assert(fixed_ < rank_);
		}

	public:
		multi_array_view<T, rank_, fixed_ + 1> operator[](size_t k)
		{
			assert(k < size());
			return multi_array_view<T, rank_, fixed_ + 1>(*this, k);
		}
		size_t size() const noexcept { return view->dims[fixed_]; }
		size_t total_size() const noexcept { return view->dim_products[fixed_]; }
		auto begin() const noexcept { return view->begin() + begin_; }
		auto end() const noexcept { return begin() + view->dim_products[fixed_]; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_view;
		friend class multi_array<T, rank_>;
	};
	template <class T, size_t rank_>
	class multi_array_view<T, rank_, 0>
	{
		multi_array<T, rank_>* view;
		const size_t begin_ = 0;
		explicit multi_array_view(multi_array<T, rank_>* a) : view(a) {}
		multi_array_view<T, rank_, 1> operator[](size_t k)
		{
			assert(k < size());
			return multi_array_view<T, rank_, 1>(*this, k);
		}

	public:
		size_t size() const noexcept { return view->dims[0]; }
		size_t total_size() const noexcept { return view->dim_products[0]; }
		auto begin() const noexcept { return view->begin() + begin_; }
		auto end() const noexcept { return begin() + view->dim_products[0]; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_view;
		friend class multi_array<T, rank_>;
	};
	template <class T>
	class multi_array_view<T, 0, 0>
	{
		multi_array<T, 0>* view;
		explicit multi_array_view(multi_array<T, 0>* a) : view(a) {}

	public:
		T* operator->() const noexcept { return &view->buf[0]; }
		T& operator*() const noexcept { return view->buf[0]; }
		size_t total_size() const noexcept { return 1; }
		auto begin() const noexcept { return view->begin(); }
		auto end() const noexcept { return begin() + 1; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_view;
		friend class multi_array<T, 0>;
	};
	template <class T, size_t rank_>
	class multi_array_view<T, rank_, rank_>
	{
		multi_array<T, rank_>* view;
		const size_t begin_;
		multi_array_view(multi_array_view<T, rank_, rank_ - 1>& a, size_t id) : view(a.view), begin_(a.begin_ + id) {}

	public:
		T* operator->() const noexcept { return &view->buf[begin_]; }
		T& operator*() const noexcept { return view->buf[begin_]; }
		size_t total_size() const noexcept { return 1; }
		auto begin() const noexcept { return view->begin() + begin_; }
		auto end() const noexcept { return begin() + 1; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_view;
		friend class multi_array<T, rank_>;
	};
	template <class T, size_t rank_, size_t fixed_>
	class multi_array_const_view
	{
		const multi_array<T, rank_>& view;
		const size_t begin_, len_;

	public:
		multi_array_const_view(const multi_array_const_view<T, rank_, fixed_ - 1>& a, size_t id)
		    : view(a.view),
		      begin_(a.begin_ + id * (a.len_ / a.view.dims[fixed_ - 1])),
		      len_(a.len_ / a.view.dims[fixed_ - 1])
		{
			static_assert(fixed_ < rank_);
		}
		multi_array_const_view<T, rank_, fixed_ + 1> operator[](size_t k) const
		{
			assert(k < size());
			return multi_array_const_view<T, rank_, fixed_ + 1>(*this, k);
		}
		size_t size() const noexcept { return view.dims[fixed_]; }
		size_t total_size() const noexcept { return len_; }
		auto begin() const noexcept { return view.begin() + begin_; }
		auto end() const noexcept { return begin() + len_; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_const_view;
	};
	template <class T, size_t rank_>
	class multi_array_const_view<T, rank_, 0>
	{
		const multi_array<T, rank_>& view;
		const size_t begin_ = 0, len_;

	public:
		explicit multi_array_const_view(const multi_array<T, rank_>& a) : view(a), len_(a.buf.size()) {}
		multi_array_const_view<T, rank_, 1> operator[](size_t k) const
		{
			assert(k < size());
			return multi_array_const_view<T, rank_, 1>(*this, k);
		}
		size_t size() const noexcept { return view.dims[0]; }
		size_t total_size() const noexcept { return len_; }
		auto begin() const noexcept { return view.begin() + begin_; }
		auto end() const noexcept { return begin() + len_; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_const_view;
	};
	template <class T>
	class multi_array_const_view<T, 0, 0>
	{
		const multi_array<T, 0>& view;

	public:
		explicit multi_array_const_view(const multi_array<T, 0>& a) : view(a) {}
		const T* operator->() const noexcept { return &view.buf[0]; }
		const T& operator*() const noexcept { return view.buf[0]; }
		size_t total_size() const noexcept { return 1; }
		auto begin() const noexcept { return view.begin(); }
		auto end() const noexcept { return begin() + 1; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_const_view;
	};
	template <class T, size_t rank_>
	class multi_array_const_view<T, rank_, rank_>
	{
		const multi_array<T, rank_>& view;
		const size_t begin_;

	public:
		multi_array_const_view(const multi_array_const_view<T, rank_, rank_ - 1>& a, size_t id)
		    : view(a.view), begin_(a.begin_ + id)
		{
		}
		const T* operator->() const noexcept { return &view.buf[begin_]; }
		const T& operator*() const noexcept { return view.buf[begin_]; }
		size_t total_size() const noexcept { return 1; }
		auto begin() const noexcept { return view.begin() + begin_; }
		auto end() const noexcept { return begin() + 1; }
		template <class T2, size_t rank2_, size_t fixed2_>
		friend class multi_array_const_view;
	};
}  // namespace util

#endif  // MULTI_ARRAY_HPP_
