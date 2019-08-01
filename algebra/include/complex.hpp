#ifndef COMPLEX_HPP_
#define COMPLEX_HPP_

#include <array>    // for array
#include <complex>  // for complex

#include "real.hpp"

namespace algebra
{
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	class complex;
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	class mpfr_complex
	{
	public:
		std::array<mpfr::real<_prec, _rnd>, 2> _val;
	};
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	class _complex_base : public mpfr_complex<_prec, _rnd>
	{
	public:
		using mpfr_complex<_prec, _rnd>::_val;
		using value_type = mpfr::real<_prec, _rnd>;

		_complex_base(const value_type& re, const value_type& im) : mpfr_complex<_prec, _rnd>{{re, im}} {}

		value_type real(const value_type& x) { return _val[0] = x; }

		value_type imag(const value_type& x) { return _val[1] = x; }

		[[nodiscard]] constexpr const value_type& real() const { return _val[0]; }

		[[nodiscard]] constexpr const value_type& imag() const { return _val[1]; }

	protected:
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		void _add(const complex<_prec2, _rnd2>& z)
		{
			_val[0] += value_type(z.real());
			_val[1] += value_type(z.imag());
		}

		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		void _sub(const complex<_prec2, _rnd2>& z)
		{
			_val[0] -= value_type(z.real());
			_val[1] -= value_type(z.imag());
		}

		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		void _mul(const complex<_prec2, _rnd2>& z)
		{
			auto re = value_type(z.real());
			auto im = value_type(z.imag());

			auto tmp = _val[0] * re - _val[1] * im;
			_val[1] = _val[0] * im + _val[1] * re;
			_val[0] = tmp;
		}

		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		void _div(const complex<_prec2, _rnd2>& z)
		{
			auto re = value_type(z.real());
			auto im = value_type(z.imag());
			int re_sg = mpfr::sgn(re), im_sg = mpfr::sgn(im);

			if (mpfr::isnan(re) || mpfr::isnan(im))
				_val[0] = _val[1] = mpfr::nan<_prec, _rnd>("");
			else if ((im_sg < 0 ? -im : im) < (re_sg < 0 ? -re : re))
			{
				value_type wr = im / re, wd = re + wr * im;

				if (mpfr::isnan(wd) || mpfr::iszero(wd))
					_val[0] = _val[1] = mpfr::nan<_prec, _rnd>("");
				else
				{
					auto tmp = (_val[0] + _val[1] * wr) / wd;
					_val[1] = (_val[1] - _val[0] * wr) / wd;
					_val[0] = tmp;
				}
			}
			else if (im_sg == 0)
				_val[0] = _val[1] = mpfr::nan<_prec, _rnd>("");
			else
			{
				value_type wr = re / im, wd = im + wr * re;

				if (mpfr::isnan(wd) || mpfr::iszero(wd))
					_val[0] = _val[1] = mpfr::nan<_prec, _rnd>("");
				else
				{
					auto tmp = (_val[0] * wr + _val[1]) / wd;
					_val[1] = (_val[1] * wr - _val[0]) / wd;
					_val[0] = tmp;
				}
			}
		}
	};
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	class complex : public _complex_base<_prec, _rnd>
	{
	public:
		using mpfr_complex<_prec, _rnd>::_val;
		using _complex_base<_prec, _rnd>::_add;
		using _complex_base<_prec, _rnd>::_mul;
		using _complex_base<_prec, _rnd>::_sub;
		using _complex_base<_prec, _rnd>::_div;
		using value_type = typename _complex_base<_prec, _rnd>::value_type;

		complex() : _complex_base<_prec, _rnd>(value_type(), value_type()) {}
		constexpr complex(const complex& z) = default;
		constexpr complex(complex&& z) = default;
		constexpr complex& operator=(const complex& z) = default;
		constexpr complex& operator=(complex&& z) = default;
		~complex() = default;

		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		explicit constexpr complex(const mpfr::real<_prec2, _rnd2>& re = mpfr::real<_prec2, _rnd2>(),
		                           const mpfr::real<_prec2, _rnd2>& im = mpfr::real<_prec2, _rnd2>())
		    : _complex_base<_prec, _rnd>(value_type(re), value_type(im))
		{
		}
		static constexpr auto nan() { return complex(mpfr::nan<_prec, _rnd>(""), mpfr::nan<_prec, _rnd>("")); }
		explicit constexpr complex(double x) : _complex_base<_prec, _rnd>(value_type(x), value_type()) {}
		explicit constexpr complex(long double x) : _complex_base<_prec, _rnd>(value_type(x), value_type()) {}
		explicit constexpr complex(int x) : _complex_base<_prec, _rnd>(value_type(x), value_type()) {}
		explicit constexpr complex(unsigned x) : _complex_base<_prec, _rnd>(value_type(x), value_type()) {}
		explicit constexpr complex(const mpfr_complex<_prec, _rnd>& z)
		    : _complex_base<_prec, _rnd>(z._val[0], z._val[1])
		{
		}
		template <class T>
		explicit constexpr complex(const std::complex<T>& z) : _complex_base<_prec, _rnd>(z.real(), z.imag())
		{
		}

		// assign operator
		complex& operator=(const value_type& x)
		{
			_val[0] = x;
			_val[1] = 0;
			return *this;
		}
		template <class T>
		complex& operator=(const std::complex<T>& z)
		{
			_val[0] = value_type(z.real());
			_val[1] = value_type(z.imag);
			return *this;
		}
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		complex& operator=(const complex<_prec2, _rnd2>& z)
		{
			_val[0] = value_type(z._val[0]);
			_val[1] = value_type(z._val[1]);
			return *this;
		}

		// add-assign operator
		complex& operator+=(const value_type& x)
		{
			_val[0] += x;
			return *this;
		}
		template <class T>
		complex& operator+=(const std::complex<T>& x)
		{
			_add(x);
			return *this;
		}
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		complex& operator+=(const complex<_prec2, _rnd2>& z)
		{
			_add(z);
			return *this;
		}

		// subtract-assign operator
		complex& operator-=(const value_type& x)
		{
			_val[0] -= x;
			return *this;
		}
		template <class T>
		complex& operator-=(const std::complex<T>& z)
		{
			_sub(z);
			return *this;
		}
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		complex& operator-=(const complex<_prec2, _rnd2>& z)
		{
			_sub(z);
			return *this;
		}

		// mutiply-assign operator
		complex& operator*=(const value_type& x)
		{
			_val[0] *= x;
			_val[1] *= x;
			return *this;
		}
		template <class T>
		complex& operator*=(const std::complex<T>& z)
		{
			_mul(z);
			return *this;
		}
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		complex& operator*=(const complex<_prec2, _rnd2>& z)
		{
			_mul(z);
			return *this;
		}

		// divide-assign operator
		complex& operator/=(const value_type& x)
		{
			_val[0] /= x;
			_val[1] /= x;
			return *this;
		}
		template <class T>
		complex& operator/=(const std::complex<T>& z)
		{
			_div(z);
			return *this;
		}
		template <mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
		complex& operator/=(const complex<_prec2, _rnd2>& z)
		{
			_div(z);
			return *this;
		}
	};
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline complex<_prec, _rnd> inv(const complex<_prec, _rnd>& z)
	{
		const auto& re = z.real();
		const auto& im = z.imag();
		int re_sg = mpfr::sgn(re), im_sg = mpfr::sgn(im);
		complex<_prec, _rnd> ans;
		if (mpfr::isnan(re) || mpfr::isnan(im))
			ans._val[0] = ans._val[1] = mpfr::nan<_prec, _rnd>("");
		else if ((im_sg < 0 ? -im : im) < (re_sg < 0 ? -re : re))
		{
			mpfr::real<_prec, _rnd> wr = im / re, wd = re + wr * im;
			if (mpfr::isnan(wd) || mpfr::iszero(wd))
				ans._val[0] = ans._val[1] = mpfr::nan<_prec, _rnd>("");
			else
			{
				ans._val[0] = 1 / wd;
				ans._val[1] = -wr / wd;
			}
		}
		else if (im_sg == 0)
			ans._val[0] = ans._val[1] = mpfr::nan<_prec, _rnd>("");
		else
		{
			mpfr::real<_prec, _rnd> wr = re / im, wd = im + wr * re;
			if (mpfr::isnan(wd) || mpfr::iszero(wd))
				ans._val[0] = ans._val[1] = mpfr::nan<_prec, _rnd>("");
			else
			{
				ans._val[0] = wr / wd;
				ans._val[1] = -1 / wd;
			}
		}
		return ans;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator==(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		return x.real() == y.real() && x.imag() == y.imag();
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator==(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		return x.real() == y && mpfr::iszero(x.imag());
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator==(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		return x == y.real() && mpfr::iszero(y.imag());
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator!=(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		return !(x == y);
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator!=(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		return !(x == y);
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr bool operator!=(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		return !(x == y);
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	constexpr auto iszero(const complex<_prec, _rnd>& z)
	{
		return mpfr::iszero(z.real()) && mpfr::iszero(z.imag());
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator+(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp += y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator+(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		return complex<res::prec, res::rnd>(x.real() + y, x.imag());
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator+(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp += y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator-(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp -= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator-(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		return complex<res::prec, res::rnd>(x.real() - y, x.imag());
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator-(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp -= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator*(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp *= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator*(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		return complex<res::prec, res::rnd>(x.real() * y, x.imag() * y);
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator*(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp *= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator/(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp /= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator/(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		return complex<res::prec, res::rnd>(x.real() / y, x.imag() / y);
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto operator/(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		complex<res::prec, res::rnd> tmp(x);
		tmp /= y;
		return tmp;
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline complex<_prec, _rnd> operator/(mpfr::mpfr_old_long x, const complex<_prec, _rnd>& y)
	{
		if (mpfr::isnan(y.real()) || mpfr::isnan(y.imag())) return complex<_prec, _rnd>::nan();
		if (x == 0) return complex<_prec, _rnd>();
		return mpfr::real<_prec, _rnd>(x) * inv(y);
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto conj(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(z.real(), -z.imag());
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto operator+(const complex<_prec, _rnd>& z)
	{
		return z;
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto operator-(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(-z.real(), -z.imag());
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto abs(const complex<_prec, _rnd>& z)
	{
		return mpfr::hypot(z.real(), z.imag());
	}
	// circ(x) := exp(I x) = cos(x) + I sin(x)
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto circ(const mpfr::real<_prec, _rnd>& x)
	{
		return complex<_prec, _rnd>(mpfr::cos(x), mpfr::sin(x));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto exp(const complex<_prec, _rnd>& z)
	{
		return circ(z.imag()) * mpfr::exp(z.real());
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto cos(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(mpfr::cos(z.real()) * mpfr::cosh(z.imag()),
		                            -mpfr::sin(z.real()) * mpfr::sinh(z.imag()));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto sin(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(mpfr::sin(z.real()) * mpfr::cosh(z.imag()),
		                            mpfr::cos(z.real()) * mpfr::sinh(z.imag()));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto cosh(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(mpfr::cosh(z.real()) * mpfr::cos(z.imag()),
		                            mpfr::sinh(z.real()) * mpfr::sin(z.imag()));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto sinh(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(mpfr::sinh(z.real()) * mpfr::cos(z.imag()),
		                            mpfr::cosh(z.real()) * mpfr::sin(z.imag()));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto log(const complex<_prec, _rnd>& z)
	{
		return complex<_prec, _rnd>(log(abs(z)), arg(z));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto log10(const complex<_prec, _rnd>& z)
	{
		static auto mul = mpfr::log10(mpfr::exp<_prec, _rnd>(1));
		return log(z) * mul;
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto _pow(const mpfr::real<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		if (mpfr::sgn(x) >= 0) return mpfr::pow(x, y);
		return exp(y * log(complex<res::prec, res::rnd>(x)));
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto pow(const complex<_prec1, _rnd1>& x, const mpfr::real<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		if (mpfr::iszero(x.imag())) return _pow(x.real(), y);
		return exp(y * log(complex<res::prec, res::rnd>(x)));
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto pow(const mpfr::real<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		if (mpfr::iszero(y.imag())) return _pow(x, y.real());
		if (mpfr::sgn(x) > 0) return exp(y * mpfr::log(mpfr::real<res::prec, res::rnd>(x)));
		return exp(y * log(complex<res::prec, res::rnd>(x)));
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	inline auto pow(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		if (mpfr::iszero(x.imag())) return pow(x.real(), y);
		return exp(y * log(complex<res::prec, res::rnd>(x)));
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	inline auto sqrt(const complex<_prec, _rnd>& z)
	{
		static auto half = 1 / mpfr::real<_prec, _rnd>(2);
		if (mpfr::iszero(z.imag()))
		{
			auto sg = mpfr::sgn(z.real());
			if (sg >= 0) return complex<_prec, _rnd>(mpfr::sqrt(z.real()));
			return complex<_prec, _rnd>(mpfr::real<_prec, _rnd>(0), mpfr::sqrt(-z.real()));
		}
		auto r = abs(z);
		auto sq = mpfr::sqrt((z.real() + r) * half);
		return complex<_prec, _rnd>(sq, z.imag() / sq * half);
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	constexpr auto real(const complex<_prec, _rnd>& z)
	{
		return z.real();
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	constexpr auto imag(const complex<_prec, _rnd>& z)
	{
		return z.imag();
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	constexpr auto arg(const complex<_prec, _rnd>& x)
	{
		return mpfr::atan2(x.imag(), x.real());
	}
	template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
	constexpr auto norm(const complex<_prec, _rnd>& z)
	{
		return real(z) * real(z) + imag(z) * imag(z);
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	constexpr auto polar(const mpfr::real<_prec1, _rnd1>& rho, const mpfr::real<_prec2, _rnd2>& theta)
	{
		using res = typename mpfr::result_type2<mpfr::real<_prec1, _rnd1>, mpfr::real<_prec2, _rnd2>, true>;
		return complex<res::prec, res::rnd>(rho * cos(mpfr::real<res::prec, res::rnd>(theta)),
		                                    rho * sin(mpfr::real<res::prec, res::rnd>(theta)));
	}
	template <mpfr::real_prec_t _prec1, mpfr::real_rnd_t _rnd1, mpfr::real_prec_t _prec2, mpfr::real_rnd_t _rnd2>
	auto inner_product(const complex<_prec1, _rnd1>& x, const complex<_prec2, _rnd2>& y)
	{
		return conj(x) * y;
	}
}  // namespace algebra

#endif  // COMPLEX_HPP_
