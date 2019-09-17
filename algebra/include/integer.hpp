#ifndef QBOOT_INTEGER_HPP_
#define QBOOT_INTEGER_HPP_

#include <istream>      // for basic_istream
#include <limits>       // for numeric_limits
#include <ostream>      // for basic_ostream
#include <stdexcept>    // for runtime_error
#include <string>       // for to_string, string_literals
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, is_integral_v, is_signed_v, conjunction_v, integral_constant
#include <utility>      // for move

#include "gmpxx.h"
#include "mpfr.h"

namespace mpfr
{
	// clang-tidy warns to use int32_t or int64_t instead of long int
	using _long = long int;            // NOLINT
	using _ulong = unsigned long int;  // NOLINT

	class integer;
	class rational;
	class real;

	template <class Tp>
	inline constexpr bool _is_mp =
	    std::is_same_v<Tp, integer> || std::is_same_v<Tp, rational> || std::is_same_v<Tp, real>;

	mpz_srcptr _take(const integer& x);
	mpq_srcptr _take(const rational& x);
	mpfr_srcptr _take(const real& x);

	// check all values in integral class I1 are included in integral class I2
	// and I1, I2 has the same signed property
	template <class I1, class I2>
	inline constexpr bool is_included_v = std::numeric_limits<I2>::min() <= std::numeric_limits<I1>::min() &&
	                                      std::numeric_limits<I1>::max() <= std::numeric_limits<I2>::max();

	template <class I1, class I2>
	struct is_included : std::integral_constant<bool, is_included_v<I1, I2>>
	{
	};

	template <class Tp>
	inline constexpr bool _ulong_convertible_v =
	    std::conjunction_v<std::is_integral<Tp>, std::is_unsigned<Tp>, is_included<Tp, _ulong>>;
	template <class Tp>
	inline constexpr bool _long_convertible_v =
	    std::conjunction_v<std::is_integral<Tp>, std::is_signed<Tp>, is_included<Tp, _long>>;

	template <class Tp>
	inline constexpr bool _mpz_is_other_operands = _long_convertible_v<Tp> || _ulong_convertible_v<Tp>;

	template <class Tp>
	struct _mp_ops;

	template <>
	struct _mp_ops<integer>
	{
		inline static void set(mpq_ptr rop, const integer& op) { mpq_set_z(rop, _take(op)); }
		inline static void set(mpfr_ptr rop, const integer& op, mpfr_rnd_t rnd) { mpfr_set_z(rop, _take(op), rnd); }
		inline static integer get(mpq_srcptr rop);
		inline static integer get(mpfr_srcptr rop, mpfr_rnd_t rnd);
		inline static void add(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			mpz_addmul(mpq_numref(rop), mpq_denref(op1), _take(op2));
		}
		inline static void add(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_add_z(rop, op1, _take(op2), rnd);
		}
		inline static void sub_a(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			mpz_submul(mpq_numref(rop), mpq_denref(op1), _take(op2));
		}
		inline static void sub_a(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_sub_z(rop, op1, _take(op2), rnd);
		}
		inline static void sub_b(mpq_ptr rop, const integer& op1, mpq_srcptr op2)
		{
			sub_a(rop, op2, op1);
			mpq_neg(rop, rop);
		}
		inline static void sub_b(mpfr_ptr rop, const integer& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_z_sub(rop, _take(op1), op2, rnd);
		}
		inline static void mul(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpz_set(mpq_denref(rop), mpq_denref(op1));
			mpz_mul(mpq_numref(rop), mpq_numref(op1), _take(op2));
		}
		inline static void mul(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_mul_z(rop, op1, _take(op2), rnd);
		}
		inline static void div_a(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpz_set(mpq_numref(rop), mpq_numref(op1));
			mpz_addmul(mpq_denref(rop), mpq_denref(op1), _take(op2));
		}
		inline static void div_a(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_div_z(rop, op1, _take(op2), rnd);
		}
		inline static void div_b(mpq_ptr rop, const integer& op1, mpq_srcptr op2)
		{
			div_a(rop, op2, op1);
			mpq_inv(rop, rop);
		}
		inline static void div_b(mpfr_ptr rop, const integer& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_ui_div(rop, 1, op2, rnd);
			mul(rop, rop, op1, rnd);
		}
		inline static int cmp(mpq_srcptr op1, const integer& op2) { return mpq_cmp_z(op1, _take(op2)); }
		inline static int cmp(mpfr_srcptr op1, const integer& op2) { return mpfr_cmp_z(op1, _take(op2)); }
	};

	template <>
	struct _mp_ops<rational>
	{
		inline static void set(mpfr_ptr rop, const rational& op, mpfr_rnd_t rnd) { mpfr_set_q(rop, _take(op), rnd); }
		inline static rational get(mpfr_srcptr rop, mpfr_rnd_t rnd);
		inline static void add(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_add_q(rop, op1, _take(op2), rnd);
		}
		inline static void sub_a(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_sub_q(rop, op1, _take(op2), rnd);
		}
		inline static void sub_b(mpfr_ptr rop, const rational& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_neg(rop, op2, rnd);
			add(rop, rop, op1, rnd);
		}
		inline static void mul(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_mul_q(rop, op1, _take(op2), rnd);
		}
		inline static void div_a(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_div_q(rop, op1, _take(op2), rnd);
		}
		inline static void div_b(mpfr_ptr rop, const rational& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_ui_div(rop, 1, op2, rnd);
			mul(rop, rop, op1, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const rational& op2) { return mpfr_cmp_q(op1, _take(op2)); }
	};

	template <>
	struct _mp_ops<_ulong>
	{
		inline static void init_set(mpz_ptr rop, _ulong op) { mpz_init_set_ui(rop, op); }
		inline static void set(mpz_ptr rop, _ulong op) { mpz_set_ui(rop, op); }
		inline static void set(mpq_ptr rop, _ulong op) { mpq_set_ui(rop, op, 1); }
		inline static void set(mpfr_ptr rop, _ulong op, mpfr_rnd_t rnd) { mpfr_set_ui(rop, op, rnd); }
		inline static _ulong get(mpz_srcptr op) { return mpz_get_ui(op); }
		inline static _ulong get(mpq_srcptr op);
		inline static _ulong get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_ui(op, rnd); }
		inline static void add(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_add_ui(rop, op1, op2); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static void sub_a(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_sub_ui(rop, op1, op2); }
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static void sub_b(mpz_ptr rop, _ulong op1, mpz_srcptr op2) { mpz_ui_sub(rop, op1, op2); }
		inline static int sub_b(mpfr_ptr rop, _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static void mul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_mul_ui(rop, op1, op2); }
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_div(rop, op1, op2, rnd);
		}
		// rop += op1 * op2
		inline static void addmul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_addmul_ui(rop, op1, op2); }
		// rop -= op1 * op2
		inline static void submul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_submul_ui(rop, op1, op2); }
		inline static int cmp(mpz_srcptr op1, _ulong op2) { return mpz_cmp_ui(op1, op2); }
		inline static int cmp(mpq_srcptr op1, _ulong op2) { return mpq_cmp_ui(op1, op2, 1); }
		inline static int cmp(mpfr_srcptr op1, _ulong op2) { return mpfr_cmp_ui(op1, op2); }
	};
	template <>
	struct _mp_ops<_long>
	{
		inline static void init_set(mpz_ptr rop, _long op) { mpz_init_set_si(rop, op); }
		inline static void set(mpz_ptr rop, _long op) { mpz_set_si(rop, op); }
		inline static void set(mpq_ptr rop, _long op) { mpq_set_si(rop, op, 1); }
		inline static void set(mpfr_ptr rop, _long op, mpfr_rnd_t rnd) { mpfr_set_si(rop, op, rnd); }
		inline static _long get(mpz_srcptr op) { return mpz_get_si(op); }
		inline static _long get(mpq_srcptr op);
		inline static _long get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_si(op, rnd); }
		inline static void add(mpz_ptr rop, mpz_srcptr op1, _long op2)
		{
			if (op2 >= 0)
				mpz_add_ui(rop, op1, _ulong(op2));
			else
				mpz_sub_ui(rop, op1, _ulong(-op2));
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_si(rop, op1, op2, rnd);
		}
		inline static void sub_a(mpz_ptr rop, mpz_srcptr op1, _long op2)
		{
			if (op2 >= 0)
				mpz_sub_ui(rop, op1, _ulong(op2));
			else
				mpz_add_ui(rop, op1, _ulong(-op2));
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static void sub_b(mpz_ptr rop, _long op1, mpz_srcptr op2)
		{
			if (op1 >= 0)
				mpz_ui_sub(rop, _ulong(op1), op2);
			else
			{
				mpz_add_ui(rop, op2, _ulong(-op1));
				mpz_neg(rop, rop);
			}
		}
		inline static int sub_b(mpfr_ptr rop, _long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static void mul(mpz_ptr rop, mpz_srcptr op1, _long op2) { mpz_mul_si(rop, op1, op2); }
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, _long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_div(rop, op1, op2, rnd);
		}
		// rop += op1 * op2
		inline static void addmul(mpz_ptr rop, mpz_srcptr op1, _long op2)
		{
			if (op2 >= 0)
				mpz_addmul_ui(rop, op1, _ulong(op2));
			else
				mpz_submul_ui(rop, op1, _ulong(-op2));
		}
		// rop -= op1 * op2
		inline static void submul(mpz_ptr rop, mpz_srcptr op1, _long op2)
		{
			if (op2 >= 0)
				mpz_submul_ui(rop, op1, _ulong(op2));
			else
				mpz_addmul_ui(rop, op1, _ulong(-op2));
		}
		inline static int cmp(mpz_srcptr op1, _long op2) { return mpz_cmp_si(op1, op2); }
		inline static int cmp(mpq_srcptr op1, _long op2) { return mpq_cmp_si(op1, op2, 1); }
		inline static int cmp(mpfr_srcptr op1, _long op2) { return mpfr_cmp_si(op1, op2); }
	};
	template <>
	struct _mp_ops<unsigned int> : _mp_ops<_ulong>
	{
	};
	template <>
	struct _mp_ops<int> : _mp_ops<_long>
	{
	};
	template <>
	struct _mp_ops<double>
	{
		inline static void init_set(mpz_ptr rop, double op) { mpz_init_set_d(rop, op); }
		inline static void set(mpz_ptr rop, double op) { mpz_set_d(rop, op); }
		inline static void set(mpq_ptr rop, double op) { mpq_set_d(rop, op); }
		inline static void set(mpfr_ptr rop, double op, mpfr_rnd_t rnd) { mpfr_set_d(rop, op, rnd); }
		inline static double get(mpz_srcptr op) { return mpz_get_d(op); }
		inline static double get(mpfr_srcptr op, mpfr_rnd_t rnd) { return mpfr_get_d(op, rnd); }
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, double op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_d(rop, op1, op2, rnd);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, double op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_d(rop, op1, op2, rnd);
		}
		inline static int sub_b(mpfr_ptr rop, double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_d_sub(rop, op1, op2, rnd);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, double op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_d(rop, op1, op2, rnd);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, double op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_d(rop, op1, op2, rnd);
		}
		inline static int div_b(mpfr_ptr rop, double op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_d_div(rop, op1, op2, rnd);
		}
		inline static int cmp(mpz_srcptr op1, double op2) { return mpz_cmp_d(op1, op2); }
		inline static int cmp(mpfr_srcptr op1, double op2) { return mpfr_cmp_d(op1, op2); }
	};

	// generate relational operators from _cmp (using ADL)

	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator==(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) == 0;
	}
	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator!=(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) != 0;
	}
	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator<(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) < 0;
	}
	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator>(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) > 0;
	}
	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator<=(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) <= 0;
	}
	template <class Tp1, class Tp2, class = std::enable_if_t<_is_mp<Tp1> || _is_mp<Tp2>>>
	inline bool operator>=(const Tp1& r1, const Tp2& r2)
	{
		return _cmp(r1, r2) >= 0;
	}

	class integer
	{
		friend class rational;
		friend class real;
		friend mpz_srcptr _take(const integer& x) { return x._x; }

	public:
		mpz_t _x;  // NOLINT

		integer() { mpz_init(_x); }
		integer(const integer& o) { mpz_init_set(_x, o._x); }
		integer(integer&& o) noexcept
		{
			mpz_init(_x);
			mpz_swap(_x, o._x);
		}
		integer& operator=(const integer& o) &
		{
			if (this != &o) mpz_set(_x, o._x);
			return *this;
		}
		integer& operator=(integer&& o) & noexcept
		{
			if (this != &o) mpz_swap(_x, o._x);
			return *this;
		}
		~integer() { mpz_clear(_x); }
		void swap(integer& o) & { mpz_swap(_x, o._x); }

		[[nodiscard]] integer clone() const { return *this; }
		[[nodiscard]] bool iszero() const { return mpz_sgn(_x) == 0; }
		void negate() & { mpz_neg(_x, _x); }

		[[nodiscard]] std::string str() const { return std::string(mpz_get_str(nullptr, 10, _x)); }

		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		explicit integer(T o)
		{
			_mp_ops<T>::init_set(_x, o);
		}
		explicit integer(std::string_view o)
		{
			using namespace std::string_literals;
			if (auto err = mpz_init_set_str(_x, o.data(), 10); err == -1)
			{
				throw std::runtime_error("in mpfr::integer(string_view):\n  invalid input format "s += o);
				mpz_clear(_x);
			}
		}
		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		integer& operator=(T o) &
		{
			_mp_ops<T>::set(_x, o);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
			return *this;
#pragma GCC diagnostic pop
		}
		integer& operator=(std::string_view o) &
		{
			using namespace std::string_literals;
			if (auto err = mpz_set_str(_x, o.data(), 10); err == -1)
				throw std::runtime_error("in mpfr::integer(string_view):\n  invalid input format "s += o);
			return *this;
		}
		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		explicit operator T() const
		{
			return _mp_ops<T>::get(_x);
		}

		// _cmp(a, b) returns the sign of a - b

		friend int _cmp(const integer& r1, const integer& r2) { return mpz_cmp(r1._x, r2._x); }
		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		friend int _cmp(const integer& r1, T r2)
		{
			return _mp_ops<T>::cmp(r1._x, r2);
		}
		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		friend int _cmp(T r1, const integer& r2)
		{
			return -_cmp(r2, r1);
		}

		template <class Tp>
		inline integer& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_add(_x, _x, o._x);
			else
				_mp_ops<Tp>::add(_x, _x, o);
			return *this;
		}

		template <class Tp>
		inline integer& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_sub(_x, _x, o._x);
			else
				_mp_ops<Tp>::sub_a(_x, _x, o);
			return *this;
		}

		template <class Tp>
		inline integer& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_mul(_x, _x, o._x);
			else
				_mp_ops<Tp>::mul(_x, _x, o);
			return *this;
		}

		friend integer mul(const integer& r1, const integer& r2) { return r1 * r2; }
		friend integer mul(integer&& r1, const integer& r2) { return std::move(r1) * r2; }
		friend integer mul(const integer& r1, integer&& r2) { return r1 * std::move(r2); }
		friend integer mul(integer&& r1, integer&& r2) { return std::move(r1) * r2; }
		friend integer mul_scalar(const integer& r1, const integer& r2) { return r1 * r2; }
		friend integer mul_scalar(integer&& r1, const integer& r2) { return std::move(r1) * r2; }
		friend integer mul_scalar(const integer& r1, integer&& r2) { return r1 * std::move(r2); }
		friend integer mul_scalar(integer&& r1, integer&& r2) { return std::move(r1) * r2; }

		friend integer operator+(const integer& a, const integer& b)
		{
			integer z;
			mpz_add(z._x, a._x, b._x);
			return z;
		}
		friend integer operator+(integer&& a, const integer& b) { return a += b; }
		friend integer operator+(const integer& a, integer&& b) { return b += a; }
		friend integer operator+(integer&& a, integer&& b) { return a += b; }

		friend integer operator-(const integer& a, const integer& b)
		{
			integer z;
			mpz_sub(z._x, a._x, b._x);
			return z;
		}
		friend integer operator-(integer&& a, const integer& b) { return a += b; }
		friend integer operator-(const integer& a, integer&& b)
		{
			mpz_sub(b._x, a._x, b._x);
			return std::move(b);
		}
		friend integer operator-(integer&& a, integer&& b) { return a += b; }

		friend integer operator*(const integer& a, const integer& b)
		{
			integer z;
			mpz_mul(z._x, a._x, b._x);
			return z;
		}
		friend integer operator*(integer&& a, const integer& b) { return a += b; }
		friend integer operator*(const integer& a, integer&& b) { return b += a; }
		friend integer operator*(integer&& a, integer&& b) { return a += b; }

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator+(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::add(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator+(integer&& r1, const Tp& r2)
		{
			return r1 += r2;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator+(const Tp& r1, const integer& r2)
		{
			return r2 + r1;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator+(const Tp& r1, integer&& r2)
		{
			return r2 += r1;
		}

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator-(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::sub_a(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator-(integer&& r1, const Tp& r2)
		{
			return r1 -= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator-(const Tp& r1, const integer& r2)
		{
			integer temp;
			_mp_ops<Tp>::sub_b(temp._x, r1, r2._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator-(const Tp& r1, integer&& r2)
		{
			_mp_ops<Tp>::sub_b(r2._x, r1, r2._x);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator*(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::mul(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator*(integer&& r1, const Tp& r2)
		{
			return r1 *= r2;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator*(const Tp& r1, const integer& r2)
		{
			return r2 * r1;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend inline integer operator*(const Tp& r1, integer&& r2)
		{
			return r2 *= r1;
		}

		inline integer operator+() const& { return *this; }
		inline integer operator+() && { return std::move(*this); }
		inline integer operator-() const&
		{
			integer temp;
			mpz_neg(temp._x, _x);
			return temp;
		}
		inline integer operator-() &&
		{
			mpz_neg(_x, _x);
			return std::move(*this);
		}
		friend inline integer mpfr::_mp_ops<integer>::get(mpq_srcptr rop);
		friend inline integer mpfr::_mp_ops<integer>::get(mpfr_srcptr rop, mpfr_rnd_t rnd);

		template <class Char, class Traits>
		friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& s, const integer& z)
		{
			return s << z._x;
		}
		template <class Char, class Traits>
		friend std::basic_istream<Char, Traits>& operator>>(std::basic_istream<Char, Traits>& s, integer& z)
		{
			return s >> z._x;
		}
	};
	inline _ulong mpfr::_mp_ops<_ulong>::get(mpq_srcptr op) { return _ulong(_mp_ops<integer>::get(op)); }
	inline _long mpfr::_mp_ops<_long>::get(mpq_srcptr op) { return _long(_mp_ops<integer>::get(op)); }
	inline integer mpfr::_mp_ops<integer>::get(mpfr_srcptr rop, mpfr_rnd_t rnd)
	{
		integer z;
		mpfr_get_z(z._x, rop, rnd);
		return z;
	}
}  // namespace mpfr

#endif  // QBOOT_INTEGER_HPP_
