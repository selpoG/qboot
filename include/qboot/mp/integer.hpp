#ifndef QBOOT_MP_INTEGER_HPP_
#define QBOOT_MP_INTEGER_HPP_

#include <istream>      // for basic_istream
#include <limits>       // for numeric_limits
#include <optional>     // for optional
#include <ostream>      // for basic_ostream
#include <stdexcept>    // for runtime_error
#include <string>       // for to_string, string_literals
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if_t, is_same_v, is_integral_v, is_signed_v, conjunction_v, integral_constant
#include <utility>      // for move

#include "gmpxx.h"
#include "mpfr.h"

namespace qboot::mp
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

	inline bool _is_even(mpz_srcptr p) { return mpz_even_p(p) != 0; }           // NOLINT
	inline bool _is_odd(mpz_srcptr p) { return mpz_odd_p(p) != 0; }             // NOLINT
	inline int _cmp_ui(mpz_srcptr p, _ulong o) { return mpz_cmp_ui(p, o); }     // NOLINT
	inline int _cmp_ui(mpq_srcptr p, _ulong o) { return mpq_cmp_ui(p, o, 1); }  // NOLINT
	inline int _cmp_si(mpz_srcptr p, _long o) { return mpz_cmp_si(p, o); }      // NOLINT
	inline int _cmp_si(mpq_srcptr p, _long o) { return mpq_cmp_si(p, o, 1); }   // NOLINT

	template <class Tp>
	inline constexpr bool _mpz_is_other_operands = _long_convertible_v<Tp> || _ulong_convertible_v<Tp>;

	template <class Tp>
	struct _mp_ops;

	template <>
	struct _mp_ops<integer>
	{
		friend class integer;
		inline static void set(mpq_ptr rop, const integer& op) { mpq_set_z(rop, data(op)); }
		inline static void set(mpq_ptr rop, const integer& numop, const integer& denop)
		{
			mpz_set(mpq_numref(rop), data(numop));
			mpz_set(mpq_denref(rop), data(denop));
		}
		inline static void set(mpfr_ptr rop, const integer& op, mpfr_rnd_t rnd) { mpfr_set_z(rop, data(op), rnd); }
		inline static integer ceil(mpq_srcptr rop);
		inline static integer truncate(mpq_srcptr rop);
		inline static integer get(mpq_srcptr rop);
		inline static integer get(mpfr_srcptr rop, mpfr_rnd_t rnd);
		inline static void add(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			addmul(mpq_numref(rop), mpq_denref(op1), op2);
		}
		inline static void add(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_add_z(rop, op1, data(op2), rnd);
		}
		inline static void sub_a(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			submul(mpq_numref(rop), mpq_denref(op1), op2);
		}
		inline static void sub_a(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_sub_z(rop, op1, data(op2), rnd);
		}
		inline static void sub_b(mpq_ptr rop, const integer& op1, mpq_srcptr op2)
		{
			sub_a(rop, op2, op1);
			mpq_neg(rop, rop);
		}
		inline static void sub_b(mpfr_ptr rop, const integer& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_z_sub(rop, data(op1), op2, rnd);
		}
		inline static void mul(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpz_set(mpq_denref(rop), mpq_denref(op1));
			mpz_mul(mpq_numref(rop), mpq_numref(op1), data(op2));
		}
		inline static void mul(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_mul_z(rop, op1, data(op2), rnd);
		}
		inline static void div_a(mpq_ptr rop, mpq_srcptr op1, const integer& op2)
		{
			if (rop != op1) mpz_set(mpq_numref(rop), mpq_numref(op1));
			mpz_mul(mpq_denref(rop), mpq_denref(op1), data(op2));
		}
		inline static void div_a(mpfr_ptr rop, mpfr_srcptr op1, const integer& op2, mpfr_rnd_t rnd)
		{
			mpfr_div_z(rop, op1, data(op2), rnd);
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
		inline static void addmul(mpz_ptr rop, mpz_srcptr op1, const integer& op2) { mpz_addmul(rop, op1, data(op2)); }
		inline static void submul(mpz_ptr rop, mpz_srcptr op1, const integer& op2) { mpz_submul(rop, op1, data(op2)); }
		inline static int cmp(mpq_srcptr op1, const integer& op2) { return mpq_cmp_z(op1, data(op2)); }
		inline static int cmp(mpfr_srcptr op1, const integer& op2) { return mpfr_cmp_z(op1, data(op2)); }

	private:
		inline static mpz_srcptr data(const integer& op);
	};

	template <>
	struct _mp_ops<rational>
	{
		friend class rational;
		inline static void set(mpfr_ptr rop, const rational& op, mpfr_rnd_t rnd) { mpfr_set_q(rop, data(op), rnd); }
		inline static rational get(mpfr_srcptr rop, mpfr_rnd_t rnd);
		inline static void add(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_add_q(rop, op1, data(op2), rnd);
		}
		inline static void sub_a(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_sub_q(rop, op1, data(op2), rnd);
		}
		inline static void sub_b(mpfr_ptr rop, const rational& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			sub_a(rop, op2, op1, rnd);
			mpfr_neg(rop, rop, rnd);
		}
		inline static void mul(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_mul_q(rop, op1, data(op2), rnd);
		}
		inline static void div_a(mpfr_ptr rop, mpfr_srcptr op1, const rational& op2, mpfr_rnd_t rnd)
		{
			mpfr_div_q(rop, op1, data(op2), rnd);
		}
		inline static void div_b(mpfr_ptr rop, const rational& op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			mpfr_ui_div(rop, 1, op2, rnd);
			mul(rop, rop, op1, rnd);
		}
		inline static int cmp(mpfr_srcptr op1, const rational& op2) { return mpfr_cmp_q(op1, data(op2)); }

	private:
		inline static mpq_srcptr data(const rational& op);
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
		inline static void add(mpq_ptr rop, mpq_srcptr op1, _ulong op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			addmul(mpq_numref(rop), mpq_denref(op1), op2);
		}
		inline static int add(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_add_ui(rop, op1, op2, rnd);
		}
		inline static void sub_a(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_sub_ui(rop, op1, op2); }
		inline static void sub_a(mpq_ptr rop, mpq_srcptr op1, _ulong op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			submul(mpq_numref(rop), mpq_denref(op1), op2);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_ui(rop, op1, op2, rnd);
		}
		inline static void sub_b(mpz_ptr rop, _ulong op1, mpz_srcptr op2) { mpz_ui_sub(rop, op1, op2); }
		inline static void sub_b(mpq_ptr rop, _ulong op1, mpq_srcptr op2)
		{
			sub_a(rop, op2, op1);
			mpq_neg(rop, rop);
		}
		inline static int sub_b(mpfr_ptr rop, _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_sub(rop, op1, op2, rnd);
		}
		inline static void mul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_mul_ui(rop, op1, op2); }
		inline static void mul(mpq_ptr rop, mpq_srcptr op1, _ulong op2)
		{
			if (rop != op1) mpz_set(mpq_denref(rop), mpq_denref(op1));
			mul(mpq_numref(rop), mpq_numref(op1), op2);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_ui(rop, op1, op2, rnd);
		}
		inline static void div_a(mpq_ptr rop, mpq_srcptr op1, _ulong op2)
		{
			if (rop != op1) mpz_set(mpq_numref(rop), mpq_numref(op1));
			mul(mpq_denref(rop), mpq_denref(op1), op2);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, _ulong op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_ui(rop, op1, op2, rnd);
		}
		inline static void div_b(mpq_ptr rop, _ulong op1, mpq_srcptr op2)
		{
			div_a(rop, op2, op1);
			mpq_inv(rop, rop);
		}
		inline static int div_b(mpfr_ptr rop, _ulong op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_ui_div(rop, op1, op2, rnd);
		}
		// rop += op1 * op2
		inline static void addmul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_addmul_ui(rop, op1, op2); }
		// rop -= op1 * op2
		inline static void submul(mpz_ptr rop, mpz_srcptr op1, _ulong op2) { mpz_submul_ui(rop, op1, op2); }
		inline static int cmp(mpz_srcptr op1, _ulong op2) { return _cmp_ui(op1, op2); }
		inline static int cmp(mpq_srcptr op1, _ulong op2) { return _cmp_ui(op1, op2); }
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
		inline static void add(mpq_ptr rop, mpq_srcptr op1, _long op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			addmul(mpq_numref(rop), mpq_denref(op1), op2);
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
		inline static void sub_a(mpq_ptr rop, mpq_srcptr op1, _long op2)
		{
			if (rop != op1) mpq_set(rop, op1);
			submul(mpq_numref(rop), mpq_denref(op1), op2);
		}
		inline static int sub_a(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_sub_si(rop, op1, op2, rnd);
		}
		inline static void sub_b(mpz_ptr rop, _long op1, mpz_srcptr op2)
		{
			sub_a(rop, op2, op1);
			mpz_neg(rop, rop);
		}
		inline static void sub_b(mpq_ptr rop, _long op1, mpq_srcptr op2)
		{
			sub_a(rop, op2, op1);
			mpq_neg(rop, rop);
		}
		inline static int sub_b(mpfr_ptr rop, _long op1, mpfr_srcptr op2, mpfr_rnd_t rnd)
		{
			return mpfr_si_sub(rop, op1, op2, rnd);
		}
		inline static void mul(mpz_ptr rop, mpz_srcptr op1, _long op2) { mpz_mul_si(rop, op1, op2); }
		inline static void mul(mpq_ptr rop, mpq_srcptr op1, _long op2)
		{
			if (rop != op1) mpz_set(mpq_denref(rop), mpq_denref(op1));
			mul(mpq_numref(rop), mpq_numref(op1), op2);
		}
		inline static int mul(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_mul_si(rop, op1, op2, rnd);
		}
		inline static void div_a(mpq_ptr rop, mpq_srcptr op1, _long op2)
		{
			if (rop != op1) mpz_set(mpq_numref(rop), mpq_numref(op1));
			mul(mpq_denref(rop), mpq_denref(op1), op2);
		}
		inline static int div_a(mpfr_ptr rop, mpfr_srcptr op1, _long op2, mpfr_rnd_t rnd)
		{
			return mpfr_div_si(rop, op1, op2, rnd);
		}
		inline static void div_b(mpq_ptr rop, _long op1, mpq_srcptr op2)
		{
			div_a(rop, op2, op1);
			mpq_inv(rop, rop);
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
		inline static int cmp(mpz_srcptr op1, _long op2) { return _cmp_si(op1, op2); }
		inline static int cmp(mpq_srcptr op1, _long op2) { return _cmp_si(op1, op2); }
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

	inline std::optional<rational> parse(std::string_view str);
	inline std::optional<rational> _parse_mantisa(std::string_view str);

	inline void _reset(mpz_ptr x)
	{
		if (x->_mp_d != nullptr) mpz_clear(x);
		x->_mp_d = nullptr;
	}

	class integer
	{
		mpz_t _x;
		void reset() { qboot::mp::_reset(_x); }

	public:
		void _reset() && { reset(); }
		integer() { mpz_init(_x); }
		integer(const integer& o) { mpz_init_set(_x, o._x); }
		integer(integer&& o) noexcept
		{
			_x->_mp_d = nullptr;
			mpz_swap(_x, o._x);
		}
		integer& operator=(const integer& o) &
		{
			if (this != &o) mpz_set(_x, o._x);
			return *this;
		}
		integer& operator=(integer&& o) & noexcept
		{
			if (this != &o)
			{
				mpz_swap(_x, o._x);
				o.reset();
			}
			return *this;
		}
		~integer()
		{
			if (_x->_mp_d != nullptr) mpz_clear(_x);
		}
		void swap(integer& o) & { mpz_swap(_x, o._x); }

		explicit integer(mpz_srcptr o) { mpz_init_set(_x, o); }

		[[nodiscard]] integer clone() const { return *this; }
		[[nodiscard]] bool iszero() const { return mpz_sgn(_x) == 0; }
		void negate() & { mpz_neg(_x, _x); }
		[[nodiscard]] bool iseven() const { return _is_even(_x); }
		[[nodiscard]] bool isodd() const { return _is_odd(_x); }
		[[nodiscard]] bool divisible_by(const integer& o) const { return mpz_divisible_p(_x, o._x) != 0; }
		[[nodiscard]] bool divisible_by(_ulong o) const { return mpz_divisible_ui_p(_x, o) != 0; }

		template <class T>
		[[nodiscard]] integer eval([[maybe_unused]] const T& x) const
		{
			return *this;
		}

		[[nodiscard]] std::string str() const
		{
			std::string s(mpz_sizeinbase(_x, 10) + 2, 0);
			mpz_get_str(s.data(), 10, _x);
			s.resize(std::strlen(s.data()));
			return s;
		}
		static std::optional<integer> _parse(std::string_view str)
		{
			std::string s(str);
			integer x;
			if (mpz_set_str(x._x, s.data(), 10) == -1) return {};
			return {std::move(x)};
		}

		template <class T, class = std::enable_if_t<_mpz_is_other_operands<T> || std::is_same_v<T, double>>>
		explicit integer(T o)
		{
			_mp_ops<T>::init_set(_x, o);
		}
		explicit integer(std::string_view o)
		{
			std::string s(o);
			using namespace std::string_literals;
			if (mpz_init_set_str(_x, s.data(), 10) == -1)
			{
				mpz_clear(_x);
				throw std::runtime_error("in qboot::mp::integer(string_view):\n  invalid input format "s += o);
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
			std::string s(o);
			using namespace std::string_literals;
			if (mpz_set_str(_x, s.data(), 10) == -1)
				throw std::runtime_error("in qboot::mp::integer(string_view):\n  invalid input format "s += o);
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
		integer& operator+=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_add(_x, _x, o._x);
			else
				_mp_ops<Tp>::add(_x, _x, o);
			return *this;
		}

		template <class Tp>
		integer& operator-=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_sub(_x, _x, o._x);
			else
				_mp_ops<Tp>::sub_a(_x, _x, o);
			return *this;
		}

		template <class Tp>
		integer& operator*=(const Tp& o) &
		{
			if constexpr (std::is_same_v<Tp, integer>)
				mpz_mul(_x, _x, o._x);
			else
				_mp_ops<Tp>::mul(_x, _x, o);
			return *this;
		}

		integer& operator/=(const integer& o) &
		{
			mpz_fdiv_q(_x, _x, o._x);
			return *this;
		}
		integer& operator%=(const integer& o) &
		{
			mpz_fdiv_r(_x, _x, o._x);
			return *this;
		}
		integer& operator/=(_ulong o) &
		{
			mpz_fdiv_q_ui(_x, _x, o);
			return *this;
		}
		integer& operator%=(_ulong o) &
		{
			mpz_fdiv_r_ui(_x, _x, o);
			return *this;
		}

		friend integer mul(const integer& r1, const integer& r2) { return r1 * r2; }
		friend integer mul(integer&& r1, const integer& r2) { return std::move(r1) * r2; }
		friend integer mul(const integer& r1, integer&& r2) { return r1 * std::move(r2); }
		friend integer mul(integer&& r1, integer&& r2) { return std::move(r1) * std::move(r2); }
		friend integer mul_scalar(const integer& r1, const integer& r2) { return r1 * r2; }
		friend integer mul_scalar(integer&& r1, const integer& r2) { return std::move(r1) * r2; }
		friend integer mul_scalar(const integer& r1, integer&& r2) { return r1 * std::move(r2); }
		friend integer mul_scalar(integer&& r1, integer&& r2) { return std::move(r1) * std::move(r2); }

		friend integer operator+(const integer& a, const integer& b)
		{
			integer z;
			mpz_add(z._x, a._x, b._x);
			return z;
		}
		friend integer operator+(integer&& a, const integer& b) { return std::move(a += b); }
		friend integer operator+(const integer& a, integer&& b) { return std::move(b += a); }
		friend integer operator+(integer&& a, integer&& b)
		{
			a += b;
			b.reset();
			return std::move(a);
		}

		friend integer operator-(const integer& a, const integer& b)
		{
			integer z;
			mpz_sub(z._x, a._x, b._x);
			return z;
		}
		friend integer operator-(integer&& a, const integer& b) { return std::move(a -= b); }
		friend integer operator-(const integer& a, integer&& b)
		{
			mpz_sub(b._x, a._x, b._x);
			return std::move(b);
		}
		friend integer operator-(integer&& a, integer&& b)
		{
			a -= b;
			b.reset();
			return std::move(a);
		}

		friend integer operator*(const integer& a, const integer& b)
		{
			integer z;
			mpz_mul(z._x, a._x, b._x);
			return z;
		}
		friend integer operator*(integer&& a, const integer& b) { return std::move(a *= b); }
		friend integer operator*(const integer& a, integer&& b) { return std::move(b *= a); }
		friend integer operator*(integer&& a, integer&& b)
		{
			a *= b;
			b.reset();
			return std::move(a);
		}

		friend integer operator/(const integer& a, const integer& b)
		{
			integer z;
			mpz_fdiv_q(z._x, a._x, b._x);
			return z;
		}
		friend integer operator/(integer&& a, const integer& b) { return std::move(a /= b); }
		friend integer operator/(const integer& a, integer&& b)
		{
			mpz_fdiv_q(b._x, a._x, b._x);
			return std::move(b);
		}
		friend integer operator/(integer&& a, integer&& b)
		{
			a /= b;
			b.reset();
			return std::move(a);
		}

		friend integer operator%(const integer& a, const integer& b)
		{
			integer z;
			mpz_fdiv_r(z._x, a._x, b._x);
			return z;
		}
		friend integer operator%(integer&& a, const integer& b) { return std::move(a %= b); }
		friend integer operator%(const integer& a, integer&& b)
		{
			mpz_fdiv_r(b._x, a._x, b._x);
			return std::move(b);
		}
		friend integer operator%(integer&& a, integer&& b)
		{
			a %= b;
			b.reset();
			return std::move(a);
		}

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator+(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::add(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator+(integer&& r1, const Tp& r2)
		{
			return std::move(r1 += r2);
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator+(const Tp& r1, const integer& r2)
		{
			return r2 + r1;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator+(const Tp& r1, integer&& r2)
		{
			return std::move(r2 += r1);
		}

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator-(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::sub_a(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator-(integer&& r1, const Tp& r2)
		{
			return std::move(r1 -= r2);
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator-(const Tp& r1, const integer& r2)
		{
			integer temp;
			_mp_ops<Tp>::sub_b(temp._x, r1, r2._x);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator-(const Tp& r1, integer&& r2)
		{
			_mp_ops<Tp>::sub_b(r2._x, r1, r2._x);
			return std::move(r2);
		}

		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator*(const integer& r1, const Tp& r2)
		{
			integer temp;
			_mp_ops<Tp>::mul(temp._x, r1._x, r2);
			return temp;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator*(integer&& r1, const Tp& r2)
		{
			return std::move(r1 *= r2);
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator*(const Tp& r1, const integer& r2)
		{
			return r2 * r1;
		}
		template <class Tp, class = std::enable_if_t<_mpz_is_other_operands<Tp>>>
		friend integer operator*(const Tp& r1, integer&& r2)
		{
			return std::move(r2 *= r1);
		}

		friend integer operator/(const integer& r1, _ulong r2)
		{
			integer temp;
			mpz_fdiv_q_ui(temp._x, r1._x, r2);
			return temp;
		}
		friend integer operator/(integer&& r1, _ulong r2) { return std::move(r1 /= r2); }
		friend integer operator/(_ulong r1, const integer& r2)
		{
			if (r2 < 0) return r1 / -r2;
			if (_cmp_ui(r2._x, r1) > 0) return integer(0);
			return integer(r1 / _ulong(r2));
		}

		friend _ulong operator%(const integer& r1, _ulong r2) { return mpz_fdiv_ui(r1._x, r2); }
		friend _ulong operator%(_ulong r1, const integer& r2)
		{
			if (mpz_cmpabs_ui(r2._x, r1) > 0) return r1;
			return r1 % _ulong(abs(r2));
		}
		friend _ulong operator%(_ulong r1, integer&& r2)
		{
			mpz_abs(r2._x, r2._x);
			if (r2 > r1) return r1;
			auto ret = r1 % _ulong(r2);
			r2.reset();
			return ret;
		}

		integer operator+() const& { return *this; }
		integer operator+() && { return std::move(*this); }
		integer operator-() const&
		{
			integer temp;
			mpz_neg(temp._x, _x);
			return temp;
		}
		integer operator-() &&
		{
			mpz_neg(_x, _x);
			return std::move(*this);
		}

		friend std::optional<rational> parse(std::string_view str);
		friend std::optional<rational> _parse_mantisa(std::string_view str);
		friend mpz_srcptr _mp_ops<integer>::data(const integer& rop);
		friend integer _mp_ops<integer>::ceil(mpq_srcptr rop);
		friend integer _mp_ops<integer>::truncate(mpq_srcptr rop);
		friend integer _mp_ops<integer>::get(mpq_srcptr rop);
		friend integer _mp_ops<integer>::get(mpfr_srcptr rop, mpfr_rnd_t rnd);

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
		friend int iszero(const integer& x);
		friend int sgn(const integer& x);
		friend integer abs(const integer& r);
		friend integer abs(integer&& r);
		friend integer pow(const integer& op1, _ulong op2);
		friend integer pow(integer&& op1, _ulong op2);
		friend integer pow(_ulong op1, _ulong op2);
		friend integer factorial(_ulong n);
		friend int cmpabs(const integer& r1, const integer& r2);
	};
	inline int iszero(const integer& x) { return mpz_sgn(x._x) == 0; }
	inline int sgn(const integer& x) { return mpz_sgn(x._x); }
	inline integer abs(const integer& r)
	{
		integer temp;
		mpz_abs(temp._x, r._x);
		return temp;
	}
	inline integer abs(integer&& r)
	{
		mpz_abs(r._x, r._x);
		return std::move(r);
	}
	inline integer pow(const integer& op1, _ulong op2)
	{
		integer temp;
		mpz_pow_ui(temp._x, op1._x, op2);
		return temp;
	}
	inline integer pow(integer&& op1, _ulong op2)
	{
		mpz_pow_ui(op1._x, op1._x, op2);
		return std::move(op1);
	}
	inline integer pow(_ulong op1, _ulong op2)
	{
		integer temp;
		mpz_ui_pow_ui(temp._x, op1, op2);
		return temp;
	}
	// return n!
	inline integer factorial(_ulong n)
	{
		integer temp;
		mpz_fac_ui(temp._x, n);
		return temp;
	}
	inline int cmpabs(const integer& r1, const integer& r2) { return mpz_cmpabs(r1._x, r2._x); }
	inline _ulong _mp_ops<_ulong>::get(mpq_srcptr op) { return _ulong(_mp_ops<integer>::get(op)); }
	inline _long _mp_ops<_long>::get(mpq_srcptr op) { return _long(_mp_ops<integer>::get(op)); }
	inline integer _mp_ops<integer>::ceil(mpq_srcptr rop)
	{
		integer z;
		mpz_cdiv_q(z._x, mpq_numref(rop), mpq_denref(rop));
		return z;
	}
	inline integer _mp_ops<integer>::truncate(mpq_srcptr rop)
	{
		integer z;
		mpz_tdiv_q(z._x, mpq_numref(rop), mpq_denref(rop));
		return z;
	}
	inline integer _mp_ops<integer>::get(mpq_srcptr rop)
	{
		integer z;
		mpz_fdiv_q(z._x, mpq_numref(rop), mpq_denref(rop));
		return z;
	}
	inline integer _mp_ops<integer>::get(mpfr_srcptr rop, mpfr_rnd_t rnd)
	{
		integer z;
		mpfr_get_z(z._x, rop, rnd);
		return z;
	}
	inline mpz_srcptr _mp_ops<integer>::data(const integer& op) { return op._x; }
}  // namespace qboot::mp

#endif  // QBOOT_MP_INTEGER_HPP_
