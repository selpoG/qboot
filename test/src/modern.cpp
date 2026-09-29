#include <algorithm>  // for ranges::sort
#include <array>      // for array
#include <compare>    // for is_eq, is_gt, is_lt, partial_ordering, strong_ordering, three_way_comparable
#include <concepts>   // for constructible_from, same_as
#include <exception>  // for exception
#include <future>     // for future, future_error, future_errc, promise
#include <iostream>   // for cerr
#include <limits>     // for numeric_limits
#include <memory>     // for make_unique
#include <span>       // for span
#include <stdexcept>  // for runtime_error
#include <string>     // for string
#include <utility>    // for as_const, declval, forward, move

#include "mpfr.h"  // for mpfr_free_cache

#include "qboot/qboot.hpp"  // for algebra types, blocks, numerical types, task queue

using qboot::algebra::Vector, qboot::algebra::Matrix, qboot::algebra::Polynomial;
using qboot::mp::integer, qboot::mp::rational, qboot::mp::real;

namespace
{
	template <class T>
	concept HasDelta = requires(const T& block)
	{
		block.delta();
	};
	template <class T>
	concept CanFixDelta = requires(const T& block, const real& delta)
	{
		block.fix_delta(delta);
	};
	template <class T>
	concept HasInnerProduct = requires(const Matrix<T>& m, const Vector<T>& v)
	{
		m.inner_product(v);
	};
	template <class T>
	concept HasTemporaryView = requires(T value)
	{
		std::move(value).view();
	};
	template <class T>
	concept HasTemporaryRow = requires(T value)
	{
		std::move(value).row_view(0);
	};
	template <class T>
	concept IntegerComparable = requires(const integer& x, const T& y)
	{
		x <=> y;
	};

	static_assert(HasDelta<qboot::ConformalBlock<qboot::PrimaryOperator>>);
	static_assert(!HasDelta<qboot::ConformalBlock<qboot::GeneralPrimaryOperator>>);
	static_assert(!CanFixDelta<qboot::ConformalBlock<qboot::PrimaryOperator>>);
	static_assert(CanFixDelta<qboot::ConformalBlock<qboot::GeneralPrimaryOperator>>);
	static_assert(HasInnerProduct<real> && !HasInnerProduct<Polynomial>);
	static_assert(!HasTemporaryView<Vector<real>> && !HasTemporaryRow<Matrix<real>>);
	static_assert(std::same_as<decltype(std::declval<const Vector<real>&>().view()), std::span<const real>>);
	static_assert(qboot::algebra::_ring<Polynomial> && qboot::algebra::_ring<Matrix<real>>);
	static_assert(!qboot::algebra::_ring<int>);
	static_assert(std::three_way_comparable<integer, std::strong_ordering>);
	static_assert(std::three_way_comparable<rational, std::strong_ordering>);
	static_assert(std::three_way_comparable<real, std::partial_ordering>);
	static_assert(IntegerComparable<int> && IntegerComparable<double> && !IntegerComparable<std::string>);
	static_assert(!std::constructible_from<real, std::array<int, 2>>);

	template <class T, class S>
	concept ScalarOperations = requires(T& x, const T& c, const S& scalar)
	{
		{ x *= scalar } -> std::same_as<T&>;
		{ x /= scalar } -> std::same_as<T&>;
		{ mul_scalar(scalar, c) } -> std::same_as<T>;
		{ mul_scalar(scalar, std::move(x)) } -> std::same_as<T>;
		{ c / scalar } -> std::same_as<T>;
		{ std::move(x) / scalar } -> std::same_as<T>;
	};
	template <class T, class S>
	concept HasDot = requires(const T& x, const S& y) { dot(x, y); };
	template <class T>
	concept PolynomialArgument = requires(const Polynomial& p, const T& x) { p.eval(x); };
	template <class F>
	concept QueueArgument = requires(qboot::_task_queue& q, F&& f) { q.push(std::forward<F>(f)); };

	using Nested = Vector<Matrix<qboot::algebra::RealFunction<Polynomial>>>;
	static_assert(ScalarOperations<Polynomial, int> && ScalarOperations<Polynomial, rational>);
	static_assert(ScalarOperations<Vector<real>, int> && ScalarOperations<Matrix<Polynomial>, double>);
	static_assert(ScalarOperations<Nested, real> && ScalarOperations<Nested, rational>);
	static_assert(!ScalarOperations<Nested, std::string>);
	static_assert(!qboot::algebra::_scale_assignable<Nested, std::string>);
	static_assert(!qboot::algebra::_scalable<Nested, std::string>);
	static_assert(!qboot::algebra::_divisible<Nested, std::string>);
	static_assert(!qboot::algebra::_divide_assignable<Nested, std::string>);
	static_assert(!qboot::algebra::_scale_assignable<Polynomial, Polynomial>);
	static_assert(!qboot::algebra::_scalable<Vector<real>, Vector<real>>);
	static_assert(!qboot::algebra::_divisible<Vector<real>, Vector<real>>);
	static_assert(!qboot::algebra::_scalable<Matrix<Polynomial>, Matrix<Polynomial>>);
	static_assert(HasDot<Matrix<real>, Vector<Polynomial>> && HasDot<Vector<Polynomial>, Matrix<real>>);
	static_assert(HasDot<Matrix<Polynomial>, Matrix<Polynomial>>);
	static_assert(!HasDot<Matrix<Polynomial>, Vector<real>>);
	static_assert(!HasDot<Matrix<Vector<real>>, Matrix<Vector<real>>>);
	static_assert(PolynomialArgument<rational> && !PolynomialArgument<std::string>);
	static_assert(QueueArgument<int (*)()> && !QueueArgument<int>);

	template <class T>
	concept BlockArgument = requires { typename qboot::ConformalBlock<T>; };
	static_assert(BlockArgument<qboot::PrimaryOperator> && BlockArgument<qboot::GeneralPrimaryOperator>);
	static_assert(!BlockArgument<int> && !BlockArgument<real>);
	static_assert(qboot::algebra::Ring<integer> && qboot::algebra::Ring<rational>);
	static_assert(qboot::algebra::Ring<Polynomial> && qboot::algebra::Ring<Vector<real>>);
	static_assert(qboot::algebra::Ring<Matrix<Polynomial>> && !qboot::algebra::Ring<int>);
	static_assert(qboot::algebra::Algebra<integer> && qboot::algebra::Algebra<Polynomial>);
	static_assert(!qboot::algebra::Algebra<Vector<real>>);

	void require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	template <class T, class U>
	void ordered(const T& low, const U& high)
	{
		require(std::is_lt(low <=> high) && std::is_gt(high <=> low), "three-way ordering");
		require(low < high && low <= high && high > low && high >= low, "rewritten ordering");
		require(low != high && high != low && !(low == high) && !(high == low), "rewritten inequality");
	}

	template <class T, class U>
	void unordered(const T& x, const U& y)
	{
		require((x <=> y) == std::partial_ordering::unordered, "NaN must be unordered");
		require((y <=> x) == std::partial_ordering::unordered, "reversed NaN comparison");
		require(x != y && y != x && !(x == y) && !(y == x), "NaN equality");
		require(!(x < y) && !(x <= y) && !(x > y) && !(x >= y), "NaN relational operators");
		require(!(y < x) && !(y <= x) && !(y > x) && !(y >= x), "reversed NaN relational operators");
	}

	void comparisons()
	{
		ordered(integer(-1), integer(2));
		ordered(integer(-1), 2U);
		ordered(-1, integer(2));
		ordered(integer(1), rational(3, 2U));
		ordered(rational(3, 2U), real(2));
		ordered(real(-1), integer(2));
		ordered(rational(-1), 2U);
		ordered(-1, rational(2));
		ordered(real(-1), 2U);
		ordered(-1.5, real(2));
		ordered(integer(1), 1.5);
		ordered(integer("9007199254740993"), real("9007199254740994"));
		const auto nan = qboot::mp::nan();
		unordered(nan, real(0));
		unordered(nan, nan);
		unordered(nan, integer(0));
		unordered(nan, rational(0));
		unordered(nan, 0);
		unordered(nan, 0.0);
		unordered(real(0), std::numeric_limits<double>::quiet_NaN());
		unordered(integer(0), std::numeric_limits<double>::quiet_NaN());
		const auto inf = std::numeric_limits<double>::infinity();
		ordered(integer(0), inf);
		ordered(-inf, integer(0));
		ordered(real(0), real(inf));
		require(real(inf) == inf && real(-inf) == -inf, "infinite equality");
		require(real("-0") == real("0") && std::is_eq(real("-0") <=> 0), "signed zero equality");
		require(integer(2) == rational(2) && rational(2) == real(2) && real(2) == 2.0, "mixed equality");
		std::array<integer, 3> numbers{integer(4), integer(-2), integer(1)};
		std::ranges::sort(numbers);
		require(numbers[0] == -2 && numbers[1] == 1 && numbers[2] == 4, "ranges comparison compatibility");
		const qboot::PrimaryOperator a(real(2), 0, rational(1)), b(real(3), 0, rational(1));
		ordered(a, b);
	}

	void views()
	{
		Vector<real> empty;
		require(empty.view().empty() && empty.iszero(), "empty span");
		Vector<real> v{real(1), real(2), real(3)};
		auto middle = v.view().subspan(1, 1);
		middle[0] = 5;
		require(v[1] == 5 && std::as_const(v).view()[1] == 5, "span must refer to original coefficients");
		Matrix<real> m(2, 3);
		m.row_view(1)[2] = 7;
		require(m.at(1, 2) == 7 && std::as_const(m).row_view(1)[2] == 7, "matrix row span");
		require(m.row_view(0).size() == 3 && m.at(0, 2) == 0, "row boundaries");
		Matrix<real> zero_columns(2, 0);
		require(zero_columns.row_view(1).empty(), "empty matrix row");
	}

	void exact_norms()
	{
		const integer n("-123456789012345678901234567890");
		require(n.norm() == n * n, "exact integer norm");
		auto moved_n = n.clone();
		require(std::move(moved_n).norm() == n * n, "rvalue integer norm");
		const rational r(integer(-7), integer(3));
		require(r.norm() == rational(integer(49), integer(9)), "exact rational norm");
		auto moved_r = r.clone();
		require(std::move(moved_r).norm() == r * r, "rvalue rational norm");
		require(integer{}.norm() == 0 && rational{}.norm() == 0, "zero norms");
		const Vector<integer> v{integer(-3), integer(4)};
		require(v.norm() == 25 && v.eval(real(2)) == v, "integer coefficient vector");
		const Vector<rational> fractions{rational(integer(1), integer(3)), rational(integer(2), integer(3))};
		require(fractions.norm() == rational(integer(5), integer(9)), "rational coefficient vector");
		require(Vector<integer>{}.norm() == 0 && Vector<rational>{}.norm() == 0, "empty exact vectors");
	}

	void native_integer_widths()
	{
		using qboot::mp::_long, qboot::mp::_ulong;
		const auto low = std::numeric_limits<_long>::min(), high = std::numeric_limits<_long>::max();
		const auto max = std::numeric_limits<_ulong>::max();
		require(_long(integer(low)) == low && _long(integer(high)) == high, "signed native integer limits");
		require(_ulong(integer(max)) == max, "unsigned native integer limit");
		require(integer(-1) % max == max - 1, "native-width remainder");
		require(qboot::mp::parse("0.000125").value() == rational(1, 8000U), "decimal place count");
	}

	void scalar_templates()
	{
		const Polynomial p{real(1), real(-2), real(3)};
		const Vector<real> v{real(2), real(-3)};
		require(mul_scalar(2, v) == mul_scalar(real(2), v), "integral scalar on const vector");
		for (const auto& c : {rational(-2), rational(3, 2U)})
		{
			auto scaled = mul_scalar(c, p);
			auto inplace = p.clone();
			inplace *= c;
			require(scaled == inplace, "const and in-place polynomial scaling");
			require(scaled / c == p, "polynomial scalar division");
			require(scaled.eval(c) == real(c) * p.eval(c), "rational polynomial evaluation");
			Matrix<Polynomial> m(1, 1);
			m.at(0, 0) = p.clone();
			require(mul_scalar(c, m).at(0, 0) == scaled, "nested scalar multiplication");
		}
	}

	void tasks()
	{
		std::future<int> active, pending;
		std::promise<void> entered, resume;
		const auto resumed = resume.get_future();
		{
			qboot::_task_queue q(1);
			auto move_only = [p = std::make_unique<int>(23)] { return *p; };
			static_assert(QueueArgument<decltype(move_only)> && !QueueArgument<decltype(move_only)&>);
			require(q.push(std::move(move_only)).get() == 23, "forwarded move-only callable");
			auto copyable = [] { return 29; };
			require(q.push(copyable).get() == 29 && copyable() == 29, "lvalue callable remains usable");
			require(q.push([p = std::make_unique<int>(42)] { return *p; }).get() == 42, "move-only task");
			active = q.push(
			    [&entered, &resumed]
			    {
				    entered.set_value();
				    resumed.wait();
				    return 17;
			    });
			entered.get_future().wait();
			pending = q.push([] { return 99; });
			q.signal_done();
			resume.set_value();
		}
		require(active.get() == 17, "shutdown must join active tasks");
		try
		{
			pending.get();
			throw std::runtime_error("shutdown executed a pending task");
		}
		catch (const std::future_error& error)
		{
			require(error.code() == std::future_errc::broken_promise, "pending task cancellation");
		}
		for (int i = 0; i < 20; ++i) { qboot::_task_queue idle(2); }
	}
}  // namespace

int main()
{
	qboot::mp::global_prec = 256;
	try
	{
		exact_norms();
		native_integer_widths();
		scalar_templates();
		comparisons();
		views();
		tasks();
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		mpfr_free_cache();
		return 1;
	}
	mpfr_free_cache();
	return 0;
}
