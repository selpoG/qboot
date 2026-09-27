#include <array>       // for array
#include <atomic>      // for atomic
#include <chrono>      // for seconds
#include <cstdint>     // for uint32_t
#include <cstdlib>     // for getenv
#include <exception>   // for exception
#include <functional>  // for function
#include <future>      // for future_status
#include <iostream>    // for cerr
#include <map>         // for map
#include <stdexcept>   // for runtime_error, logic_error
#include <string>      // for string
#include <vector>      // for vector

#include "mpfr.h"  // for mpfr_free_cache

#include "qboot/qboot.hpp"  // for algebra types, task queues, numerical types

using qboot::algebra::RealFunction, qboot::algebra::ComplexFunction, qboot::algebra::FunctionSymmetry;
using qboot::algebra::Vector, qboot::algebra::Matrix, qboot::algebra::Polynomial;
using qboot::mp::real, qboot::mp::rational;

namespace
{
	void require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	template <class T, class F>
	void require_throws(F f)
	{
		try
		{
			f();
		}
		catch (const T&)
		{
			return;
		}
		throw std::runtime_error("Expected exception was not thrown");
	}

	void scalar()
	{
		const Polynomial p{real(2), real(3), real(-5)};
		for (const auto& c : {real(-4), real(0), real(1), real("0.5")})
		{
			auto q = mul_scalar(c, p);
			for (const auto& x : {real(-2), real(0), real(3)})
				require(q.eval(x) == c * p.eval(x), "polynomial scalar multiplication");
			auto r = p.clone();
			r *= c;
			require(q == r, "scalar multiplication disagrees with in-place multiplication");
		}
		Matrix<Polynomial> m(1, 1);
		m.at(0, 0) = p.clone();
		require(mul_scalar(real(4), m).at(0, 0).eval(real(2)) == 4 * p.eval(real(2)), "matrix scalar multiplication");
	}

	void subtraction()
	{
		const std::array sym{FunctionSymmetry::Even, FunctionSymmetry::Odd, FunctionSymmetry::Mixed};
		for (auto sx : sym)
			for (auto sy : sym)
			{
				ComplexFunction<real> x(3, sx), y(3, sy);
				for (uint32_t dy = 0; dy <= 1; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= 3; ++dx)
					{
						if (qboot::algebra::_matches(sx, dx)) x.at(dx, dy) = real(2 + dx + dy);
						if (qboot::algebra::_matches(sy, dx)) y.at(dx, dy) = real(7 + dx + dy);
					}
				auto z = x - y;
				for (uint32_t dy = 0; dy <= 1; ++dy)
					for (uint32_t dx = 0; dx + 2 * dy <= 3; ++dx)
					{
						if (!qboot::algebra::_matches(z.symmetry(), dx)) continue;
						auto a = qboot::algebra::_matches(sx, dx) ? x.at(dx, dy) : real(0);
						auto b = qboot::algebra::_matches(sy, dx) ? y.at(dx, dy) : real(0);
						require(z.at(dx, dy) == a - b, "subtraction across function symmetries");
					}
			}
	}

	void norm()
	{
		require(Polynomial().norm() == 0, "zero polynomial norm");
		require(Vector<real>().abs() == 0, "empty vector norm");
		require(Matrix<real>(0, 3).norm() == 0, "empty matrix norm");
		require(Vector<Vector<real>>().norm() == 0, "empty nested vector norm");
		require(Polynomial{real(3), real(4)}.norm() == 25, "nonzero polynomial norm");
	}

	void linear()
	{
		Polynomial p;
		p._mul_linear(real(2));
		require(p.iszero(), "linear multiplication of zero polynomial");
		Polynomial q{real(2), real(3)};
		q._mul_linear(real(4));
		require(q == Polynomial{real(8), real(14), real(3)}, "linear multiplication of nonzero polynomial");
	}

	void shift()
	{
		for (uint32_t lambda = 0; lambda <= 3; ++lambda)
			for (uint32_t p = 0; p <= 5; ++p)
			{
				RealFunction<real> f(lambda);
				for (uint32_t i = 0; i <= lambda; ++i) f.at(i) = real(i + 1);
				f.shift(p);
				require(f.lambda() == lambda, "shift changed truncation order");
				for (uint32_t i = 0; i <= lambda; ++i)
					require(f.at(i) == (i < p ? real(0) : real(i - p + 1)), "truncated power series shift");
			}
		qboot::Context context(4, 1, rational(3));
		require(context.rho_to_z()._total_memory() == 4, "low-order context construction");
	}

	void parsing()
	{
		for (const auto* s : {"1/0", "0/0", "-7/00"}) require(!qboot::mp::parse(s), "zero denominator accepted");
		require(qboot::mp::parse("2/4").value() == rational("1/2"), "rational canonicalization");
		require(qboot::mp::parse("0/5").value() == 0, "zero rational parsing");
	}

	void parallel()
	{
		for (uint32_t p : {0u, 1u, 2u, 4u})
		{
			std::vector<std::function<int()>> ints{[]() -> int { throw std::runtime_error("task failed"); }};
			std::vector<std::function<bool()>> bools{[]() -> bool { throw std::runtime_error("task failed"); }};
			std::vector<std::function<void()>> voids{[] { throw std::runtime_error("task failed"); }};
			require_throws<std::runtime_error>([&] { qboot::_parallel_evaluate(ints, p); });
			require_throws<std::runtime_error>([&] { qboot::_parallel_evaluate(bools, p); });
			require_throws<std::runtime_error>([&] { qboot::_parallel_evaluate(voids, p); });
			ints = {[] { return 42; }};
			require(qboot::_parallel_evaluate(ints, p).at(0) == 42, "parallel evaluation after failure");
		}
	}

	void queue()
	{
		for (uint32_t p : {0u, 1u, 2u})
		{
			qboot::_task_queue q(p);
			auto value = q.push([] { return 42; });
			require(value.wait_for(std::chrono::seconds(2)) == std::future_status::ready,
			        "task queue made no progress");
			require(value.get() == 42, "queued task result");
			auto failure = q.push([]() -> int { throw std::runtime_error("task failed"); });
			require_throws<std::runtime_error>([&] { failure.get(); });
			auto void_failure = q.push([] { throw std::runtime_error("task failed"); });
			require_throws<std::runtime_error>([&] { void_failure.get(); });
			q.push([] {}).get();
			require(q.push([] { return 43; }).get() == 43, "task queue stopped after a failed task");
			q.signal_done();
			require_throws<std::logic_error>([&] { q.push([] {}); });
		}
		std::atomic<uint32_t> calls{0};
		qboot::_memoized<int(int)> memo(
		    [&calls](int) -> int
		    {
			    ++calls;
			    throw std::runtime_error("memoized failure");
		    });
		for (uint32_t i = 0; i < 2; ++i) require_throws<std::runtime_error>([&] { memo(0); });
		require(calls == 1, "memoized failure was not cached");
	}
}  // namespace

int main()
{
	qboot::mp::global_prec = 256;
	qboot::mp::global_rnd = MPFR_RNDN;
	try
	{
		const std::map<std::string, std::function<void()>> tests{
		    {"scalar", scalar}, {"subtraction", subtraction}, {"norm", norm},         {"linear", linear},
		    {"shift", shift},   {"parsing", parsing},         {"parallel", parallel}, {"queue", queue}};
		if (const auto* name = std::getenv("QBOOT_REGRESSION"))
			tests.at(name)();
		else
			for (const auto& entry : tests) entry.second();
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
