#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include "qboot/algebra/complex_function.hpp"
#include "qboot/algebra/matrix.hpp"
#include "qboot/algebra/polynomial.hpp"
#include "qboot/algebra/real_function.hpp"
#include "qboot/block.hpp"
#include "qboot/bootstrap_equation.hpp"
#include "qboot/conformal_scale.hpp"
#include "qboot/context.hpp"
#include "qboot/hor_formula.hpp"
#include "qboot/hor_recursion.hpp"
#include "qboot/mp/real.hpp"
#include "qboot/my_filesystem.hpp"
#include "qboot/polynomial_program.hpp"
#include "qboot/primary_op.hpp"
#include "qboot/sdpb_input.hpp"
#include "qboot/task_queue.hpp"

using qboot::gBlock, qboot::Context, qboot::algebra::ComplexFunction, qboot::algebra::RealFunction,
    qboot::ConformalScale, qboot::PolynomialProgram, qboot::Sector;
using qboot::algebra::Vector, qboot::algebra::Matrix, qboot::algebra::Polynomial;
using qboot::mp::real, qboot::mp::rational, qboot::mp::parse;
using std::array, std::unique_ptr, std::cout, std::endl, std::map, std::optional, std::make_unique, std::move,
    std::pair, std::vector, std::string;
namespace fs = qboot::fs;

template <class T>
using dict = map<string, T, std::less<>>;
using Op = qboot::PrimaryOperator;
using Eq = qboot::Equation;
constexpr auto ContinuousType = qboot::SectorType::Continuous;
constexpr auto Odd = qboot::algebra::FunctionSymmetry::Odd;
constexpr auto Even = qboot::algebra::FunctionSymmetry::Even;

class WatchScope : public qboot::_event_base
{
#ifndef NDEBUG
	vector<string> tags_{};
	dict<std::chrono::system_clock::time_point> start_{};
	dict<std::chrono::system_clock::time_point> end_{};
	std::mutex mutex_{};
	static constexpr double to_seconds(std::chrono::system_clock::duration dur)
	{
		return double(std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count()) / 1e9;
	}

public:
	WatchScope() = default;
	WatchScope(const WatchScope&) = delete;
	WatchScope(WatchScope&&) noexcept = delete;
	WatchScope& operator=(const WatchScope&) = delete;
	WatchScope& operator=(WatchScope&&) noexcept = delete;
	~WatchScope() override;
	void on_begin(std::string_view tag) override
	{
		std::lock_guard<std::mutex> guard(mutex_);
		tags_.emplace_back(tag);
		start_.emplace(tag, std::chrono::system_clock::now());
	}
	void on_end(std::string_view tag) override
	{
		std::lock_guard<std::mutex> guard(mutex_);
		end_.emplace(tag, std::chrono::system_clock::now());
		auto ss = std::string(tag);
		std::cout << std::setw(25) << std::left << ("[" + ss + "] ") << to_seconds(end_[ss] - start_[ss]) << std::endl;
	}
#endif
};
#ifndef NDEBUG
WatchScope::~WatchScope() = default;
#endif

[[maybe_unused]] static string name_single(const dict<rational>& deltas);
[[maybe_unused]] static string name_mixed(const dict<rational>& deltas);
[[maybe_unused]] static void single_ising(const Context& c, const dict<rational>& deltas, uint32_t numax,
                                          uint32_t maxspin, qboot::_event_base* event);
[[maybe_unused]] static void mixed_ising(const Context& c, const dict<rational>& deltas, uint32_t numax,
                                         uint32_t maxspin, qboot::_event_base* event);

string name_single(const dict<rational>& deltas)
{
	return string("single-ising-") + deltas.at("s").str('#') + "-" + deltas.at("e").str('#');
}
string name_mixed(const dict<rational>& deltas)
{
	return string("mixed-ising-") + deltas.at("s").str('#') + "-" + deltas.at("e").str('#') + "-" +
	       deltas.at("e1").str('#');
}

void mixed_ising(const Context& c, const dict<rational>& deltas, uint32_t numax = 20, uint32_t maxspin = 24,
                 qboot::_event_base* event = nullptr)
{
	dict<Op> ops;
	ops.emplace("s", Op(real(deltas.at("s")), 0, c));
	ops.emplace("e", Op(real(deltas.at("e")), 0, c));
	ops.emplace("e1", Op(real(deltas.at("e1")), 0, c));
	auto ext = [&ops](auto o1, auto o2, auto o3, auto o4) {
		return array{ops.at(o1), ops.at(o2), ops.at(o3), ops.at(o4)};
	};
	Op one{c}, T{2, c};
	vector<Sector> secs{
	    {"unit", 1, {real(1)}}, {"scalar", 2}, {"e1", 2}, {"T", 2, {ops.at("s").delta(), ops.at("e").delta()}}};
	{
		Sector even("even", 2, ContinuousType);
		// even.add_op(0, real("3.815"), real("3.845"));
		even.add_op(0, real(6));
		even.add_op(2, real(5));
		for (uint32_t spin = 4; spin <= maxspin; spin += 2) even.add_op(spin);
		secs.push_back(even);
		Sector oddp("odd+", 1, ContinuousType);
		oddp.add_op(0, real(3));
		// oddp.add_op(0, real("5.25"), real("5.33"));
		for (uint32_t spin = 2; spin <= maxspin; spin += 2) oddp.add_op(spin);
		secs.push_back(oddp);
		Sector oddm("odd-", 1, ContinuousType);
		for (uint32_t spin = 1; spin <= maxspin; spin += 2) oddm.add_op(spin);
		secs.push_back(oddm);
	}
	qboot::BootstrapEquation boot(c, secs, numax);
	{
		Eq eq(boot, Odd);
		eq.add("unit", one, ext("s", "s", "s", "s"));
		eq.add("scalar", 0, 0, ops.at("e"), ext("s", "s", "s", "s"));
		eq.add("e1", 0, 0, ops.at("e1"), ext("s", "s", "s", "s"));
		eq.add("T", 0, 0, T, ext("s", "s", "s", "s"));
		eq.add("even", 0, 0, ext("s", "s", "s", "s"));
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("unit", one, ext("e", "e", "e", "e"));
		eq.add("scalar", 1, 1, ops.at("e"), ext("e", "e", "e", "e"));
		eq.add("e1", 1, 1, ops.at("e1"), ext("e", "e", "e", "e"));
		eq.add("T", 1, 1, T, ext("e", "e", "e", "e"));
		eq.add("even", 1, 1, ext("e", "e", "e", "e"));
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("scalar", 0, 0, ops.at("s"), ext("e", "s", "e", "s"));
		eq.add("odd+", ext("e", "s", "e", "s"));
		eq.add("odd-", ext("e", "s", "e", "s"));
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("unit", one, ext("e", "e", "s", "s"));
		eq.add("scalar", 0, 0, ops.at("s"), ext("e", "s", "s", "e"));
		eq.add("scalar", 0, 1, ops.at("e"), ext("e", "e", "s", "s"));
		eq.add("e1", 0, 1, ops.at("e1"), ext("e", "e", "s", "s"));
		eq.add("T", 0, 1, T, ext("e", "e", "s", "s"));
		eq.add("odd+", ext("e", "s", "s", "e"));
		eq.add("odd-", real(-1), ext("e", "s", "s", "e"));
		eq.add("even", 0, 1, ext("e", "e", "s", "s"));
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Even);
		eq.add("unit", one, ext("e", "e", "s", "s"));
		eq.add("scalar", 0, 0, real(-1), ops.at("s"), ext("e", "s", "s", "e"));
		eq.add("scalar", 0, 1, ops.at("e"), ext("e", "e", "s", "s"));
		eq.add("e1", 0, 1, ops.at("e1"), ext("e", "e", "s", "s"));
		eq.add("T", 0, 1, T, ext("e", "e", "s", "s"));
		eq.add("odd+", real(-1), ext("e", "s", "s", "e"));
		eq.add("odd-", ext("e", "s", "s", "e"));
		eq.add("even", 0, 1, ext("e", "e", "s", "s"));
		boot.add_equation(eq);
	}
	boot.finish();
	auto root = fs::current_path() / name_mixed(deltas);
	auto pmp = boot.ope_maximize("T", "unit", c.parallel(), event);
	// auto pmp = boot.find_contradiction("unit", 8, event);
	move(pmp).create_input(c.parallel(), event).write(root, c.parallel());
}

void single_ising(const Context& c, const dict<rational>& deltas, uint32_t numax = 20, uint32_t maxspin = 24,
                  qboot::_event_base* event = nullptr)
{
	dict<Op> ops;
	ops.emplace("s", Op(real(deltas.at("s")), 0, c));
	ops.emplace("e", Op(real(deltas.at("e")), 0, c));
	auto ext = [&ops](auto o1, auto o2, auto o3, auto o4) {
		return array{ops.at(o1), ops.at(o2), ops.at(o3), ops.at(o4)};
	};
	Op u{c}, T{2, c};
	vector<Sector> secs{{"unit", 1, {real(1)}}, {"e", 1}, {"T", 1}};
	{
		Sector even("even", 1, ContinuousType);
		even.add_op(0, real("1.409"), real("1.4135"));
		even.add_op(0, real(6));
		even.add_op(2, real(5));
		for (uint32_t spin = 4; spin <= maxspin; spin += 2) even.add_op(spin);
		secs.push_back(move(even));
	}
	qboot::BootstrapEquation boot(c, move(secs), numax);
	Eq eq(boot, Odd);
	eq.add("unit", u, ext("s", "s", "s", "s"));
	eq.add("e", ops.at("e"), ext("s", "s", "s", "s"));
	eq.add("T", T, ext("s", "s", "s", "s"));
	eq.add("even", ext("s", "s", "s", "s"));
	boot.add_equation(move(eq));
	boot.finish();
	auto root = fs::current_path() / name_single(deltas);
	auto pmp = boot.find_contradiction("unit", c.parallel(), event);
	move(pmp).create_input(c.parallel(), event).write(root, c.parallel());
}

int main(int argc, char* argv[])
{
	qboot::mp::global_prec = 1000;
	qboot::mp::global_rnd = MPFR_RNDN;
	constexpr uint32_t n_Max = 400, lambda = 14, dim = 3, maxspin = 24, parallel = 8;
	[[maybe_unused]] constexpr uint32_t numax = 6;
	if (argc != 1 && argc != 4)
	{
		cout << "usage: ./program [delta_s delta_e delta_e1]" << endl;
		return 1;
	}
	dict<rational> deltas;
	deltas["s"] = parse("0.5181489").value();
	deltas["e"] = parse("1.412625").value();
	deltas["e1"] = parse("3.83").value();
	if (argc == 4)
	{
		unique_ptr<char*[]> args(argv);
		deltas["s"] = parse(args[1]).value();
		deltas["e"] = parse(args[2]).value();
		deltas["e1"] = parse(args[3]).value();
		args.release();
	}
	Context c(n_Max, lambda, dim, parallel);
#ifndef NDEBUG
	auto stopwatch = std::make_unique<WatchScope>();
#else
	unique_ptr<qboot::_event_base> stopwatch{};
#endif
	mixed_ising(c, deltas, numax, maxspin, stopwatch.get());
	return 0;
}
