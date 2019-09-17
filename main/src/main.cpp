#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include "block.hpp"
#include "bootstrap_equation.hpp"
#include "complex_function.hpp"
#include "conformal_scale.hpp"
#include "context.hpp"
#include "hor_formula.hpp"
#include "hor_recursion.hpp"
#include "matrix.hpp"
#include "my_filesystem.hpp"
#include "polynomial.hpp"
#include "polynomial_program.hpp"
#include "primary_op.hpp"
#include "real.hpp"
#include "real_function.hpp"
#include "sdpb_input.hpp"

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using qboot::gBlock, qboot::Context, algebra::ComplexFunction, algebra::RealFunction, qboot::ConformalScale,
    qboot::PolynomialProgram, qboot::Sector;
using std::array, std::unique_ptr, std::cout, std::endl, std::map, std::optional, std::make_unique, std::move,
    std::pair, std::vector;
namespace fs = qboot::fs;

using R = mpfr::real;
using Op = qboot::PrimaryOperator;
using Eq = qboot::Equation;
using PolIneq = qboot::PolynomialInequalityWithCoeffs;
constexpr auto ContinuousType = qboot::SectorType::Continuous;
constexpr auto Odd = algebra::FunctionSymmetry::Odd;
constexpr auto Even = algebra::FunctionSymmetry::Even;

[[maybe_unused]] static void single_ising(const Context& c, const Op& s, const Op& e, uint32_t numax, uint32_t maxspin);
[[maybe_unused]] static void mixed_ising(const Context& c, const Op& s, const Op& e, const Op& e1, uint32_t numax,
                                         uint32_t maxspin);
[[maybe_unused]] static void test_sdpb();

void mixed_ising(const Context& c, const Op& s, const Op& e, const Op& e1, uint32_t numax = 20, uint32_t maxspin = 24)
{
	Op u{c}, T{2, c};
	vector<Sector> secs{{"unit", 1, {R(1)}}, {"scalar", 2}, {"e1", 2}, {"T", 2, {s.delta(), e.delta()}}};
	{
		Sector even("even", 2, ContinuousType);
		// even.add_op(0, R("3.815"), R("3.845"));
		even.add_op(0, R(6));
		even.add_op(2, R(5));
		for (uint32_t spin = 4; spin <= maxspin; spin += 2) even.add_op(spin);
		secs.push_back(even);
		Sector oddp("odd+", 1, ContinuousType);
		oddp.add_op(0, R(3));
		// oddp.add_op(0, R("5.25"), R("5.33"));
		for (uint32_t spin = 2; spin <= maxspin; spin += 2) oddp.add_op(spin);
		secs.push_back(oddp);
		Sector oddm("odd-", 1, ContinuousType);
		for (uint32_t spin = 1; spin <= maxspin; spin += 2) oddm.add_op(spin);
		secs.push_back(oddm);
	}
	qboot::BootstrapEquation boot(c, secs, numax);
	{
		Eq eq(boot, Odd);
		eq.add("unit", u, {s, s, s, s});
		eq.add("scalar", 0, 0, e, {s, s, s, s});
		eq.add("e1", 0, 0, e1, {s, s, s, s});
		eq.add("T", 0, 0, T, {s, s, s, s});
		eq.add("even", 0, 0, {s, s, s, s});
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("unit", u, {e, e, e, e});
		eq.add("scalar", 1, 1, e, {e, e, e, e});
		eq.add("e1", 1, 1, e1, {e, e, e, e});
		eq.add("T", 1, 1, T, {e, e, e, e});
		eq.add("even", 1, 1, {e, e, e, e});
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("scalar", 0, 0, s, {e, s, e, s});
		eq.add("odd+", {e, s, e, s});
		eq.add("odd-", {e, s, e, s});
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Odd);
		eq.add("unit", u, {e, e, s, s});
		eq.add("scalar", 0, 0, s, {e, s, s, e});
		eq.add("scalar", 0, 1, e, {e, e, s, s});
		eq.add("e1", 0, 1, e1, {e, e, s, s});
		eq.add("T", 0, 1, T, {e, e, s, s});
		eq.add("odd+", {e, s, s, e});
		eq.add("odd-", R(-1), {e, s, s, e});
		eq.add("even", 0, 1, {e, e, s, s});
		boot.add_equation(eq);
	}
	{
		Eq eq(boot, Even);
		eq.add("unit", u, {e, e, s, s});
		eq.add("scalar", 0, 0, R(-1), s, {e, s, s, e});
		eq.add("scalar", 0, 1, e, {e, e, s, s});
		eq.add("e1", 0, 1, e1, {e, e, s, s});
		eq.add("T", 0, 1, T, {e, e, s, s});
		eq.add("odd+", R(-1), {e, s, s, e});
		eq.add("odd-", {e, s, s, e});
		eq.add("even", 0, 1, {e, e, s, s});
		boot.add_equation(eq);
	}
	boot.finish();
	auto root =
	    fs::current_path() / ("mixed-ising-" + s.delta().str(8) + "-" + e.delta().str(8) + "-" + e1.delta().str(8));
	auto pmp = boot.ope_maximize("T", "unit", true);
	// auto pmp = boot.find_contradiction("unit", true);
	move(pmp).create_input().write(root);
}

void single_ising(const Context& c, const Op& s, const Op& e, uint32_t numax = 20, uint32_t maxspin = 24)
{
	Op u{c}, T{2, c};
	vector<Sector> secs{{"unit", 1, {R(1)}}, {"e", 1}, {"T", 1}};
	{
		Sector even("even", 1, ContinuousType);
		even.add_op(0, R("1.409"), R("1.4135"));
		even.add_op(0, R(6));
		even.add_op(2, R(5));
		for (uint32_t spin = 4; spin <= maxspin; spin += 2) even.add_op(spin);
		secs.push_back(move(even));
	}
	qboot::BootstrapEquation boot(c, move(secs), numax);
	Eq eq(boot, Odd);
	eq.add("unit", u, {s, s, s, s});
	eq.add("e", e, {s, s, s, s});
	eq.add("T", T, {s, s, s, s});
	eq.add("even", {s, s, s, s});
	boot.add_equation(move(eq));
	boot.finish();
	auto root = fs::current_path() / ("single-ising-" + s.delta().str(8) + "-" + e.delta().str(8));
	auto pmp = boot.find_contradiction("unit");
	move(pmp).create_input().write(root);
}

class TestScale : public qboot::ScaleFactor
{
	static constexpr uint32_t deg_ = 4;
	Vector<R> xs_;

public:
	TestScale() : xs_(qboot::sample_points(deg_)) {}
	TestScale(const TestScale&) = delete;
	TestScale& operator=(const TestScale&) = delete;
	TestScale(TestScale&&) noexcept = default;
	TestScale& operator=(TestScale&&) noexcept = default;
	~TestScale() override;
	[[nodiscard]] uint32_t max_degree() const override { return deg_; }
	[[nodiscard]] R eval(const R& v) const override { return mpfr::exp(-v); }
	[[nodiscard]] Vector<R> sample_scalings() override
	{
		Vector<R> sc(deg_ + 1);
		for (uint32_t i = 0; i <= deg_; ++i) sc[i] = eval(xs_[i]);
		return sc;
	}
	[[nodiscard]] R sample_point(uint32_t k) override { return xs_[k]; }
	[[nodiscard]] Vector<R> sample_points() override { return xs_.clone(); }
	[[nodiscard]] Polynomial bilinear_base(uint32_t m) override
	{
		assert(m <= deg_ / 2);
		switch (m)
		{
		case 0: return {R(1)};                  // 1
		case 1: return {R(-1), R(1)};           // -1 + x
		default: return {R(1), R(-2), R(0.5)};  // 1 - 2 x + x ^ 2 / 2
		}
	}
	[[nodiscard]] Vector<Polynomial> bilinear_bases() override
	{
		Vector<Polynomial> q(deg_ / 2 + 1);
		for (uint32_t i = 0; i <= deg_ / 2; ++i) q[i] = bilinear_base(i);
		return q;
	}
};
TestScale::~TestScale() = default;

void test_sdpb()
{
	// https://github.com/davidsd/sdpb/blob/master/test/test.xml
	// {1 + x ^ 4, x ^ 2 + x ^ 4 / 12}
	Vector<Polynomial> elm{{R(1), R(0), R(0), R(0), R(1)}, {R(0), R(0), R(1), R(0), R(1) / 12}};

	PolynomialProgram prg(1);
	prg.objective_constant() = R(0);
	prg.objectives({R(-1)});
	auto ineq = make_unique<PolIneq>(1u, make_unique<TestScale>(), Vector{elm[1].clone()}, -elm[0]);
	prg.add_inequality(move(ineq));

	auto sdpb = move(prg).create_input();
	auto root = fs::current_path() / "test";
	sdpb.write(root);
}

int main(int argc, char* argv[])
{
	mpfr::global_prec = 1000;
	mpfr::global_rnd = MPFR_RNDN;
	constexpr uint32_t n_Max = 400, lambda = 14, dim_ = 3, maxspin = 24;
	[[maybe_unused]] constexpr uint32_t numax = 6;
	assert(argc == 4);
	unique_ptr<char*[]> args(argv);
	auto d_s = R(args[1]);
	auto d_e = R(args[2]);
	auto d_e1 = R(args[3]);
	Context c(n_Max, lambda, dim_);
	Op sigma(d_s, 0, c);
	Op epsilon(d_e, 0, c);
	Op e1(d_e1, 0, c);
	mixed_ising(c, sigma, epsilon, e1, numax, maxspin);
	args.release();
	return 0;
}
