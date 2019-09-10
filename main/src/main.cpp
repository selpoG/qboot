#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "block.hpp"
#include "bootstrap_equation.hpp"
#include "complex_function.hpp"
#include "conformal_scale.hpp"
#include "context_variables.hpp"
#include "hor_formula.hpp"
#include "hor_recursion.hpp"
#include "matrix.hpp"
#include "polynomial.hpp"
#include "polynomial_program.hpp"
#include "primary_op.hpp"
#include "real.hpp"
#include "real_function.hpp"
#include "sdpb_input.hpp"

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using qboot::h_asymptotic, qboot::gBlock, qboot::Context, algebra::ComplexFunction, algebra::RealFunction,
    qboot::ConformalScale, qboot::PolynomialProgram;
using FunctionSymmetry = algebra::FunctionSymmetry;
using std::array, std::unique_ptr, std::cout, std::endl, std::map, std::optional, std::make_unique, std::move,
    std::pair;
namespace fs = std::filesystem;

using R = mpfr::real;
using Op = qboot::PrimaryOperator;
using GOp = qboot::GeneralPrimaryOperator;
using Block = qboot::ConformalBlock<Op>;
using GBlock = qboot::ConformalBlock<GOp>;
using GGBlock = qboot::GeneralConformalBlock;
using PolIneq = qboot::PolynomialInequalityWithCoeffs;
using EvalIneq = qboot::PolynomialInequalityEvaluated;

[[maybe_unused]] static void single_ising(const Context& c, const Op& s, const Op& e, uint32_t numax, uint32_t maxspin);
[[maybe_unused]] static void mixed_ising(const Context& c, const Op& s, const Op& e, const Op& e1, uint32_t numax,
                                         uint32_t maxspin);
[[maybe_unused]] static void test_sdpb();

void mixed_ising(const Context& c, const Op& s, const Op& e, const Op& e1, uint32_t numax = 20, uint32_t maxspin = 24)
{
	constexpr FunctionSymmetry Od = FunctionSymmetry::Odd;
	constexpr FunctionSymmetry Ev = FunctionSymmetry::Even;
	using E = qboot::Entry;
	using B = Block;
	using G = GGBlock;
	qboot::BootstrapEquation boot(c, array{pair{"unit", 1u}, pair{"scalar", 2u}, pair{"e1", 2u}, pair{"T", 2u},
	                                       pair{"even", 2u}, pair{"odd+", 1u}, pair{"odd-", 1u}});
	auto num_poles = [numax]([[maybe_unused]] uint32_t spin) { return numax; };
	// spectrum
	Op u{c.epsilon()}, T{2, c.epsilon()};  // do not register point-like spectrums ("unit", "scalar", "T")
	boot.register_ope("T", {R(s.delta()), R(e.delta())});  // from Ward identity
	// boot.register_operator("even", {0, num_poles(0), c.epsilon(), R("3.815"), R("3.845")});
	boot.register_operator("even", {0, num_poles(0), c.epsilon(), R("6")});
	// boot.register_operator("odd+", {0, num_poles(0), c.epsilon(), R("5.25"), R("5.33")});
	boot.register_operator("odd+", {0, num_poles(0), c.epsilon(), R(3)});
	boot.register_operator("even", {2, num_poles(2), c.epsilon(), R("5")});
	for (uint32_t spin = 1; spin <= maxspin; ++spin)
		if (spin % 2 == 0)
		{
			if (spin != 2) boot.register_operator("even", {spin, num_poles(spin), c.epsilon()});
			boot.register_operator("odd+", {spin, num_poles(spin), c.epsilon()});
		}
		else
			boot.register_operator("odd-", {spin, num_poles(spin), c.epsilon()});
	// equations
	boot.add_equation(Od, array{pair{"unit", E(B{u, s, s, s, s, Od})}, pair{"scalar", E(0, 0, B{e, s, s, s, s, Od})},
	                            pair{"e1", E(0, 0, B{e1, s, s, s, s, Od})}, pair{"T", E(0, 0, B{T, s, s, s, s, Od})},
	                            pair{"even", E(0, 0, G{s, s, s, s, Od})}});
	boot.add_equation(Od, array{pair{"unit", E(B{u, e, e, e, e, Od})}, pair{"scalar", E(1, 1, B{e, e, e, e, e, Od})},
	                            pair{"e1", E(1, 1, B{e1, e, e, e, e, Od})}, pair{"T", E(1, 1, B{T, e, e, e, e, Od})},
	                            pair{"even", E(1, 1, G{e, e, e, e, Od})}});
	boot.add_equation(Od, array{pair{"scalar", E(1, 1, B{s, e, s, e, s, Od})}, pair{"odd+", E(G{e, s, e, s, Od})},
	                            pair{"odd-", E(G{e, s, e, s, Od})}});
	boot.add_equation(Od,
	                  array{pair{"unit", E(B{u, e, e, s, s, Od})}, pair{"scalar", E(0, 0, B{s, e, s, s, e, Od})},
	                        pair{"scalar", E(0, 1, B{e, e, e, s, s, Od})}, pair{"e1", E(0, 1, B{e1, e, e, s, s, Od})},
	                        pair{"T", E(0, 1, B{T, e, e, s, s, Od})}, pair{"odd+", E(G{e, s, s, e, Od})},
	                        pair{"odd-", E(R(-1), G{e, s, s, e, Od})}, pair{"even", E(0, 1, G{e, e, s, s, Od})}});
	boot.add_equation(Ev,
	                  array{pair{"unit", E(B{u, e, e, s, s, Ev})}, pair{"scalar", E(0, 0, R(-1), B{s, e, s, s, e, Ev})},
	                        pair{"scalar", E(0, 1, B{e, e, e, s, s, Ev})}, pair{"e1", E(0, 1, B{e1, e, e, s, s, Ev})},
	                        pair{"T", E(0, 1, B{T, e, e, s, s, Ev})}, pair{"odd+", E(R(-1), G{e, s, s, e, Ev})},
	                        pair{"odd-", E(G{e, s, s, e, Ev})}, pair{"even", E(0, 1, G{e, e, s, s, Ev})}});
	auto root =
	    fs::current_path() / ("mixed-ising-" + s.delta().str(8) + "-" + e.delta().str(8) + "-" + e1.delta().str(8));
	auto pmp = boot.ope_minimize("T", "unit", true);
	move(pmp).create_input().write_all(root);
}

void single_ising(const Context& c, const Op& s, const Op& e, uint32_t numax = 20, uint32_t maxspin = 24)
{
	constexpr FunctionSymmetry Od = FunctionSymmetry::Odd;
	using E = qboot::Entry;
	using B = Block;
	using G = GGBlock;
	qboot::BootstrapEquation boot(c, array{pair{"unit", 1u}, pair{"e", 1u}, pair{"T", 1u}, pair{"even", 1u}});
	auto num_poles = [numax](uint32_t spin) { return numax + std::min(numax, spin) / 2; };
	Op u{c.epsilon()}, T{2, c.epsilon()};
	boot.register_operator("even", {0, num_poles(0), c.epsilon(), R("1.409"), R("1.4135")});
	boot.register_operator("even", {0, num_poles(0), c.epsilon(), R(6)});
	boot.register_operator("even", {2, num_poles(2), c.epsilon(), R(5)});
	for (uint32_t spin = 4; spin <= maxspin; spin += 2)
		boot.register_operator("even", {spin, num_poles(spin), c.epsilon()});
	boot.add_equation(Od, array{pair{"unit", E(B{u, s, s, s, s, Od})}, pair{"e", E(B{e, s, s, s, s, Od})},
	                            pair{"T", E(B{T, s, s, s, s, Od})}, pair{"even", E(G{s, s, s, s, Od})}});
	auto root = fs::current_path() / ("single-ising-" + s.delta().str(8) + "-" + e.delta().str(8));
	auto pmp = boot.find_contradiction("unit");
	move(pmp).create_input().write_all(root);
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
	auto ineq = make_unique<PolIneq>(1u, std::make_unique<TestScale>(), Vector{elm[1].clone()}, -elm[0]);
	prg.add_inequality(move(ineq));

	auto sdpb = move(prg).create_input();
	auto root = fs::current_path() / "test";
	sdpb.write_all(root);
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
	Op sigma(d_s, 0, c.epsilon());
	Op epsilon(d_e, 0, c.epsilon());
	Op e1(d_e1, 0, c.epsilon());
	mixed_ising(c, sigma, epsilon, e1, numax, maxspin);
	args.release();
	return 0;
}
