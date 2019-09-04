#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "block.hpp"
#include "bootstrap_equation.hpp"
#include "chol_and_inverse.hpp"
#include "complex_function.hpp"
#include "context_variables.hpp"
#include "damped_rational.hpp"
#include "hor_formula.hpp"
#include "hor_recursion.hpp"
#include "integral_decomp.hpp"
#include "matrix.hpp"
#include "partial_fraction.hpp"
#include "pole_data.hpp"
#include "polynomial.hpp"
#include "polynomial_program.hpp"
#include "primary_op.hpp"
#include "real.hpp"
#include "real_function.hpp"
#include "sdpb_input.hpp"

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using qboot::h_asymptotic, qboot::gBlock, qboot::Context, algebra::ComplexFunction, algebra::RealFunction,
    qboot::ConformalScale, qboot::PolynomialProgramming;
using FunctionSymmetry = algebra::FunctionSymmetry;
using std::array, std::unique_ptr, std::cout, std::endl, std::map, std::optional, std::make_unique, std::move,
    std::pair;
namespace fs = std::filesystem;

using R = mpfr::real<1000, MPFR_RNDN>;
using Op = qboot::PrimaryOperator<R>;
using GOp = qboot::GeneralPrimaryOperator<R>;
using Block = qboot::ConformalBlock<R, Op>;
using GBlock = qboot::ConformalBlock<R, GOp>;
using GGBlock = qboot::GeneralConformalBlock<R>;
using PolIneq = qboot::PolynomialInequalityWithCoeffs<R>;
using EvalIneq = qboot::PolynomialInequalityEvaluated<R>;

[[maybe_unused]] static Vector<R> to_vec(const unique_ptr<mpfr_t[]>& a, uint32_t s);
[[maybe_unused]] static Matrix<R> to_mat(const unique_ptr<mpfr_t[]>& a, uint32_t r, uint32_t c);
[[maybe_unused]] static ComplexFunction<R> to_hol(const unique_ptr<mpfr_t[]>& a, uint32_t l);
static RealFunction<R> to_func(const unique_ptr<mpfr_t[]>& a, uint32_t l);
static void test_rec(uint32_t nMax, const Op& op, const R& d12, const R& d34, const R& very_small);
static void test_real(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
                      const R& very_small);
static void test_g(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
                   const R& very_small);
static void test_op(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
                    const R& very_small);
static void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34,
                   const R& very_small);
[[maybe_unused]] static void single_ising(const Context<R>& c, const Op& s, const Op& e, uint32_t numax,
                                          uint32_t maxspin);
[[maybe_unused]] static void mixed_ising(const Context<R>& c, const Op& s, const Op& e, uint32_t numax,
                                         uint32_t maxspin);
[[maybe_unused]] static void test_sdpb();

template <class T, class CallBack_T>
void check(const T& p, const T& q, const R& very_small, CallBack_T on_err)
{
	auto err = (p - q).norm();
	if (err > very_small)
	{
		on_err();
		cout << "p = " << p << endl;
		cout << "q = " << q << endl;
		cout << "p - q = " << (p - q) << endl;
		cout << "err = " << err << endl;
		assert(false);
	}
}

Vector<R> to_vec(const unique_ptr<mpfr_t[]>& a, uint32_t s)
{
	Vector<R> v(s);
	for (uint32_t i = 0; i < s; ++i) v[i] = R(a[i]);
	return v;
}

Matrix<R> to_mat(const unique_ptr<mpfr_t[]>& a, uint32_t r, uint32_t c)
{
	Matrix<R> v(r, c);
	for (uint32_t i = 0; i < r; ++i)
		for (uint32_t j = 0; j < c; ++j) v.at(i, j) = R(a[i * c + j]);
	return v;
}

ComplexFunction<R> to_hol(const unique_ptr<mpfr_t[]>& a, uint32_t l)
{
	ComplexFunction<R> f(l);
	uint32_t i = 0;
	for (uint32_t dy = 0; dy <= l / 2; ++dy)
		for (uint32_t dx = 0; dx + 2 * dy <= l; ++dx) f.at(dx, dy) = R(a[i++]);
	return f;
}

RealFunction<R> to_func(const unique_ptr<mpfr_t[]>& a, uint32_t l)
{
	RealFunction<R> f(l);
	for (uint32_t k = 0; k <= l; ++k) f.at(k) = R(a[k]);
	return f;
}

void test_rec(uint32_t nMax, const Op& op, const R& d12, const R& d34, const R& very_small)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = qboot::hBlock_shifted<R>(op, S, P, nMax);
	auto q_ptr =
	    op.spin() == 0
	        ? qboot2::recursionSpinZeroVector(nMax, op.epsilon()._x, &delta._x, S._x, P._x, R::prec, R::rnd)
	        : qboot2::recursionNonZeroVector(nMax, op.epsilon()._x, ell._x, delta._x, S._x, P._x, R::prec, R::rnd);
	auto q = to_func(q_ptr, p.lambda());
	check(p, q, very_small, [&]() {
		cout << "test_rec(nMax=" << nMax << "op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")" << endl;
	});
}

void test_real(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
               const R& very_small)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = qboot::gBlock_real<R>(op, S, P, c);
	auto q_ptr = qboot2::real_axis_result(op.epsilon()._x, ell._x, &delta._x, S._x, P._x, cb);
	auto q = to_func(q_ptr, p.lambda());
	check(p, q, very_small, [&]() {
		cout << "test_real(c=" << c.str() << ", cb, op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")"
		     << endl;
	});
}

void test_g(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
            const R& very_small)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = gBlock<R>(op, S, P, c);
	auto q_ptr = qboot2::gBlock_full(op.epsilon()._x, ell._x, &delta._x, S._x, P._x, cb);
	auto q = to_hol(q_ptr, cb.lambda);
	check(p, q, very_small, [&]() {
		cout << "test_g(c=" << c.str() << ", cb, op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")" << endl;
	});
}

void test_op(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34,
             const R& very_small)
{
	test_rec(c.n_Max(), op, d12, d34, very_small);
	test_real(c, cb, op, d12, d34, very_small);
	test_g(c, cb, op, d12, d34, very_small);
}

void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34, const R& very_small)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = h_asymptotic(S, cb1);
	auto q_ptr = qboot2::h_asymptotic(cb1.epsilon()._x, S._x, cb2);
	auto q = to_func(q_ptr, p.lambda());
	check(p, q, very_small,
	      [&]() { cout << "test_h(cb1=" << cb1.str() << ", cb2, d12=" << d12 << ", d34=" << d34 << ")" << endl; });
}

void mixed_ising(const Context<R>& c, const Op& s, const Op& e, uint32_t numax = 20, uint32_t maxspin = 24)
{
	constexpr FunctionSymmetry Od = FunctionSymmetry::Odd;
	constexpr FunctionSymmetry Ev = FunctionSymmetry::Even;
	using E = qboot::Entry<R>;
	using B = Block;
	using G = GGBlock;
	qboot::BootstrapEquation<R> boot(c, numax,
	                                 array{pair{"unit", 1u}, pair{"scalar", 2u}, pair{"T", 2u}, pair{"even", 2u},
	                                       pair{"odd+", 1u}, pair{"odd-", 1u}});
	// spectrum
	Op u{c.epsilon()}, T{2, c.epsilon()};  // do not register point-like spectrums ("unit", "scalar", "T")
	boot.register_operator("even", {0, c.epsilon(), R("3.8"), R("3.88")});
	boot.register_operator("even", {0, c.epsilon(), R(6)});
	boot.register_operator("odd+", {0, c.epsilon(), R("5.25"), R("5.33")});
	boot.register_operator("odd+", {0, c.epsilon(), R(8)});
	boot.register_operator("even", {2, c.epsilon(), R(5)});
	for (uint32_t spin = 1; spin <= maxspin; ++spin)
		if (spin % 2 == 0)
		{
			if (spin != 2) boot.register_operator("even", {spin, c.epsilon()});
			boot.register_operator("odd+", {spin, c.epsilon()});
		}
		else
			boot.register_operator("odd-", {spin, c.epsilon()});
	// equations
	boot.add_equation(Od, array{pair{"unit", E(B{u, s, s, s, s, Od})}, pair{"scalar", E(0, 0, B{e, s, s, s, s, Od})},
	                            pair{"T", E(0, 0, B{T, s, s, s, s, Od})}, pair{"even", E(0, 0, G{s, s, s, s, Od})}});
	boot.add_equation(Od, array{pair{"unit", E(B{u, e, e, e, e, Od})}, pair{"scalar", E(1, 1, B{e, e, e, e, e, Od})},
	                            pair{"T", E(1, 1, B{T, e, e, e, e, Od})}, pair{"even", E(1, 1, G{e, e, e, e, Od})}});
	boot.add_equation(Od, array{pair{"scalar", E(1, 1, B{s, e, s, e, s, Od})}, pair{"odd+", E(G{e, s, e, s, Od})},
	                            pair{"odd-", E(G{e, s, e, s, Od})}});
	boot.add_equation(Od, array{pair{"unit", E(B{u, e, e, s, s, Od})}, pair{"scalar", E(0, 0, B{s, e, s, s, e, Od})},
	                            pair{"scalar", E(0, 1, B{e, e, e, s, s, Od})}, pair{"T", E(0, 1, B{T, e, e, s, s, Od})},
	                            pair{"odd+", E(G{e, s, s, e, Od})}, pair{"odd-", E(R(-1), G{e, s, s, e, Od})},
	                            pair{"even", E(0, 1, G{e, e, s, s, Od})}});
	boot.add_equation(Ev,
	                  array{pair{"unit", E(B{u, e, e, s, s, Ev})}, pair{"scalar", E(0, 0, R(-1), B{s, e, s, s, e, Ev})},
	                        pair{"scalar", E(0, 1, B{e, e, e, s, s, Ev})}, pair{"T", E(0, 1, B{T, e, e, s, s, Ev})},
	                        pair{"odd+", E(R(-1), G{e, s, s, e, Ev})}, pair{"odd-", E(G{e, s, s, e, Ev})},
	                        pair{"even", E(0, 1, G{e, e, s, s, Ev})}});
	auto root = fs::current_path() / ("mixed-ising-" + s.delta().str(8) + "-" + e.delta().str(8));
	auto pmp = boot.create_pmp("unit", [numax](auto spin) { return numax + std::min(numax, spin) / 2; });
	move(pmp).create_input().write_all(root);
}

void single_ising(const Context<R>& c, const Op& s, const Op& e, uint32_t numax = 20, uint32_t maxspin = 24)
{
	constexpr FunctionSymmetry Od = FunctionSymmetry::Odd;
	using E = qboot::Entry<R>;
	using B = Block;
	using G = GGBlock;
	qboot::BootstrapEquation<R> boot(c, numax, array{pair{"unit", 1u}, pair{"e", 1u}, pair{"T", 1u}, pair{"even", 1u}});
	Op u{c.epsilon()}, T{2, c.epsilon()};
	boot.register_operator("even", {0, c.epsilon(), R("1.409"), R("1.4135")});
	boot.register_operator("even", {0, c.epsilon(), R(6)});
	boot.register_operator("even", {2, c.epsilon(), R(5)});
	for (uint32_t spin = 4; spin <= maxspin; spin += 2) boot.register_operator("even", {spin, c.epsilon()});
	boot.add_equation(Od, array{pair{"unit", E(B{u, s, s, s, s, Od})}, pair{"e", E(B{e, s, s, s, s, Od})},
	                            pair{"T", E(B{T, s, s, s, s, Od})}, pair{"even", E(G{s, s, s, s, Od})}});
	auto root = fs::current_path() / ("single-ising-" + s.delta().str(8) + "-" + e.delta().str(8));
	auto pmp = boot.create_pmp("unit", [numax](auto spin) { return numax + std::min(numax, spin) / 2; });
	move(pmp).create_input().write_all(root);
}

class TestScale : public qboot::ScaleFactor<R>
{
	static constexpr uint32_t deg_ = 4;
	Vector<R> xs_;

public:
	TestScale() : xs_(qboot::sample_points<R>(deg_)) {}
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
	[[nodiscard]] Polynomial<R> bilinear_base(uint32_t m) override
	{
		assert(m <= deg_ / 2);
		switch (m)
		{
		case 0: return {R(1)};                  // 1
		case 1: return {R(-1), R(1)};           // -1 + x
		default: return {R(1), R(-2), R(0.5)};  // 1 - 2 x + x ^ 2 / 2
		}
	}
	[[nodiscard]] Vector<Polynomial<R>> bilinear_bases() override
	{
		Vector<Polynomial<R>> q(deg_ / 2 + 1);
		for (uint32_t i = 0; i <= deg_ / 2; ++i) q[i] = bilinear_base(i);
		return q;
	}
};
TestScale::~TestScale() = default;

void test_sdpb()
{
	// https://github.com/davidsd/sdpb/blob/master/test/test.xml
	// {1 + x ^ 4, x ^ 2 + x ^ 4 / 12}
	Vector<Polynomial<R>> elm{{R(1), R(0), R(0), R(0), R(1)}, {R(0), R(0), R(1), R(0), R(1) / 12}};

	PolynomialProgramming<R> prg(1);
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
	constexpr uint32_t n_Max = 400, lambda = 15, dim_ = 3, maxdim = 10, maxspin = 28;
	[[maybe_unused]] constexpr uint32_t numax = 16;
	if (argc > 1)
	{
		assert(argc == 3);
		unique_ptr<char*[]> args(argv);
		auto d_s = R(args[1]);
		auto d_e = R(args[2]);
		Context<R> c(n_Max, lambda, dim_);
		Op sigma(d_s, 0, c.epsilon());
		Op epsilon(d_e, 0, c.epsilon());
		single_ising(c, sigma, epsilon, numax, maxspin);
		args.release();
		return 0;
	}
	R very_small = R("3e-568");
	R d12 = R::sqrt(3), d34 = R::sqrt(5) - 1;
	R S = (d34 - d12) / 2, P = -d12 * d34 / 2, d23h = R(0.7);
	auto c2 = qboot2::context_construct(n_Max, R::prec, lambda);
	for (uint32_t dim = 3; dim <= maxdim; dim += 2)
	{
		Context<R> c(n_Max, lambda, dim);
		for (uint32_t spin = 0; spin <= maxspin; ++spin)
		{
			test_op(c, c2, Op(c.unitarity_bound(spin) + R::sqrt(2), spin, c.epsilon()), d12, d34, very_small);
			test_op(c, c2, Op(c.unitarity_bound(spin), spin, c.epsilon()), d12, d34, very_small);
		}
		test_op(c, c2, Op(c.unitarity_bound(0) + R("0.5"), 0, c.epsilon()), d12, d34, very_small);
		test_h(c, c2, d12, d34, very_small);
	}
	qboot2::clear_cb_context(&c2);
	return 0;
}
