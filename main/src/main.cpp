#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "block.hpp"
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
#include "primary_op.hpp"
#include "real.hpp"
#include "real_function.hpp"
#include "real_io.hpp"
#include "sdpb_input.hpp"

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using qboot::h_asymptotic, qboot::gBlock, qboot::Context, algebra::ComplexFunction, algebra::RealFunction,
    qboot::RationalApproxData;
using std::array;

using R = mpfr::real<1000, MPFR_RNDN>;
using Op = qboot::PrimaryOperator<R>;
using GOp = qboot::GeneralPrimaryOperator<R>;
using GBlock = qboot::ConformalBlock<R, GOp>;
static R very_small = R("2e-591");

[[maybe_unused]] static Vector<R> to_vec(const mpfr_t* a, uint32_t s);
[[maybe_unused]] static Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c);
static ComplexFunction<R> to_hol(const mpfr_t* a, uint32_t l);
static RealFunction<R> to_func(const mpfr_t* a, uint32_t l);
static void test_rec(uint32_t nMax, const Op& op, const R& d12, const R& d34);
static void test_real(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34);
static void test_g(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34);
static void test_op(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34);
static void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34);
[[maybe_unused]] static void solve_ising(const Context<R>& c, const R& ds, const R& de, uint32_t numax,
                                         uint32_t maxspin);
[[maybe_unused]] static void test_sdpb();

template <class T>
void check(const T& p, const T& q)
{
	auto err = (p - q).norm();
	if (err > very_small)
	{
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << (p - q) << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}
template <class T, class CallBack_T>
void check(const T& p, const T& q, CallBack_T on_err)
{
	auto err = (p - q).norm();
	if (err > very_small)
	{
		on_err();
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << (p - q) << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}

Vector<R> to_vec(const mpfr_t* a, uint32_t s)
{
	Vector<R> v(s);
	for (uint32_t i = 0; i < s; ++i) v[i] = R(a[i]);
	return v;
}

Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c)
{
	Matrix<R> v(r, c);
	for (uint32_t i = 0; i < r; ++i)
		for (uint32_t j = 0; j < c; ++j) v.at(i, j) = R(a[i * c + j]);
	return v;
}

ComplexFunction<R> to_hol(const mpfr_t* a, uint32_t l)
{
	ComplexFunction<R> f(l);
	uint32_t i = 0;
	for (uint32_t dy = 0; dy <= l / 2; ++dy)
		for (uint32_t dx = 0; dx + 2 * dy <= l; ++dx) f.at(dx, dy) = R(a[i++]);
	return f;
}

RealFunction<R> to_func(const mpfr_t* a, uint32_t l)
{
	RealFunction<R> f(l);
	for (uint32_t k = 0; k <= l; ++k) f.at(k) = R(a[k]);
	return f;
}

void test_rec(uint32_t nMax, const Op& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = qboot::hBlock_shifted<R>(op, S, P, nMax);
	std::unique_ptr<mpfr_t[]> q_ptr =
	    op.spin() == 0
	        ? qboot2::recursionSpinZeroVector(nMax, op.epsilon()._x, delta._x, S._x, P._x, 1000, MPFR_RNDN)
	        : qboot2::recursionNonZeroVector(nMax, op.epsilon()._x, ell._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto q = to_func(q_ptr.get(), p.lambda());
	check(p, q, [&]() {
		std::cout << "test_rec(nMax=" << nMax << "op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")"
		          << std::endl;
	});
}

void test_real(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = qboot::gBlock_real<R>(op, S, P, c);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::real_axis_result(op.epsilon()._x, ell._x, delta._x, S._x, P._x, cb));
	auto q = to_func(q_ptr.get(), p.lambda());
	check(p, q, [&]() {
		std::cout << "test_real(c=" << c.str() << ", cb, op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")"
		          << std::endl;
	});
}

void test_g(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = gBlock<R>(op, S, P, c);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::gBlock_full(op.epsilon()._x, ell._x, delta._x, S._x, P._x, cb));
	auto q = to_hol(q_ptr.get(), cb.lambda);
	check(p, q, [&]() {
		std::cout << "test_g(c=" << c.str() << ", cb, op=" << op.str() << ", d12=" << d12 << ", d34=" << d34 << ")"
		          << std::endl;
	});
}

void test_op(const Context<R>& c, const qboot2::cb_context& cb, const Op& op, const R& d12, const R& d34)
{
	test_rec(c.n_Max, op, d12, d34);
	test_real(c, cb, op, d12, d34);
	test_g(c, cb, op, d12, d34);
}

void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = h_asymptotic(S, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::h_asymptotic(cb1.epsilon._x, S._x, cb2));
	auto q = to_func(q_ptr.get(), p.lambda());
	check(p, q, [&]() {
		std::cout << "test_h(cb1=" << cb1.str() << ", cb2, d12=" << d12 << ", d34=" << d34 << ")" << std::endl;
	});
}

void solve_ising(const Context<R>& c, const R& ds, const R& de, uint32_t numax = 20, uint32_t maxspin = 24)
{
	for (uint32_t spin = 0; spin <= maxspin; spin += 2)
	{
		R gap = spin == 0 ? de : c.unitary_bound(spin);
		auto op = GOp(spin, c.epsilon);
		auto block = GBlock(op, ds, ds, ds, ds, algebra::FunctionSymmetry::Odd);
		RationalApproxData<R> ag(numax + std::min(numax, spin) / 2, spin, c, ds, ds, ds, ds);
		auto sp = ag.sample_points();
		auto q = ag.get_bilinear_basis(gap);
		Vector<ComplexFunction<R>> bls(sp.size());
		Vector<R> scales(sp.size());
		for (uint32_t i = 0; i < sp.size(); ++i)
		{
			R delta = gap + sp[i];
			scales[i] = ag.get_scale(delta);
			bls[i] = c.evaluate(block, delta) / scales[i];
		}
		std::cout << "F_{-} = " << algebra::polynomial_interpolate(bls, sp) << std::endl;
		std::cout << "basis = " << q << std::endl;
		std::cout << "scales = " << scales << std::endl;
	}
}

void test_sdpb()
{
	// https://github.com/davidsd/sdpb/blob/master/test/test.xml
	constexpr uint32_t deg = 4, d0 = deg / 2, d1 = deg == 0 ? 0 : (deg - 1) / 2;
	Vector<Polynomial<R>> bilinear{{R(1)}, {R(-1), R(1)}, {R(1), R(-2), R(0.5)}};
	Vector<Polynomial<R>> elm{{R(1), R(0), R(0), R(0), R(1)}, {R(0), R(0), R(1), R(0), R(1) / 12}};
	auto sample_points = qboot::sample_points<R>(deg);
	Vector<R> sample_scale(deg + 1);
	for (uint32_t i = 0; i <= deg; ++i) sample_scale[i] = mpfr::exp(-sample_points[i]);
	Vector<R> c(deg + 1);
	for (uint32_t i = 0; i <= deg; ++i) c.at(i) = elm[0].eval(sample_points[i]) * sample_scale[i];
	Matrix<R> B(deg + 1, elm.size() - 1);
	for (uint32_t i = 0; i <= deg; ++i) B.at(i, 0) = -elm[1].eval(sample_points[i]) * sample_scale[i];
	Matrix<R> qb0(d0 + 1, deg + 1);
	for (uint32_t i = 0; i <= d0; ++i)
		for (uint32_t j = 0; j <= deg; ++j)
			qb0.at(i, j) = bilinear[i].eval(sample_points[j]) * mpfr::sqrt(sample_scale[j]);
	Matrix<R> qb1(d1 + 1, deg + 1);
	for (uint32_t i = 0; i <= d1; ++i)
		for (uint32_t j = 0; j <= deg; ++j)
			qb1.at(i, j) = bilinear[i].eval(sample_points[j]) * mpfr::sqrt(sample_scale[j] * sample_points[j]);
	constexpr uint32_t num_cons = 1;
	qboot::SDPBInput<R> sdpb(R(0), Vector<R>{R(-1)}, num_cons);
	sdpb.register_constraint(0, {1, 4, std::move(B), std::move(c), {std::move(qb0), std::move(qb1)}});
	sdpb.write_all(std::cout);
}

int main()
{
	constexpr uint32_t n_Max = 100, lambda = 5, dim_ = 3, maxdim = 10, maxspin = 10;
	[[maybe_unused]] constexpr uint32_t numax = 5;
	std::map<uint32_t, std::optional<Context<R>>> cs;
	for (uint32_t dim = 3; dim <= maxdim; dim += 2) cs.emplace(dim, Context<R>(n_Max, lambda, dim));
	R d12 = mpfr::sqrt(R(3)), d34 = mpfr::sqrt(R(5)) - 1;
	R S = (d34 - d12) / 2, P = -d12 * d34 / 2, d23h = R(0.7);
	auto c2 = qboot2::context_construct(n_Max, R::prec, lambda);
	{
		const auto& c = cs[dim_].value();
		R d_s = R(0.5181475), d_e = R(1.412617);
		for (uint32_t spin = 0; spin <= 2; ++spin)
		{
			R gap = R(spin == 0 ? 3 : spin + 1);
			auto op = GOp(spin, c.epsilon);
			auto block = GBlock(op, d_s, d_e, d_s, d_e, algebra::FunctionSymmetry::Odd);
			std::cout << "block = " << block.str() << std::endl;
			RationalApproxData<R> ag(numax, spin, c, d_s, d_e, d_s, d_e);
			auto sp = ag.sample_points();
			auto q = ag.get_bilinear_basis(gap);
			const auto& pol = ag.get_poles();
			std::cout << "pol = " << pol << std::endl;
			Vector<ComplexFunction<R>> bls(sp.size());
			Vector<R> scales(sp.size());
			std::cout << "ps = " << gap << " + " << sp << std::endl;
			for (uint32_t i = 0; i < sp.size(); ++i)
			{
				auto delta = gap + sp[i];
				scales[i] = ag.get_scale(delta);
				bls[i] = c.evaluate(block, delta) / scales[i];
			}
			const auto& ip = algebra::polynomial_interpolate(bls, sp);
			std::cout << "F_{-} = " << ip << std::endl;
			std::cout << "basis = " << q << std::endl;
			std::cout << "scales = " << scales << std::endl;
			std::cout << "error = " << (evals(ip, sp) - bls).norm() << std::endl;
		}
	}
	for (uint32_t dim = 3; dim <= maxdim; dim += 2)
	{
		const auto& c = cs[dim].value();
		for (uint32_t spin = 0; spin <= maxspin; ++spin)
		{
			test_op(c, c2, Op(c.unitary_bound(spin) + mpfr::sqrt(R(2)), spin, c.epsilon), d12, d34);
			test_op(c, c2, Op(c.unitary_bound(spin), spin, c.epsilon), d12, d34);
		}
		test_op(c, c2, Op(c.unitary_bound(0) + R("0.5"), 0, c.epsilon), d12, d34);
		test_h(c, c2, d12, d34);
	}
	qboot2::clear_cb_context(c2);
	return 0;
}
