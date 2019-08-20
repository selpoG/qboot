#include <iomanip>
#include <iostream>

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

using algebra::Vector, algebra::Matrix, algebra::Tensor, algebra::Polynomial;
using qboot::h_asymptotic, qboot::gBlock, qboot::Context, qboot::PrimaryOperator, qboot::ComplexFunction,
    qboot::RealFunction, qboot::RationalApproxData;
using std::array;

using R = mpfr::real<1000, MPFR_RNDN>;
static R very_small = R("5e-568");

[[maybe_unused]] static Vector<R> to_vec(const mpfr_t* a, uint32_t s);
[[maybe_unused]] static Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c);
static ComplexFunction<R> to_hol(const mpfr_t* a, uint32_t l);
static RealFunction<R> to_func(const mpfr_t* a, uint32_t l);
static void test_rec(const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_real(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_g(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_op(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34);

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
		for (uint32_t j = 0; j < c; ++j) v.get(i, j) = R(a[i * c + j]);
	return v;
}

ComplexFunction<R> to_hol(const mpfr_t* a, uint32_t l)
{
	ComplexFunction<R> f(l);
	uint32_t i = 0;
	for (uint32_t dy = 0; dy <= l / 2; ++dy)
		for (uint32_t dx = 0; dx + 2 * dy <= l; ++dx) f.get(dx, dy) = R(a[i++]);
	return f;
}

RealFunction<R> to_func(const mpfr_t* a, uint32_t l)
{
	RealFunction<R> f(l);
	for (uint32_t k = 0; k <= l; ++k) f.get(k) = R(a[k]);
	return f;
}

void test_rec(const PrimaryOperator<R>& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto nMax = op.context().n_Max;
	auto p = qboot::hBlock_shifted(op, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr =
	    op.spin() == 0
	        ? qboot2::recursionSpinZeroVector(nMax, op.epsilon()._x, delta._x, S._x, P._x, 1000, MPFR_RNDN)
	        : qboot2::recursionNonZeroVector(nMax, op.epsilon()._x, ell._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto q = to_func(q_ptr.get(), p.lambda());
	auto err = (q - p).norm();
	if (err > very_small)
	{
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << p - q << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}

void test_real(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = qboot::gBlock_real(op, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::real_axis_result(op.epsilon()._x, ell._x, delta._x, S._x, P._x, cb));
	auto q = to_func(q_ptr.get(), p.lambda());
	auto err = (q - p).norm();
	if (err > very_small)
	{
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << p - q << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}

void test_g(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto p = gBlock(op, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::gBlock_full(op.epsilon()._x, ell._x, delta._x, S._x, P._x, cb));
	auto q = to_hol(q_ptr.get(), cb.lambda);
	auto err = (q - p).norm();
	if (err > very_small)
	{
		std::cout << "cb1 = " << op.context().str() << ", eps = " << op.epsilon() << ", delta = " << delta
		          << ", d12 = " << d12 << ", d34 = " << d34 << ", spin = " << op.spin() << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << p - q << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}

void test_op(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34)
{
	test_rec(op, d12, d34);
	test_real(cb, op, d12, d34);
	test_g(cb, op, d12, d34);
}

void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = h_asymptotic(S, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::h_asymptotic(cb1.epsilon._x, S._x, cb2));
	auto q = to_func(q_ptr.get(), p.lambda());
	auto err = (q - p).norm();
	if (err > very_small)
	{
		std::cout << "p = " << p << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "p - q = " << p - q << std::endl;
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
}

int main()
{
	constexpr uint32_t n_Max = 100, lambda = 4, dim_ = 3;
	[[maybe_unused]] constexpr uint32_t numax = 5;
	std::map<uint32_t, std::optional<Context<R>>> cs;
	for (uint32_t dim = 3; dim < 10; dim += 2) cs.emplace(dim, Context<R>(n_Max, lambda, dim));
	R d12 = mpfr::sqrt(R(3)), d34 = mpfr::sqrt(R(5)) - 1;
	R S = (d34 - d12) / 2, P = -d12 * d34 / 2, d23h = R(0.7);
	auto c2 = qboot2::context_construct(n_Max, R::prec, lambda);
	{
		const auto& c = cs[dim_].value();
		R d_s = R(0.5181475), d_e = R(1.412617);
		for (uint32_t spin = 0; spin <= 2; spin++)
		{
			R gap = R(spin == 0 ? 3 : spin + 1);
			RationalApproxData<R> ag(numax, spin, c, d_s, d_e, d_s, d_e);
			auto sp = ag.sample_points();
			auto q = ag.get_bilinear_basis(gap);
			const auto& pol = ag.get_poles();
			std::cout << "pol = " << pol << std::endl;
			for (uint32_t i = 0; i < sp.size(); i++)
			{
				R delta = gap + sp[i];
				auto op = c.get_primary(delta, spin);
				std::cout << "op = " << op.str() << std::endl;
				// auto g = gBlock(op, d_s, d_e, d_s, d_e);
				std::cout << "F_{-} = " << c.F_block(op, d_s, d_e, d_s, d_e) << std::endl;
				// std::cout << "F_{+} = " << c.H_block(op, d_s, d_e, d_s, d_e) << std::endl;
				// std::cout << "scale = " << ag.get_scale(delta) << std::endl;
				// Vector<R> q_eval(q.size());
				// for (uint32_t j = 0; j < q.size(); j++) q_eval[j] = q[j].eval(sp[j]);
				// std::cout << "q = " << q_eval << std::endl;
			}
		}
	}
	for (uint32_t dim = 3; dim < 10; dim += 2)
	{
		const auto& c = cs[dim].value();
		for (uint32_t spin = 0; spin < 10; ++spin)
		{
			test_op(c2, c.get_primary(c.unitary_bound(spin) + mpfr::sqrt(R(2)), spin), d12, d34);
			test_op(c2, c.get_primary(c.unitary_bound(spin), spin), d12, d34);
		}
		test_op(c2, c.get_primary(c.epsilon + R("0.5"), 0), d12, d34);
		test_h(c, c2, d12, d34);
	}
	qboot2::clear_cb_context(c2);
	return 0;
}
