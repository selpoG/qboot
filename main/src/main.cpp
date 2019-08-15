#include <iomanip>
#include <iostream>

#include "chol_and_inverse.hpp"
#include "complex_function.hpp"
#include "context_variables.hpp"
#include "damped_rational.hpp"
#include "hor_formula.hpp"
#include "hor_recursion.hpp"
#include "integral_decomp.hpp"
#include "io.hpp"
#include "matrix.hpp"
#include "multi_array.hpp"
#include "my_hash.hpp"
#include "partial_fraction.hpp"
#include "pole_data.hpp"
#include "polynomial.hpp"
#include "primary_op.hpp"
#include "real.hpp"
#include "real_io.hpp"

using algebra::Vector, algebra::Matrix, algebra::Tensor, algebra::Polynomial;
using qboot::h_asymptotic, qboot::gBlock_full, qboot::Context, qboot::PrimaryOperator, qboot::ComplexFunction;
using std::array;

using R = mpfr::real<1000, MPFR_RNDN>;
static R very_small = R("5e-568");

static Vector<R> to_vec(const mpfr_t* a, uint32_t s);
[[maybe_unused]] static Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c);
static ComplexFunction<R> to_hol(const mpfr_t* a, uint32_t l);
static void test_rec(const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_real(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_g(const qboot2::cb_context& cb, const PrimaryOperator<R>& op, const R& d12, const R& d34);
static void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12, const R& d34);

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
	ComplexFunction f(l);
	uint32_t i = 0;
	for (uint32_t dy = 0; dy <= l / 2; dy++)
		for (uint32_t dx = 0; dx + 2 * dy <= l; dx++) f.get(dx, dy) = R(a[i++]);
	return f;
}

void test_rec(const PrimaryOperator<R>& op, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(op.spin()), delta = op.delta();
	auto nMax = op.context().n_Max;
	auto p = qboot::power_series_in_rho(op, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr =
	    op.spin() == 0
	        ? qboot2::recursionSpinZeroVector(nMax, op.epsilon()._x, delta._x, S._x, P._x, 1000, MPFR_RNDN)
	        : qboot2::recursionNonZeroVector(nMax, op.epsilon()._x, ell._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto q = to_vec(q_ptr.get(), p.size());
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
	auto p = qboot::hBlock_powered(op, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::real_axis_result(op.epsilon()._x, ell._x, delta._x, S._x, P._x, cb));
	auto q = to_vec(q_ptr.get(), p.size());
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
	auto p = gBlock_full(op, S, P);
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

void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = h_asymptotic(S, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::h_asymptotic(eps._x, S._x, cb2));
	auto q = to_vec(q_ptr.get(), p.size());
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
	constexpr uint32_t n_Max = 50, lambda = 6, dim_ = 3;
	[[maybe_unused]] constexpr uint32_t numax = 2;
	std::map<uint32_t, std::optional<Context<R>>> cs;
	for (uint32_t dim = 3; dim < 10; dim += 2) cs.emplace(dim, Context<R>(n_Max, lambda, dim));
	R delta = mpfr::sqrt(R(2)), d12 = mpfr::sqrt(R(3)), d34 = mpfr::sqrt(R(5)) - 1;
	auto c2 = qboot2::context_construct(n_Max, R::prec, lambda);
	{
		const auto& c = cs[dim_].value();
		auto op = c.get_primary(R(0.81), 0);
		std::cout << c.F_minus_matrix(R(0.7)) * gBlock_full(op, R(0), R(0)).as_vector() << std::endl;
		qboot::RationalApproxData<R> ag(numax, 0, c, d12, d34);
		std::cout << ag.approx_g() << std::endl;
	}
	for (uint32_t dim = 3; dim < 10; dim += 2)
	{
		const auto& c = cs[dim].value();
		for (uint32_t spin = 0; spin < 10; ++spin)
		{
			auto op = c.get_primary(delta, spin);
			test_rec(op, d12, d34);
			test_real(c2, op, d12, d34);
			test_g(c2, op, d12, d34);
		}
		test_h(c, c2, c.epsilon, d12, d34);
	}
	std::cout << "check divergent cases" << std::endl;
	for (uint32_t dim = 3; dim < 10; dim += 2)
	{
		const auto& c = cs[dim].value();
		for (uint32_t spin = 0; spin < 10; ++spin)
		{
			auto op = c.get_primary(c.unitary_bound(spin), spin);
			test_real(c2, op, d12, d34);
			test_rec(op, d12, d34);
			test_g(c2, op, d12, d34);
		}
	}
	for (uint32_t dim = 3; dim < 10; dim += 2)
	{
		const auto& c = cs[dim].value();
		auto op = c.get_primary(c.epsilon + R("0.5"), 0);
		test_real(c2, op, d12, d34);
		test_rec(op, d12, d34);
		test_g(c2, op, d12, d34);
	}
	qboot2::clear_cb_context(c2);
	return 0;
}
