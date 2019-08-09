#include <iomanip>
#include <iostream>

#include "chol_and_inverse.hpp"
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
#include "real.hpp"
#include "real_io.hpp"

using algebra::Vector, algebra::Matrix, algebra::Tensor, algebra::Polynomial;
using qboot::hBlock_times_rho_n, qboot::h_asymptotic, qboot::gBlock_full, qboot::Context;
using std::array;

using R = mpfr::real<1000, MPFR_RNDN>;
static R very_small = R("5e-568");

static Vector<R> to_vec(const mpfr_t* a, uint32_t s);
[[maybe_unused]] static Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c);
template <size_t r>
static Vector<R> to_vec(const array<mpfr_t, r>& a);
template <size_t r, size_t c>
static Matrix<R> to_mat(const array<array<mpfr_t, c>, r>& a);
static void test(const R& eps, const R& delta, const R& d12, const R& d34, uint32_t n);
static void test(const R& eps, const R& delta, const R& d12, const R& d34, uint32_t n, int32_t spin);
static void test_rec(uint32_t nMax, const R& eps, R* delta, const R& d12, const R& d34);
static void test_rec(uint32_t nMax, const R& eps, R* delta, const R& d12, const R& d34, int32_t spin);
static void test_real(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12,
                      const R& d34, int32_t spin);
static void test_g(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12,
                   const R& d34, int32_t spin);
static void test_hB(const Context<R>& cb1, const qboot2::cb_context& cb2, uint32_t n, const R& eps, R* delta,
                    const R& d12, const R& d34, int32_t spin);
static void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12, const R& d34);

Vector<R> to_vec(const mpfr_t* a, uint32_t s)
{
	Vector<R> v(s);
	for (uint32_t i = 0; i < s; i++) v[i] = R(a[i]);
	return v;
}

Matrix<R> to_mat(const mpfr_t* a, uint32_t r, uint32_t c)
{
	Matrix<R> v(r, c);
	for (uint32_t i = 0; i < r; i++)
		for (uint32_t j = 0; j < c; j++) v.get(i, j) = R(a[i * c + j]);
	return v;
}

template <size_t r>
Vector<R> to_vec(const array<mpfr_t, r>& a)
{
	Vector<R> v(r);
	for (uint32_t i = 0; i < r; i++) v[i] = R(a.at(i));
	return v;
}

template <size_t r, size_t c>
Matrix<R> to_mat(const array<array<mpfr_t, c>, r>& a)
{
	Matrix<R> m(r, c);
	for (uint32_t i = 0; i < r; i++)
		for (uint32_t j = 0; j < c; j++) m.get(i, j) = R(a.at(i).at(j));
	return m;
}

void test(const R& eps, const R& delta, const R& d12, const R& d34, uint32_t n)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b;
	auto p = qboot::_get_rec_coeffs(eps, R(0), delta, S, P);
	array<array<mpfr_t, 4>, 6> q{};
	qboot2::initialize_spin_zero_coeffs_folder(q, 1000);
	qboot2::set_zero_spin_rec_coeffs(q, eps._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto err = (to_mat(q) - p).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	auto v = qboot::_evaluate_at_n(p, n);
	array<mpfr_t, 6> w{};
	for (uint32_t i = 0; i < 6; i++) mpfr_init2(w.at(i), 1000);
	qboot2::spin_zero_evaluate_at_n(w, q, int32_t(n), MPFR_RNDN);
	err = (to_vec(w) - v).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	qboot2::deallocate_spin_zero_coeffs_folder(q);
}

void test(const R& eps, const R& delta, const R& d12, const R& d34, uint32_t n, int32_t spin)
{
	if (spin == 0)
	{
		test(eps, delta, d12, d34, n);
		return;
	}
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::_get_rec_coeffs(eps, ell, delta, S, P);
	array<array<mpfr_t, 5>, 8> q{};
	qboot2::initialize_spin_nonzero_coeffs_folder(q, 1000);
	qboot2::set_nonzero_spin_rec_coeffs(q, eps._x, ell._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto err = (to_mat(q) - p).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	auto v = qboot::_evaluate_at_n(p, n);
	array<mpfr_t, 8> w{};
	for (uint32_t i = 0; i < 8; i++) mpfr_init2(w.at(i), 1000);
	qboot2::spin_nonzero_evaluate_at_n(w, q, int32_t(n), MPFR_RNDN);
	err = (to_vec(w) - v).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	qboot2::deallocate_spin_nonzero_coeffs_folder(q);
}

void test_rec(uint32_t nMax, const R& eps, R* delta, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b;
	auto p = qboot::_recursion_vector(nMax, eps, R(0), *delta, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(
	    qboot2::recursionSpinZeroVector(nMax, eps._x, delta->_x, S._x, P._x, 1000, MPFR_RNDN));
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

void test_rec(uint32_t nMax, const R& eps, R* delta, const R& d12, const R& d34, int32_t spin)
{
	if (spin == 0)
	{
		test_rec(nMax, eps, delta, d12, d34);
		return;
	}
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::_recursion_vector(nMax, eps, ell, *delta, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(
	    qboot2::recursionNonZeroVector(nMax, eps._x, ell._x, delta->_x, S._x, P._x, 1000, MPFR_RNDN));
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

void test_real(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12, const R& d34,
               int32_t spin)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::_real_axis_result(eps, ell, *delta, S, P, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::real_axis_result(eps._x, ell._x, delta->_x, S._x, P._x, cb2));
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

void test_g(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12, const R& d34,
            int32_t spin)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = gBlock_full(eps, ell, *delta, S, P, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::gBlock_full(eps._x, ell._x, delta->_x, S._x, P._x, cb2));
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

void test_hB(const Context<R>& cb1, const qboot2::cb_context& cb2, uint32_t n, const R& eps, R* delta, const R& d12,
             const R& d34, int32_t spin)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = hBlock_times_rho_n(n, eps, ell, *delta, S, P, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::hBlock_times_rho_n(n, eps._x, ell._x, delta->_x, S._x, P._x, cb2));
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

void test_h(const Context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = h_asymptotic(eps, S, cb1);
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
	auto c = Context<R>(250, 8, 3);
	auto c2 = qboot2::context_construct(250, 1000, 8);
	R delta = mpfr::sqrt(R(2)), d12 = mpfr::sqrt(R(3)), d34 = mpfr::sqrt(R(5)) - 1;
	std::cout << c.F_minus_matrix(R(0.7)) * gBlock_full(R(0.5), R(0), R(0.81), R(0), R(0), c) << std::endl;
	for (int32_t dim = 3; dim < 8; dim += 2)
	{
		R eps = R(dim) / 2 - 1;
		for (int32_t spin = 0; spin < 6; spin++)
		{
			test_rec(c.n_Max, eps, &delta, d12, d34, spin);
			for (uint32_t n = 0; n < 6; n++)
			{
				test(eps, delta, d12, d34, n, spin);
				test_hB(c, c2, n, eps, &delta, d12, d34, spin);
			}
			test_real(c, c2, eps, &delta, d12, d34, spin);
			test_g(c, c2, eps, &delta, d12, d34, spin);
		}
		test_h(c, c2, eps, d12, d34);
	}
	for (int32_t dim = 3; dim < 8; dim += 2)
		for (int32_t spin = 0; spin < 6; spin++)
		{
			R eps = R(dim) / 2 - 1, unit = spin == 0 ? eps : 2 * eps + spin;
			for (uint32_t n = 0; n < 6; n++)
			{
				test(eps, unit, d12, d34, n, spin);
				test_hB(c, c2, n, eps, &unit, d12, d34, spin);
			}
			test_real(c, c2, eps, &unit, d12, d34, spin);
			test_rec(c.n_Max, eps, &unit, d12, d34, spin);
			test_g(c, c2, eps, &unit, d12, d34, spin);
		}
	for (int32_t dim = 3; dim < 8; dim += 2)
	{
		R eps = R(dim) / 2 - 1, unit = eps + R("0.5");
		for (uint32_t n = 0; n < 6; n++)
		{
			test(eps, unit, d12, d34, n, 0);
			test_hB(c, c2, n, eps, &unit, d12, d34, 0);
		}
		test_real(c, c2, eps, &unit, d12, d34, 0);
		test_rec(c.n_Max, eps, &unit, d12, d34, 0);
		test_g(c, c2, eps, &unit, d12, d34, 0);
	}
	qboot2::clear_cb_context(c2);
	return 0;
}
