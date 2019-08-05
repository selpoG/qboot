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
#include "polynomial.hpp"
#include "real.hpp"
#include "real_io.hpp"

using R = mpfr::real<1000, MPFR_RNDN>;
static R very_small = R("5e-568");

static algebra::Vector<R> to_vec(const mpfr_t* a, size_t s);
static algebra::Matrix<R> to_mat(const mpfr_t* a, size_t r, size_t c);
template <size_t r>
static algebra::Vector<R> to_vec(const std::array<mpfr_t, r>& a);
template <size_t r, size_t c>
static algebra::Matrix<R> to_mat(const std::array<std::array<mpfr_t, c>, r>& a);
static void test(const R& eps, const R& delta, const R& d12, const R& d34, size_t n);
static void test(const R& eps, const R& delta, const R& d12, const R& d34, size_t n, int spin);
static void test_rec(size_t nMax, const R& eps, R* delta, const R& d12, const R& d34);
static void test_rec(size_t nMax, const R& eps, R* delta, const R& d12, const R& d34, int spin);
static void test_real(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta,
                      const R& d12, const R& d34, int spin);
static void test_g(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12,
                   const R& d34, int spin);
static void test_hB(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, size_t n, const R& eps, R* delta,
                    const R& d12, const R& d34, int spin);
static void test_h(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12,
                   const R& d34);

algebra::Vector<R> to_vec(const mpfr_t* a, size_t s)
{
	algebra::Vector<R> v(s);
	for (size_t i = 0; i < s; i++) v[i] = R(a[i]);
	return v;
}

algebra::Matrix<R> to_mat(const mpfr_t* a, size_t r, size_t c)
{
	algebra::Matrix<R> v(r, c);
	for (size_t i = 0; i < r; i++)
		for (size_t j = 0; j < c; j++) v.get(i, j) = R(a[i * c + j]);
	return v;
}

template <size_t r>
algebra::Vector<R> to_vec(const std::array<mpfr_t, r>& a)
{
	algebra::Vector<R> v(r);
	for (size_t i = 0; i < r; i++) v[i] = R(a.at(i));
	return v;
}

template <size_t r, size_t c>
algebra::Matrix<R> to_mat(const std::array<std::array<mpfr_t, c>, r>& a)
{
	algebra::Matrix<R> m(r, c);
	for (size_t i = 0; i < r; i++)
		for (size_t j = 0; j < c; j++) m.get(i, j) = R(a.at(i).at(j));
	return m;
}

void test(const R& eps, const R& delta, const R& d12, const R& d34, size_t n)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b;
	auto p = qboot::_get_rec_coeffs(eps, R(0), delta, S, P);
	std::array<std::array<mpfr_t, 4>, 6> q{};
	qboot2::initialize_spin_zero_coeffs_folder(q, 1000);
	qboot2::set_zero_spin_rec_coeffs(q, eps._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto err = (to_mat(q) - p).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	auto v = qboot::_evaluate_at_n(p, n);
	std::array<mpfr_t, 6> w{};
	for (size_t i = 0; i < 6; i++) mpfr_init2(w[i], 1000);
	qboot2::spin_zero_evaluate_at_n(w, q, long(n), MPFR_RNDN);
	err = (to_vec(w) - v).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	qboot2::deallocate_spin_zero_coeffs_folder(q);
}

void test(const R& eps, const R& delta, const R& d12, const R& d34, size_t n, int spin)
{
	if (spin == 0)
	{
		test(eps, delta, d12, d34, n);
		return;
	}
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::_get_rec_coeffs(eps, ell, delta, S, P);
	std::array<std::array<mpfr_t, 5>, 8> q{};
	qboot2::initialize_spin_nonzero_coeffs_folder(q, 1000);
	qboot2::set_nonzero_spin_rec_coeffs(q, eps._x, ell._x, delta._x, S._x, P._x, 1000, MPFR_RNDN);
	auto err = (to_mat(q) - p).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	auto v = qboot::_evaluate_at_n(p, n);
	std::array<mpfr_t, 8> w{};
	for (size_t i = 0; i < 8; i++) mpfr_init2(w[i], 1000);
	qboot2::spin_nonzero_evaluate_at_n(w, q, long(n), 1000, MPFR_RNDN);
	err = (to_vec(w) - v).norm();
	if (err > very_small)
	{
		std::cout << "err = " << err << std::endl;
		assert(false);
	}
	qboot2::deallocate_spin_nonzero_coeffs_folder(q);
}

void test_rec(size_t nMax, const R& eps, R* delta, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b;
	auto p = qboot::_recursion_vector(nMax, eps, R(0), *delta, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::recursionSpinZeroVector(static_cast<unsigned long>(nMax), eps._x, delta->_x,
	                                                                S._x, P._x, 1000, MPFR_RNDN));
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

void test_rec(size_t nMax, const R& eps, R* delta, const R& d12, const R& d34, int spin)
{
	if (spin == 0)
	{
		test_rec(nMax, eps, delta, d12, d34);
		return;
	}
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::_recursion_vector(nMax, eps, ell, *delta, S, P);
	std::unique_ptr<mpfr_t[]> q_ptr(qboot2::recursionNonZeroVector(static_cast<unsigned long>(nMax), eps._x, ell._x,
	                                                               delta->_x, S._x, P._x, 1000, MPFR_RNDN));
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

void test_real(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12,
               const R& d34, int spin)
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

void test_g(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, R* delta, const R& d12,
            const R& d34, int spin)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::gBlock_full(eps, ell, *delta, S, P, cb1);
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

void test_hB(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, size_t n, const R& eps, R* delta,
             const R& d12, const R& d34, int spin)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b, P = 2 * a * b, ell = R(spin);
	auto p = qboot::hBlock_times_rho_n(n, eps, ell, *delta, S, P, cb1);
	std::unique_ptr<mpfr_t[]> q_ptr(
	    qboot2::hBlock_times_rho_n(static_cast<unsigned long>(n), eps._x, ell._x, delta->_x, S._x, P._x, cb2));
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

void test_h(const qboot::cb_context<R>& cb1, const qboot2::cb_context& cb2, const R& eps, const R& d12, const R& d34)
{
	R a = -d12 / 2, b = d34 / 2, S = a + b;
	auto p = qboot::h_asymptotic(eps, S, cb1);
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
	auto c = qboot::cb_context<R>(250, 10);
	auto c2 = qboot2::context_construct(250, 1000, 10);
	R delta = mpfr::sqrt(R(2)), d12 = mpfr::sqrt(R(3)), d34 = mpfr::sqrt(R(5)) - 1;
	auto hoge =
	    algebra::Polynomial<R>::linear_power(R(1), R(1), 5) * algebra::Polynomial<R>::linear_power(R(1), R(-1), 7);
	auto fuga =
	    algebra::Polynomial<R>::linear_power(R(1), R(1), 7) * algebra::Polynomial<R>::linear_power(R(1), R(-1), 5);
	std::cout << hoge << std::endl;
	std::cout << fuga << std::endl;
	std::cout << hoge + fuga << std::endl;
	for (int dim = 3; dim < 8; dim += 2)
	{
		R eps = R(dim) / 2 - 1;
		for (int spin = 0; spin < 6; spin++)
		{
			for (size_t n = 0; n < 6; n++)
			{
				test(eps, delta, d12, d34, n, spin);
				test_hB(c, c2, n, eps, &delta, d12, d34, spin);
			}
			test_real(c, c2, eps, &delta, d12, d34, spin);
			test_rec(c.n_Max, eps, &delta, d12, d34, spin);
			test_g(c, c2, eps, &delta, d12, d34, spin);
		}
		test_h(c, c2, eps, d12, d34);
	}
	for (int dim = 3; dim < 8; dim += 2)
		for (int spin = 0; spin < 6; spin++)
		{
			R eps = R(dim) / 2 - 1, unit = spin == 0 ? eps : 2 * eps + spin;
			for (size_t n = 0; n < 6; n++)
			{
				test(eps, unit, d12, d34, n, spin);
				test_hB(c, c2, n, eps, &unit, d12, d34, spin);
			}
			test_real(c, c2, eps, &unit, d12, d34, spin);
			test_rec(c.n_Max, eps, &unit, d12, d34, spin);
			test_g(c, c2, eps, &unit, d12, d34, spin);
		}
	for (int dim = 3; dim < 8; dim += 2)
	{
		R eps = R(dim) / 2 - 1, unit = eps + R("0.5");
		for (size_t n = 0; n < 6; n++)
		{
			test(eps, unit, d12, d34, n, 0);
			test_hB(c, c2, n, eps, &unit, d12, d34, 0);
		}
		test_real(c, c2, eps, &unit, d12, d34, 0);
		test_rec(c.n_Max, eps, &unit, d12, d34, 0);
		test_g(c, c2, eps, &unit, d12, d34, 0);
	}
	return 0;
}
