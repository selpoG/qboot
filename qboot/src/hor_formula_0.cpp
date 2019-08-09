#include "hor_formula.hpp"
using std::array;

namespace qboot2
{
	void initialize_spin_zero_coeffs_folder(array<array<mpfr_t, 4>, 6>& a, mpfr_prec_t prec)
	{
		for (uint32_t i = 0; i <= 5; i++)
		{
			for (uint32_t j = 0; j <= 3; j++) { mpfr_init2(a.at(i).at(j), prec); }
		}
	}

	void deallocate_spin_zero_coeffs_folder(array<array<mpfr_t, 4>, 6>& a)
	{
		for (uint32_t i = 0; i <= 5; i++)
		{
			for (uint32_t j = 0; j <= 3; j++) { mpfr_clear(a.at(i).at(j)); }
		}
	}

	void spin_zero_evaluate_at_n(array<mpfr_t, 6>& a, array<array<mpfr_t, 4>, 6>& rho, int32_t n, mpfr_rnd_t rnd)
	{
		for (uint32_t i = 0; i <= 5; i++)
		{
			mpfr_set(a.at(i), rho.at(i).at(3), rnd);
			for (uint32_t j = 0; j <= 2; j++)
			{
				mpfr_mul_si(a.at(i), a.at(i), n, rnd);
				mpfr_add(a.at(i), a.at(i), rho.at(i).at(2 - j), rnd);
			}
		}
	}

	void set_zero_spin_rec_coeffs(array<array<mpfr_t, 4>, 6>& a, const mpfr_t& epsilon, const mpfr_t& Delta,
	                              const mpfr_t& S, const mpfr_t& P, mpfr_prec_t prec, mpfr_rnd_t rnd)
	{
		mpfr_t _t_0;
		mpfr_init2(_t_0, prec);
		mpfr_t _t_1;
		mpfr_init2(_t_1, prec);
		mpfr_t _t_2;
		mpfr_init2(_t_2, prec);
		mpfr_t _t_3;
		mpfr_init2(_t_3, prec);
		mpfr_set_si(a[0][0], 0, MPFR_RNDN);
		mpfr_neg(_t_0, Delta, rnd);
		mpfr_sub_si(_t_2, Delta, 1, rnd);
		mpfr_add(_t_0, _t_0, epsilon, rnd);
		mpfr_mul_si(_t_1, _t_2, -2, rnd);
		mpfr_add_si(_t_0, _t_0, 1, rnd);
		mpfr_mul(a[0][1], _t_1, _t_0, rnd);
		mpfr_mul_si(_t_0, Delta, 3, rnd);
		mpfr_mul_si(_t_1, epsilon, 2, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[0][2], _t_0, 3, rnd);
		mpfr_set_si(a[0][3], 1, MPFR_RNDN);
		mpfr_mul_si(_t_0, Delta, 4, rnd);
		mpfr_mul(_t_0, _t_0, S, rnd);
		mpfr_add(_t_0, _t_0, Delta, rnd);
		mpfr_mul_si(_t_2, S, 8, rnd);
		mpfr_mul_si(_t_1, Delta, -2, rnd);
		mpfr_mul_si(_t_3, epsilon, 2, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_2, P, 4, rnd);
		mpfr_add(_t_1, _t_1, _t_3, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_add_si(_t_1, _t_1, 3, rnd);
		mpfr_sub_si(_t_0, _t_0, 2, rnd);
		mpfr_mul(a[1][0], _t_1, _t_0, rnd);
		mpfr_mul_si(_t_0, Delta, 2, rnd);
		mpfr_mul_si(_t_1, Delta, 24, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_mul_si(_t_2, Delta, 2, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 10, rnd);
		mpfr_mul_si(_t_2, S, 16, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, S, 36, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, P, 8, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 6, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[1][1], _t_0, 11, rnd);
		mpfr_mul_si(_t_0, Delta, 3, rnd);
		mpfr_mul_si(_t_1, S, 12, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 2, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[1][2], _t_0, 6, rnd);
		mpfr_set_si(a[1][3], 1, MPFR_RNDN);
		mpfr_mul_si(_t_0, Delta, -8, rnd);
		mpfr_mul_si(_t_1, Delta, 4, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, Delta, rnd);
		mpfr_mul(_t_0, _t_0, S, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_mul_si(_t_2, Delta, 4, rnd);
		mpfr_mul_si(_t_3, Delta, 32, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, Delta, rnd);
		mpfr_mul(_t_3, _t_3, S, rnd);
		mpfr_mul_si(_t_1, Delta, 8, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_3, _t_3, S, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_add(_t_0, _t_0, _t_3, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_mul_si(_t_2, Delta, 52, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, Delta, 8, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, P, rnd);
		mpfr_mul_si(_t_2, Delta, 12, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 16, rnd);
		mpfr_mul_si(_t_2, S, 64, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, S, 32, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, P, rnd);
		mpfr_mul_si(_t_2, S, 32, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, S, 72, rnd);
		mpfr_mul_si(_t_2, P, 24, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, P, 12, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 8, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[2][0], _t_0, 16, rnd);
		mpfr_mul_si(_t_0, Delta, -4, rnd);
		mpfr_mul_si(_t_1, Delta, 24, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_mul_si(_t_2, Delta, 4, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 20, rnd);
		mpfr_mul_si(_t_2, S, 32, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, S, 16, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, S, 60, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, P, 8, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 4, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[2][1], _t_0, 24, rnd);
		mpfr_mul_si(_t_0, Delta, -6, rnd);
		mpfr_mul_si(_t_1, S, 12, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[2][2], _t_0, 12, rnd);
		mpfr_set_si(a[2][3], -2, MPFR_RNDN);
		mpfr_mul_si(_t_0, Delta, 8, rnd);
		mpfr_mul_si(_t_1, Delta, 4, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, Delta, rnd);
		mpfr_mul(_t_0, _t_0, S, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_mul_si(_t_2, Delta, 8, rnd);
		mpfr_mul_si(_t_3, Delta, 32, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, Delta, rnd);
		mpfr_mul(_t_3, _t_3, S, rnd);
		mpfr_mul_si(_t_1, Delta, 24, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_3, _t_3, S, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_add(_t_0, _t_0, _t_3, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_mul_si(_t_2, Delta, 52, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, Delta, 8, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, P, rnd);
		mpfr_mul_si(_t_2, Delta, 16, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 42, rnd);
		mpfr_mul_si(_t_2, S, 96, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, S, 32, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, P, rnd);
		mpfr_mul_si(_t_2, S, 48, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, S, 72, rnd);
		mpfr_mul_si(_t_2, P, 24, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, P, 28, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 12, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[3][0], _t_0, 54, rnd);
		mpfr_mul_si(_t_0, Delta, -4, rnd);
		mpfr_mul_si(_t_1, Delta, 24, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_mul_si(_t_2, Delta, 4, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 32, rnd);
		mpfr_mul_si(_t_2, S, 32, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, S, rnd);
		mpfr_mul_si(_t_1, S, 16, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, S, 60, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, P, 8, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 4, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[3][1], _t_0, 54, rnd);
		mpfr_mul_si(_t_0, Delta, -6, rnd);
		mpfr_mul_si(_t_1, S, 12, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[3][2], _t_0, 18, rnd);
		mpfr_set_si(a[3][3], -2, MPFR_RNDN);
		mpfr_mul_si(_t_0, Delta, 8, rnd);
		mpfr_mul_si(_t_1, Delta, 4, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, Delta, rnd);
		mpfr_mul(_t_0, _t_0, S, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_mul_si(_t_2, Delta, 4, rnd);
		mpfr_mul_si(_t_3, Delta, 24, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, Delta, rnd);
		mpfr_mul(_t_3, _t_3, S, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_3, _t_3, epsilon, rnd);
		mpfr_mul_si(_t_1, Delta, 76, rnd);
		mpfr_add(_t_0, _t_0, _t_3, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_mul_si(_t_2, Delta, 8, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, P, rnd);
		mpfr_mul_si(_t_1, Delta, 20, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, Delta, 20, rnd);
		mpfr_mul_si(_t_2, S, 64, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, S, 144, rnd);
		mpfr_mul_si(_t_2, P, 8, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_sub(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, P, 28, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 24, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[4][0], _t_0, 24, rnd);
		mpfr_mul_si(_t_0, Delta, 2, rnd);
		mpfr_mul_si(_t_1, Delta, 24, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, S, rnd);
		mpfr_mul_si(_t_2, Delta, 6, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, Delta, 16, rnd);
		mpfr_mul_si(_t_2, S, 16, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul(_t_2, _t_2, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_2, rnd);
		mpfr_mul_si(_t_1, S, 84, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, P, 8, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 14, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[4][1], _t_0, 26, rnd);
		mpfr_mul_si(_t_0, Delta, 3, rnd);
		mpfr_mul_si(_t_1, S, 12, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 2, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[4][2], _t_0, 9, rnd);
		mpfr_set_si(a[4][3], 1, MPFR_RNDN);
		mpfr_mul_si(_t_0, Delta, 2, rnd);
		mpfr_sub_si(_t_1, Delta, 4, rnd);
		mpfr_sub_si(_t_0, _t_0, 5, rnd);
		mpfr_mul_si(_t_2, epsilon, 2, rnd);
		mpfr_mul(_t_0, _t_1, _t_0, rnd);
		mpfr_sub_si(_t_2, _t_2, 3, rnd);
		mpfr_mul(a[5][0], _t_0, _t_2, rnd);
		mpfr_mul_si(_t_0, Delta, 2, rnd);
		mpfr_mul_si(_t_1, Delta, 6, rnd);
		mpfr_mul(_t_0, _t_0, Delta, rnd);
		mpfr_mul(_t_1, _t_1, epsilon, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, Delta, 22, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_mul_si(_t_1, epsilon, 18, rnd);
		mpfr_sub(_t_0, _t_0, _t_1, rnd);
		mpfr_add_si(a[5][1], _t_0, 47, rnd);
		mpfr_mul_si(_t_0, Delta, 3, rnd);
		mpfr_mul_si(_t_1, epsilon, 2, rnd);
		mpfr_add(_t_0, _t_0, _t_1, rnd);
		mpfr_sub_si(a[5][2], _t_0, 12, rnd);
		mpfr_set_si(a[5][3], 1, MPFR_RNDN);
		mpfr_clear(_t_0);
		mpfr_clear(_t_1);
		mpfr_clear(_t_2);
		mpfr_clear(_t_3);
	}
}  // namespace qboot2
