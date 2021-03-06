#include "qboot/hor_formula.hpp"

using qboot::algebra::Vector, qboot::algebra::Polynomial;
using qboot::mp::real;

namespace qboot
{
	Vector<Polynomial> _get_zero_spin_rec_coeffs(const PrimaryOperator& op, const real& S, const real& P)
	{
		assert(op.spin() == 0);
		const auto& epsilon = op.epsilon();
		const auto& Delta = op.delta();
		Vector<real> b(4);
		Vector<Polynomial> ps(6);
		real t0, t1, t2, t3;
		// \sum_i a[0, i] n ^ i = n (n + Delta - 1) (n + 2 Delta - 2 epsilon - 2)
		ps[0] = Polynomial{1};
		ps[0]._mul_linear(Delta - 1);
		ps[0]._mul_linear(2 * Delta - (2 * epsilon + 2));
		t0 = 4 * Delta;
		t0 *= S;
		t0 += Delta;
		t2 = 8 * S;
		t1 = Delta * (-2);
		t3 = 2 * epsilon;
		t0 -= t2;
		t2 = 4 * P;
		t1 += t3;
		t0 += t2;
		t1 += 3;
		t0 -= 2;
		// a[1, 0] = (3 - 2 * Delta + 2 * epsilon) * (-2 + Delta + 4 * P - 8 * S + 4 * Delta * S);
		b[0] = t1 * t0;
		t0 = 2 * Delta;
		t1 = 24 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 2 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 10 * Delta;
		t2 = 16 * S;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 36 * S;
		t0 += t1;
		t1 = 8 * P;
		t0 -= t1;
		t1 = 6 * epsilon;
		t0 += t1;
		// a[1, 1] = 2 * Delta * (Delta - epsilon - 12 * S - 5) + (8 * S + 3) * (4 * epsilon + 9) / 2 - 8 * P - 5 / 2;
		b[1] = t0 + 11;
		t0 = 3 * Delta;
		t1 = 12 * S;
		t0 -= t1;
		t1 = 2 * epsilon;
		t0 -= t1;
		// a[1, 2] = 3 * Delta - 2 * epsilon - 12 * S - 6;
		b[2] = t0 - 6;
		// a[1, 3] = 1;
		b[3] = 1;
		ps[1] = {b[0], b[1], b[2], b[3]};
		t0 = Delta * (-8);
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 = 4 * Delta;
		t3 = 32 * Delta;
		t0 -= t1;
		t2 *= Delta;
		t3 *= S;
		t1 = 8 * Delta;
		t0 += t2;
		t3 *= S;
		t1 *= S;
		t0 += t3;
		t1 *= epsilon;
		t2 = 52 * Delta;
		t0 += t1;
		t2 *= S;
		t1 = 8 * Delta;
		t0 += t2;
		t1 *= P;
		t2 = 12 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 16 * Delta;
		t2 = 64 * S;
		t0 -= t1;
		t2 *= S;
		t1 = 32 * S;
		t0 -= t2;
		t1 *= P;
		t2 = 32 * S;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 72 * S;
		t2 = 24 * P;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 12 * P;
		t0 -= t1;
		t1 = 8 * epsilon;
		t0 -= t1;
		// a[2, 0] = 4 * ((4 - 3 * P - S * (18 - 8 * P + 16 * S) - epsilon * (2 + 6 P + 8 S))
		// + Delta (2 * P + epsilon * (3 + 2 * S) + S * (13 + 8 * S) - 4 - Delta * (epsilon + 2 * S - 1)));
		b[0] = t0 + 16;
		t0 = Delta * (-4);
		t1 = 24 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 4 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 20 * Delta;
		t2 = 32 * S;
		t0 += t1;
		t2 *= S;
		t1 = 16 * S;
		t0 += t2;
		t1 *= epsilon;
		t0 += t1;
		t1 = 60 * S;
		t0 += t1;
		t1 = 8 * P;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 += t1;
		// a[2, 1] = -4 * (6 - epsilon - 2 * P - 15 * S - 4 * S * (epsilon + 2 * S) + Delta * (epsilon + 6 * S - 5 +
		// Delta));
		b[1] = t0 - 24;
		t0 = Delta * (-6);
		t1 = 12 * S;
		t0 -= t1;
		// a[2, 2] = 6 * (2 - Delta - 2 * S);
		b[2] = t0 + 12;
		// a[2, 3] = -2;
		b[3] = -2;
		ps[2] = {b[0], b[1], b[2], b[3]};
		t0 = 8 * Delta;
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 = 8 * Delta;
		t3 = 32 * Delta;
		t0 -= t1;
		t2 *= Delta;
		t3 *= S;
		t1 = 24 * Delta;
		t0 += t2;
		t3 *= S;
		t1 *= S;
		t0 += t3;
		t1 *= epsilon;
		t2 = 52 * Delta;
		t0 += t1;
		t2 *= S;
		t1 = 8 * Delta;
		t0 -= t2;
		t1 *= P;
		t2 = 16 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 42 * Delta;
		t2 = 96 * S;
		t0 -= t1;
		t2 *= S;
		t1 = 32 * S;
		t0 -= t2;
		t1 *= P;
		t2 = 48 * S;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 72 * S;
		t2 = 24 * P;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 28 * P;
		t0 -= t1;
		t1 = 12 * epsilon;
		t0 -= t1;
		b[0] = t0 + 54;
		t0 = Delta * (-4);
		t1 = 24 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 4 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 32 * Delta;
		t2 = 32 * S;
		t0 += t1;
		t2 *= S;
		t1 = 16 * S;
		t0 += t2;
		t1 *= epsilon;
		t0 += t1;
		t1 = 60 * S;
		t0 -= t1;
		t1 = 8 * P;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 += t1;
		b[1] = t0 - 54;
		t0 = Delta * (-6);
		t1 = 12 * S;
		t0 += t1;
		b[2] = t0 + 18;
		b[3] = -2;
		ps[3] = {b[0], b[1], b[2], b[3]};
		t0 = 8 * Delta;
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 = 4 * Delta;
		t3 = 24 * Delta;
		t0 += t1;
		t2 *= Delta;
		t3 *= S;
		t0 -= t2;
		t3 *= epsilon;
		t1 = 76 * Delta;
		t0 += t3;
		t1 *= S;
		t2 = 8 * Delta;
		t0 -= t1;
		t2 *= P;
		t1 = 20 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t0 -= t1;
		t1 = 20 * Delta;
		t2 = 64 * S;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 144 * S;
		t2 = 8 * P;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 28 * P;
		t0 += t1;
		t1 = 24 * epsilon;
		t0 += t1;
		b[0] = t0 - 24;
		t0 = 2 * Delta;
		t1 = 24 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 6 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 16 * Delta;
		t2 = 16 * S;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 84 * S;
		t0 -= t1;
		t1 = 8 * P;
		t0 -= t1;
		t1 = 14 * epsilon;
		t0 -= t1;
		b[1] = t0 + 26;
		t0 = 3 * Delta;
		t1 = 12 * S;
		t0 += t1;
		t1 = 2 * epsilon;
		t0 += t1;
		b[2] = t0 - 9;
		b[3] = 1;
		ps[4] = {b[0], b[1], b[2], b[3]};
		// \sum_i a[5, i] n ^ i = (n + Delta - 4) (n + 2 Delta - 5) (n + 2 epsilon - 3)
		ps[5] = {Delta - 4, real(1)};
		ps[5]._mul_linear(2 * Delta - 5);
		ps[5]._mul_linear(real(2 * epsilon - 3));
		return ps;
	}

	Vector<Polynomial> _get_nonzero_spin_rec_coeffs(const PrimaryOperator& op, const real& S, const real& P)  // NOLINT
	{
		const auto& epsilon = op.epsilon();
		const auto& Delta = op.delta();
		uint32_t spin = op.spin();
		Vector<real> b(5);
		Vector<Polynomial> ps(8);
		real t0, t1, t2, t3, t4;
		// \sum_i a[0, i] n ^ i
		//   = -n (n + 2 Delta - 2 epsilon - 2) (n + Delta - 2 epsilon - spin - 1) (n + Delta + spin - 1)
		ps[0] = {real(-1), 1};
		ps[0]._mul_linear(Delta + (spin - 1));
		ps[0]._mul_linear(2 * Delta - 2 * epsilon - 2);
		ps[0]._mul_linear(Delta - 2 * epsilon - spin - 1);
		t0 = Delta * (-4);
		t0 *= Delta;
		t1 = 8 * Delta;
		t0 *= S;
		t2 = Delta * Delta;
		t1 *= S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 16 * Delta;
		t0 += t1;
		t2 *= S;
		t1 = 4 * Delta;
		t0 += t2;
		t1 *= P;
		t2 = 2 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 4 * S;
		t0 += t2;
		t2 = 4 * Delta;
		t1 *= spin;
		t3 = 8 * S;
		t0 += t2;
		t1 *= spin;
		t3 *= spin;
		t0 += t1;
		t3 *= epsilon;
		t1 = 16 * S;
		t0 += t3;
		t1 *= epsilon;
		t0 -= t1;
		t1 = 16 * S;
		t2 = 8 * P;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 4 * P;
		t0 += t1;
		t1 = spin * spin;
		t2 = 2 * spin;
		t0 += t1;
		t2 *= epsilon;
		t3 = Delta * (-2);
		t1 = 2 * epsilon;
		t0 += t2;
		t2 = 4 * epsilon;
		t1 = t3 + t1;
		t0 -= t2;
		t1 += 3;
		t0 -= 4;
		b[0] = t1 * t0;
		t0 = Delta * (-2);
		t1 = 40 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 6 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 += t1;
		t2 *= epsilon;
		t1 = 16 * Delta;
		t3 = 80 * Delta;
		t0 += t2;
		t1 *= Delta;
		t3 *= S;
		t0 += t1;
		t3 *= epsilon;
		t1 = 128 * Delta;
		t0 -= t3;
		t1 *= S;
		t2 = 16 * Delta;
		t3 = 2 * Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 4 * Delta;
		t0 += t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 4 * Delta;
		t0 += t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 32 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 8 * S;
		t0 -= t1;
		t1 = 38 * Delta;
		t2 *= spin;
		t3 = 16 * S;
		t0 -= t1;
		t2 *= spin;
		t3 *= spin;
		t1 = 32 * S;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= epsilon;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 128 * S;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 104 * S;
		t2 = 24 * P;
		t0 += t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 -= t2;
		t2 = 20 * P;
		t1 *= spin;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 4 * spin;
		t3 = 4 * spin;
		t0 -= t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 -= t2;
		t3 *= epsilon;
		t1 = 8 * spin;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 12 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 38 * epsilon;
		t0 += t1;
		b[1] = t0 + 28;
		t0 = Delta * (-5);
		t1 = 48 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 10 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 22 * Delta;
		t2 = 48 * S;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 72 * S;
		t0 -= t1;
		t1 = 8 * P;
		t0 += t1;
		t1 = spin * spin;
		t2 = 2 * spin;
		t0 += t1;
		t2 *= epsilon;
		t1 = 4 * epsilon;
		t0 += t2;
		t1 *= epsilon;
		t0 -= t1;
		t1 = 22 * epsilon;
		t0 -= t1;
		b[2] = t0 - 23;
		t0 = Delta * (-4);
		t1 = 16 * S;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 += t1;
		b[3] = t0 + 8;
		b[4] = -1;
		ps[1] = {b[0], b[1], b[2], b[3], b[4]};
		t0 = 8 * Delta;
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 8 * Delta;
		t3 = 64 * Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 *= Delta;
		t3 *= Delta;
		t4 = 24 * Delta;
		t0 += t1;
		t2 *= Delta;
		t3 *= S;
		t4 *= Delta;
		t0 -= t2;
		t3 *= S;
		t4 *= S;
		t1 = 84 * Delta;
		t0 -= t3;
		t4 *= epsilon;
		t1 *= Delta;
		t2 = 8 * Delta;
		t3 = 8 * Delta;
		t0 -= t4;
		t1 *= S;
		t2 *= Delta;
		t3 *= Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= epsilon;
		t1 = 12 * Delta;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= Delta;
		t2 = 96 * Delta;
		t0 -= t3;
		t1 *= epsilon;
		t3 = 60 * Delta;
		t2 *= S;
		t0 -= t1;
		t3 *= Delta;
		t2 *= S;
		t1 = 304 * Delta;
		t0 += t3;
		t2 *= epsilon;
		t1 *= S;
		t3 = 64 * Delta;
		t4 = 8 * Delta;
		t0 += t2;
		t1 *= S;
		t3 *= S;
		t4 *= S;
		t2 = 16 * Delta;
		t0 += t1;
		t3 *= P;
		t4 *= spin;
		t2 *= S;
		t1 = 16 * Delta;
		t0 -= t3;
		t4 *= spin;
		t2 *= spin;
		t1 *= S;
		t0 -= t4;
		t2 *= epsilon;
		t1 *= epsilon;
		t3 = 168 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t3 *= S;
		t0 += t1;
		t3 *= epsilon;
		t1 = 256 * Delta;
		t2 = 40 * Delta;
		t0 += t3;
		t1 *= S;
		t2 *= P;
		t3 = 4 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t1 = 12 * Delta;
		t3 *= spin;
		t0 += t2;
		t1 *= P;
		t3 *= spin;
		t2 = 8 * Delta;
		t4 = 8 * Delta;
		t0 += t1;
		t3 *= epsilon;
		t2 *= spin;
		t4 *= spin;
		t0 -= t3;
		t2 *= spin;
		t4 *= epsilon;
		t1 = 16 * Delta;
		t0 += t2;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 32 * Delta;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 4 * Delta;
		t3 = 16 * S;
		t0 += t2;
		t1 *= epsilon;
		t3 *= S;
		t2 = 32 * S;
		t0 -= t1;
		t1 = 156 * Delta;
		t3 *= spin;
		t2 *= S;
		t0 -= t1;
		t3 *= spin;
		t2 *= spin;
		t1 = 256 * S;
		t0 += t3;
		t2 *= epsilon;
		t1 *= S;
		t0 += t2;
		t1 *= epsilon;
		t2 = 352 * S;
		t3 = 96 * S;
		t0 -= t1;
		t2 *= S;
		t3 *= P;
		t1 = 8 * S;
		t0 -= t2;
		t3 *= epsilon;
		t2 = 112 * S;
		t1 *= spin;
		t0 += t3;
		t2 *= P;
		t1 *= spin;
		t3 = 20 * S;
		t4 = 16 * S;
		t0 += t2;
		t1 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 += t1;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 40 * S;
		t0 += t3;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 64 * S;
		t0 += t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 256 * S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 48 * P;
		t0 -= t1;
		t1 = 240 * S;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 40 * P;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 4 * spin;
		t0 -= t1;
		t1 = 8 * P;
		t2 *= spin;
		t0 += t1;
		t2 *= epsilon;
		t1 = 16 * spin;
		t3 = 8 * spin;
		t0 += t2;
		t1 *= spin;
		t3 *= epsilon;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 32 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 40 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 36 * epsilon;
		t0 += t1;
		b[0] = t0 + 136;
		t0 = 6 * Delta;
		t1 = 40 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 2 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 += t1;
		t2 *= epsilon;
		t1 = 58 * Delta;
		t3 = 160 * Delta;
		t0 += t2;
		t1 *= Delta;
		t3 *= S;
		t2 = 80 * Delta;
		t0 -= t1;
		t3 *= S;
		t2 *= S;
		t0 -= t3;
		t2 *= epsilon;
		t1 = 224 * Delta;
		t0 -= t2;
		t1 *= S;
		t2 = 16 * Delta;
		t3 = 6 * Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 12 * Delta;
		t0 -= t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 12 * Delta;
		t0 -= t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 12 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 128 * S;
		t0 += t1;
		t1 = 186 * Delta;
		t2 *= S;
		t0 += t1;
		t2 *= epsilon;
		t1 = 336 * S;
		t0 += t2;
		t1 *= S;
		t2 = 64 * S;
		t3 = 8 * S;
		t0 += t1;
		t2 *= P;
		t3 *= spin;
		t1 = 16 * S;
		t0 -= t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 32 * S;
		t0 -= t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 224 * S;
		t0 += t2;
		t1 *= epsilon;
		t0 += t1;
		t1 = 296 * S;
		t2 = 40 * P;
		t0 += t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 += t2;
		t2 = 12 * P;
		t1 *= spin;
		t0 += t2;
		t1 *= epsilon;
		t2 = 14 * spin;
		t3 = 4 * spin;
		t0 += t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 += t2;
		t3 *= epsilon;
		t1 = 28 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 28 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 46 * epsilon;
		t0 -= t1;
		b[1] = t0 - 194;
		t0 = 15 * Delta;
		t1 = 48 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 6 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 78 * Delta;
		t2 = 80 * S;
		t0 -= t1;
		t2 *= S;
		t1 = 48 * S;
		t0 -= t2;
		t1 *= epsilon;
		t0 -= t1;
		t1 = 120 * S;
		t0 -= t1;
		t1 = 8 * P;
		t2 = 3 * spin;
		t0 -= t1;
		t2 *= spin;
		t1 = 6 * spin;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 4 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 22 * epsilon;
		t0 += t1;
		b[2] = t0 + 107;
		t0 = 12 * Delta;
		t1 = 16 * S;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 -= t1;
		b[3] = t0 - 28;
		b[4] = 3;
		ps[2] = {b[0], b[1], b[2], b[3], b[4]};
		t0 = Delta * (-16);
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 14 * Delta;
		t3 = 64 * Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 *= Delta;
		t3 *= Delta;
		t4 = 32 * Delta;
		t0 += t1;
		t2 *= Delta;
		t3 *= S;
		t4 *= Delta;
		t0 -= t2;
		t3 *= S;
		t4 *= S;
		t1 = 208 * Delta;
		t0 -= t3;
		t4 *= epsilon;
		t1 *= Delta;
		t2 = 16 * Delta;
		t3 = 8 * Delta;
		t0 -= t4;
		t1 *= S;
		t2 *= Delta;
		t3 *= Delta;
		t0 += t1;
		t2 *= P;
		t3 *= epsilon;
		t1 = 14 * Delta;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= Delta;
		t2 = 128 * Delta;
		t0 -= t3;
		t1 *= epsilon;
		t3 = 133 * Delta;
		t2 *= S;
		t4 = 96 * Delta;
		t0 -= t1;
		t3 *= Delta;
		t2 *= S;
		t4 *= S;
		t0 += t3;
		t2 *= S;
		t4 *= S;
		t1 = 464 * Delta;
		t0 += t2;
		t4 *= epsilon;
		t1 *= S;
		t2 = 64 * Delta;
		t3 = 16 * Delta;
		t0 += t4;
		t1 *= S;
		t2 *= S;
		t3 *= S;
		t4 = 32 * Delta;
		t0 += t1;
		t2 *= P;
		t3 *= spin;
		t4 *= S;
		t1 = 64 * Delta;
		t0 += t2;
		t3 *= spin;
		t4 *= spin;
		t1 *= S;
		t0 += t3;
		t4 *= epsilon;
		t1 *= epsilon;
		t2 = 96 * Delta;
		t0 += t4;
		t1 *= epsilon;
		t2 *= S;
		t0 += t1;
		t2 *= epsilon;
		t1 = 816 * Delta;
		t3 = 16 * Delta;
		t0 += t2;
		t1 *= S;
		t3 *= P;
		t2 = 4 * Delta;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 120 * Delta;
		t2 *= spin;
		t0 -= t3;
		t1 *= P;
		t2 *= spin;
		t3 = 14 * Delta;
		t4 = 8 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 -= t2;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 28 * Delta;
		t0 += t3;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 44 * Delta;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 22 * Delta;
		t0 += t2;
		t1 *= epsilon;
		t2 = 384 * S;
		t0 -= t1;
		t1 = 432 * Delta;
		t2 *= S;
		t3 = 128 * S;
		t4 = 16 * S;
		t0 -= t1;
		t2 *= S;
		t3 *= S;
		t4 *= S;
		t1 = 32 * S;
		t0 -= t2;
		t3 *= P;
		t4 *= spin;
		t1 *= S;
		t0 += t3;
		t4 *= spin;
		t1 *= spin;
		t2 = 384 * S;
		t0 += t4;
		t1 *= epsilon;
		t2 *= S;
		t0 += t1;
		t2 *= epsilon;
		t1 = 768 * S;
		t3 = 160 * S;
		t0 -= t2;
		t1 *= S;
		t3 *= P;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 80 * S;
		t2 = 48 * S;
		t0 -= t3;
		t1 *= P;
		t2 *= spin;
		t3 = 96 * S;
		t0 -= t1;
		t2 *= spin;
		t3 *= spin;
		t1 = 192 * S;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= epsilon;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 64 * P;
		t0 -= t1;
		t1 = 1008 * S;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 48 * P;
		t0 += t2;
		t1 *= epsilon;
		t2 = 2 * spin;
		t0 -= t1;
		t1 = 216 * P;
		t2 *= spin;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 33 * spin;
		t3 = 4 * spin;
		t0 += t2;
		t1 *= spin;
		t3 *= epsilon;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 66 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 72 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 108 * epsilon;
		t0 += t1;
		b[0] = t0 + 468;
		t0 = 6 * Delta;
		t1 = 80 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 2 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 88 * Delta;
		t3 = 160 * Delta;
		t0 += t2;
		t1 *= Delta;
		t3 *= S;
		t2 = 32 * Delta;
		t0 -= t1;
		t3 *= S;
		t2 *= S;
		t0 -= t3;
		t2 *= epsilon;
		t1 = 544 * Delta;
		t0 -= t2;
		t1 *= S;
		t2 = 32 * Delta;
		t3 = 6 * Delta;
		t0 += t1;
		t2 *= P;
		t3 *= spin;
		t1 = 12 * Delta;
		t0 -= t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 12 * Delta;
		t0 -= t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 24 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 128 * S;
		t0 += t1;
		t1 = 378 * Delta;
		t2 *= S;
		t3 = 128 * S;
		t0 += t1;
		t2 *= S;
		t3 *= S;
		t0 += t2;
		t3 *= epsilon;
		t1 = 496 * S;
		t0 += t3;
		t1 *= S;
		t2 = 64 * S;
		t3 = 16 * S;
		t0 += t1;
		t2 *= P;
		t3 *= spin;
		t1 = 32 * S;
		t0 += t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 64 * S;
		t0 += t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 912 * S;
		t2 = 16 * P;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 -= t2;
		t2 = 120 * P;
		t1 *= spin;
		t0 += t2;
		t1 *= epsilon;
		t2 = 20 * spin;
		t3 = 4 * spin;
		t0 += t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 += t2;
		t3 *= epsilon;
		t1 = 40 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 36 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 102 * epsilon;
		t0 -= t1;
		b[1] = t0 - 504;
		t0 = 15 * Delta;
		t1 = 96 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 6 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 114 * Delta;
		t2 = 80 * S;
		t0 -= t1;
		t2 *= S;
		t0 -= t2;
		t1 = 288 * S;
		t0 += t1;
		t1 = 16 * P;
		t2 = 3 * spin;
		t0 -= t1;
		t2 *= spin;
		t1 = 6 * spin;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 4 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 34 * epsilon;
		t0 += t1;
		b[2] = t0 + 209;
		t0 = 12 * Delta;
		t1 = 32 * S;
		t0 -= t1;
		t1 = 4 * epsilon;
		t0 -= t1;
		b[3] = t0 - 40;
		b[4] = 3;
		ps[3] = {b[0], b[1], b[2], b[3], b[4]};
		t0 = Delta * (-16);
		t1 = 8 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 16 * Delta;
		t3 = 64 * Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 *= Delta;
		t3 *= Delta;
		t4 = 32 * Delta;
		t0 -= t1;
		t2 *= Delta;
		t3 *= S;
		t4 *= Delta;
		t0 += t2;
		t3 *= S;
		t4 *= S;
		t1 = 288 * Delta;
		t0 += t3;
		t4 *= epsilon;
		t1 *= Delta;
		t2 = 16 * Delta;
		t0 -= t4;
		t1 *= S;
		t2 *= Delta;
		t3 = 80 * Delta;
		t0 += t1;
		t2 *= P;
		t3 *= Delta;
		t1 = 128 * Delta;
		t0 += t2;
		t3 *= epsilon;
		t2 = 164 * Delta;
		t1 *= S;
		t4 = 160 * Delta;
		t0 += t3;
		t2 *= Delta;
		t1 *= S;
		t4 *= S;
		t0 -= t2;
		t1 *= S;
		t4 *= S;
		t2 = 592 * Delta;
		t0 += t1;
		t4 *= epsilon;
		t2 *= S;
		t1 = 64 * Delta;
		t3 = 16 * Delta;
		t0 += t4;
		t2 *= S;
		t1 *= S;
		t3 *= S;
		t4 = 32 * Delta;
		t0 -= t2;
		t1 *= P;
		t3 *= spin;
		t4 *= S;
		t2 = 64 * Delta;
		t0 += t1;
		t3 *= spin;
		t4 *= spin;
		t2 *= S;
		t0 += t3;
		t4 *= epsilon;
		t2 *= epsilon;
		t1 = 128 * Delta;
		t0 += t4;
		t2 *= epsilon;
		t1 *= S;
		t0 += t2;
		t1 *= epsilon;
		t2 = 1456 * Delta;
		t3 = 16 * Delta;
		t0 += t1;
		t2 *= S;
		t3 *= P;
		t1 = 8 * Delta;
		t0 -= t2;
		t3 *= epsilon;
		t2 = 104 * Delta;
		t1 *= spin;
		t0 -= t3;
		t2 *= P;
		t1 *= spin;
		t3 = 16 * Delta;
		t4 = 16 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 += t1;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 32 * Delta;
		t0 -= t3;
		t4 *= epsilon;
		t1 *= spin;
		t0 += t4;
		t1 *= epsilon;
		t2 = 280 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 512 * S;
		t0 -= t2;
		t2 = 572 * Delta;
		t1 *= S;
		t3 = 128 * S;
		t4 = 16 * S;
		t0 += t2;
		t1 *= S;
		t3 *= S;
		t4 *= S;
		t2 = 32 * S;
		t0 -= t1;
		t3 *= P;
		t4 *= spin;
		t2 *= S;
		t0 -= t3;
		t4 *= spin;
		t2 *= spin;
		t1 = 512 * S;
		t0 -= t4;
		t2 *= epsilon;
		t1 *= S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 1216 * S;
		t3 = 160 * S;
		t0 -= t1;
		t2 *= S;
		t3 *= P;
		t0 += t2;
		t3 *= epsilon;
		t1 = 368 * S;
		t2 = 64 * S;
		t0 += t3;
		t1 *= P;
		t2 *= spin;
		t3 = 128 * S;
		t0 -= t1;
		t2 *= spin;
		t3 *= spin;
		t1 = 256 * S;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= epsilon;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 64 * P;
		t0 -= t1;
		t1 = 2240 * S;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 160 * P;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 16 * spin;
		t0 += t1;
		t1 = 160 * P;
		t2 *= spin;
		t0 += t1;
		t2 *= epsilon;
		t1 = 40 * spin;
		t3 = 32 * spin;
		t0 -= t2;
		t1 *= spin;
		t3 *= epsilon;
		t0 += t1;
		t3 *= epsilon;
		t1 = 80 * spin;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 16 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 312 * epsilon;
		t0 += t1;
		b[0] = t0 - 664;
		t0 = Delta * (-6);
		t1 = 80 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 22 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 98 * Delta;
		t3 = 160 * Delta;
		t0 -= t2;
		t1 *= Delta;
		t3 *= S;
		t2 = 32 * Delta;
		t0 += t1;
		t3 *= S;
		t2 *= S;
		t0 += t3;
		t2 *= epsilon;
		t1 = 736 * Delta;
		t0 -= t2;
		t1 *= S;
		t2 = 32 * Delta;
		t3 = 6 * Delta;
		t0 += t1;
		t2 *= P;
		t3 *= spin;
		t1 = 12 * Delta;
		t0 += t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 4 * Delta;
		t0 += t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 140 * Delta;
		t0 += t2;
		t1 *= epsilon;
		t2 = 128 * S;
		t0 += t1;
		t1 = 458 * Delta;
		t2 *= S;
		t3 = 128 * S;
		t0 -= t1;
		t2 *= S;
		t3 *= S;
		t0 += t2;
		t3 *= epsilon;
		t1 = 624 * S;
		t0 += t3;
		t1 *= S;
		t2 = 64 * S;
		t3 = 16 * S;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 32 * S;
		t0 += t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 64 * S;
		t0 += t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 1584 * S;
		t2 = 16 * P;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 -= t2;
		t2 = 104 * P;
		t1 *= spin;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 22 * spin;
		t3 = 4 * spin;
		t0 += t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 -= t2;
		t3 *= epsilon;
		t1 = 44 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 20 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 214 * epsilon;
		t0 -= t1;
		b[1] = t0 + 658;
		t0 = Delta * (-15);
		t1 = 96 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 18 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 126 * Delta;
		t2 = 80 * S;
		t0 += t1;
		t2 *= S;
		t0 += t2;
		t1 = 384 * S;
		t0 += t1;
		t1 = 16 * P;
		t2 = 3 * spin;
		t0 += t1;
		t2 *= spin;
		t1 = 6 * spin;
		t0 += t2;
		t1 *= epsilon;
		t2 = 4 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 50 * epsilon;
		t0 += t1;
		b[2] = t0 - 251;
		t0 = Delta * (-12);
		t1 = 32 * S;
		t0 -= t1;
		t1 = 4 * epsilon;
		t0 -= t1;
		b[3] = t0 + 44;
		b[4] = -3;
		ps[4] = {b[0], b[1], b[2], b[3], b[4]};
		t0 = 8 * Delta;
		t1 = 8 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 22 * Delta;
		t3 = 64 * Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 *= Delta;
		t3 *= Delta;
		t4 = 56 * Delta;
		t0 -= t1;
		t2 *= Delta;
		t3 *= S;
		t4 *= Delta;
		t0 += t2;
		t3 *= S;
		t4 *= S;
		t1 = 164 * Delta;
		t0 += t3;
		t4 *= epsilon;
		t1 *= Delta;
		t2 = 8 * Delta;
		t0 += t4;
		t1 *= S;
		t2 *= Delta;
		t3 = 102 * Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= Delta;
		t1 = 160 * Delta;
		t0 += t2;
		t3 *= epsilon;
		t2 = 277 * Delta;
		t1 *= S;
		t0 += t3;
		t2 *= Delta;
		t1 *= S;
		t3 = 752 * Delta;
		t0 -= t2;
		t1 *= epsilon;
		t3 *= S;
		t2 = 64 * Delta;
		t4 = 8 * Delta;
		t0 += t1;
		t3 *= S;
		t2 *= S;
		t4 *= S;
		t1 = 16 * Delta;
		t0 -= t3;
		t2 *= P;
		t4 *= spin;
		t1 *= S;
		t3 = 48 * Delta;
		t0 -= t2;
		t4 *= spin;
		t1 *= spin;
		t3 *= S;
		t0 -= t4;
		t1 *= epsilon;
		t3 *= epsilon;
		t2 = 504 * Delta;
		t0 -= t1;
		t3 *= epsilon;
		t2 *= S;
		t0 += t3;
		t2 *= epsilon;
		t1 = 896 * Delta;
		t3 = 40 * Delta;
		t0 -= t2;
		t1 *= S;
		t3 *= P;
		t2 = 8 * Delta;
		t0 += t1;
		t3 *= epsilon;
		t1 = 100 * Delta;
		t2 *= spin;
		t0 += t3;
		t1 *= P;
		t2 *= spin;
		t3 = 22 * Delta;
		t4 = 16 * Delta;
		t0 -= t1;
		t2 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 += t2;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 44 * Delta;
		t0 -= t3;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 4 * Delta;
		t0 += t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 438 * Delta;
		t3 = 16 * S;
		t0 -= t2;
		t1 *= epsilon;
		t3 *= S;
		t2 = 32 * S;
		t0 -= t1;
		t1 = 1168 * Delta;
		t3 *= spin;
		t2 *= S;
		t0 += t1;
		t3 *= spin;
		t2 *= spin;
		t1 = 640 * S;
		t0 -= t3;
		t2 *= epsilon;
		t1 *= S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 1920 * S;
		t3 = 96 * S;
		t0 -= t1;
		t2 *= S;
		t3 *= P;
		t1 = 8 * S;
		t0 += t2;
		t3 *= epsilon;
		t2 = 336 * S;
		t1 *= spin;
		t0 -= t3;
		t2 *= P;
		t1 *= spin;
		t3 = 36 * S;
		t4 = 16 * S;
		t0 += t2;
		t1 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 -= t1;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 72 * S;
		t0 += t3;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 160 * S;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 1040 * S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 48 * P;
		t0 += t1;
		t1 = 1440 * S;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 240 * P;
		t0 += t2;
		t1 *= epsilon;
		t2 = 18 * spin;
		t0 -= t1;
		t1 = 300 * P;
		t2 *= spin;
		t0 += t1;
		t2 *= epsilon;
		t1 = 65 * spin;
		t3 = 36 * spin;
		t0 -= t2;
		t1 *= spin;
		t3 *= epsilon;
		t0 += t1;
		t3 *= epsilon;
		t1 = 130 * spin;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 40 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 580 * epsilon;
		t0 += t1;
		b[0] = t0 - 1620;
		t0 = Delta * (-6);
		t1 = 40 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 22 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 += t1;
		t2 *= epsilon;
		t1 = 128 * Delta;
		t3 = 160 * Delta;
		t0 -= t2;
		t1 *= Delta;
		t3 *= S;
		t2 = 112 * Delta;
		t0 += t1;
		t3 *= S;
		t2 *= S;
		t0 += t3;
		t2 *= epsilon;
		t1 = 416 * Delta;
		t0 += t2;
		t1 *= S;
		t2 = 16 * Delta;
		t3 = 6 * Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 12 * Delta;
		t0 += t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 4 * Delta;
		t0 += t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 176 * Delta;
		t0 += t2;
		t1 *= epsilon;
		t2 = 128 * S;
		t0 += t1;
		t1 = 746 * Delta;
		t2 *= S;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 784 * S;
		t0 += t2;
		t1 *= S;
		t2 = 64 * S;
		t3 = 8 * S;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 16 * S;
		t0 -= t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 32 * S;
		t0 -= t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 448 * S;
		t0 += t2;
		t1 *= epsilon;
		t0 -= t1;
		t1 = 968 * S;
		t2 = 40 * P;
		t0 += t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 += t2;
		t2 = 100 * P;
		t1 *= spin;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 28 * spin;
		t3 = 4 * spin;
		t0 += t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 -= t2;
		t3 *= epsilon;
		t1 = 56 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 28 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 326 * epsilon;
		t0 -= t1;
		b[1] = t0 + 1304;
		t0 = Delta * (-15);
		t1 = 48 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 18 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 162 * Delta;
		t2 = 80 * S;
		t0 += t1;
		t2 *= S;
		t1 = 48 * S;
		t0 += t2;
		t1 *= epsilon;
		t0 += t1;
		t1 = 216 * S;
		t0 -= t1;
		t1 = 8 * P;
		t2 = 3 * spin;
		t0 += t1;
		t2 *= spin;
		t1 = 6 * spin;
		t0 += t2;
		t1 *= epsilon;
		t2 = 4 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 62 * epsilon;
		t0 += t1;
		b[2] = t0 - 401;
		t0 = Delta * (-12);
		t1 = 16 * S;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 -= t1;
		b[3] = t0 + 56;
		b[4] = -3;
		ps[5] = {b[0], b[1], b[2], b[3], b[4]};
		t0 = 8 * Delta;
		t1 = 4 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 8 * Delta;
		t3 = 56 * Delta;
		t0 *= S;
		t1 *= epsilon;
		t2 *= Delta;
		t3 *= Delta;
		t0 += t1;
		t2 *= Delta;
		t3 *= S;
		t1 = 204 * Delta;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= Delta;
		t2 = 8 * Delta;
		t4 = 8 * Delta;
		t0 += t3;
		t1 *= S;
		t2 *= Delta;
		t4 *= Delta;
		t0 -= t1;
		t2 *= P;
		t4 *= epsilon;
		t1 = 68 * Delta;
		t0 -= t2;
		t4 *= epsilon;
		t1 *= Delta;
		t2 = 8 * Delta;
		t0 += t4;
		t1 *= epsilon;
		t3 = 104 * Delta;
		t2 *= S;
		t4 = 16 * Delta;
		t0 -= t1;
		t3 *= Delta;
		t2 *= spin;
		t4 *= S;
		t1 = 48 * Delta;
		t0 += t3;
		t2 *= spin;
		t4 *= spin;
		t1 *= S;
		t0 -= t2;
		t4 *= epsilon;
		t1 *= epsilon;
		t2 = 616 * Delta;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= S;
		t0 += t1;
		t2 *= epsilon;
		t1 = 1360 * Delta;
		t3 = 24 * Delta;
		t0 -= t2;
		t1 *= S;
		t3 *= P;
		t2 = 4 * Delta;
		t0 += t1;
		t3 *= epsilon;
		t1 = 92 * Delta;
		t2 *= spin;
		t0 -= t3;
		t1 *= P;
		t2 *= spin;
		t3 = 8 * Delta;
		t4 = 8 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t3 *= spin;
		t4 *= spin;
		t0 -= t2;
		t3 *= spin;
		t4 *= epsilon;
		t1 = 16 * Delta;
		t0 += t3;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 64 * Delta;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 348 * Delta;
		t3 = 8 * S;
		t0 -= t2;
		t1 *= epsilon;
		t3 *= spin;
		t0 += t1;
		t1 = 440 * Delta;
		t3 *= spin;
		t2 = 44 * S;
		t4 = 16 * S;
		t0 -= t1;
		t3 *= epsilon;
		t2 *= spin;
		t4 *= spin;
		t0 -= t3;
		t2 *= spin;
		t4 *= epsilon;
		t1 = 88 * S;
		t0 += t2;
		t4 *= epsilon;
		t1 *= spin;
		t2 = 192 * S;
		t0 -= t4;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 += t1;
		t2 *= epsilon;
		t1 = 1536 * S;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 16 * P;
		t0 += t1;
		t1 = 2640 * S;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 136 * P;
		t0 -= t2;
		t1 *= epsilon;
		t2 = 12 * spin;
		t0 += t1;
		t1 = 264 * P;
		t2 *= spin;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 24 * spin;
		t3 = 24 * spin;
		t0 += t2;
		t1 *= spin;
		t3 *= epsilon;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 48 * spin;
		t0 += t3;
		t1 *= epsilon;
		t2 = 120 * epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 540 * epsilon;
		t0 -= t1;
		b[0] = t0 + 600;
		t0 = 2 * Delta;
		t1 = 40 * Delta;
		t0 *= Delta;
		t1 *= Delta;
		t2 = 14 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 *= Delta;
		t0 += t1;
		t2 *= epsilon;
		t1 = 46 * Delta;
		t3 = 112 * Delta;
		t0 += t2;
		t1 *= Delta;
		t3 *= S;
		t0 -= t1;
		t3 *= epsilon;
		t1 = 512 * Delta;
		t0 += t3;
		t1 *= S;
		t2 = 16 * Delta;
		t3 = 2 * Delta;
		t0 -= t1;
		t2 *= P;
		t3 *= spin;
		t1 = 4 * Delta;
		t0 -= t2;
		t3 *= spin;
		t1 *= spin;
		t2 = 12 * Delta;
		t0 -= t3;
		t1 *= epsilon;
		t2 *= epsilon;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 140 * Delta;
		t0 += t2;
		t1 *= epsilon;
		t2 = 8 * S;
		t0 -= t1;
		t1 = 278 * Delta;
		t2 *= spin;
		t3 = 16 * S;
		t0 += t1;
		t2 *= spin;
		t3 *= spin;
		t1 = 32 * S;
		t0 -= t2;
		t3 *= epsilon;
		t1 *= epsilon;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 544 * S;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 1448 * S;
		t2 = 24 * P;
		t0 += t1;
		t2 *= epsilon;
		t1 = 2 * spin;
		t0 -= t2;
		t2 = 92 * P;
		t1 *= spin;
		t0 += t2;
		t1 *= epsilon;
		t2 = 10 * spin;
		t3 = 4 * spin;
		t0 -= t1;
		t2 *= spin;
		t3 *= epsilon;
		t0 += t2;
		t3 *= epsilon;
		t1 = 20 * spin;
		t0 -= t3;
		t1 *= epsilon;
		t2 = 44 * epsilon;
		t0 += t1;
		t2 *= epsilon;
		t0 -= t2;
		t1 = 318 * epsilon;
		t0 += t1;
		b[1] = t0 - 490;
		t0 = 5 * Delta;
		t1 = 48 * Delta;
		t0 *= Delta;
		t1 *= S;
		t2 = 14 * Delta;
		t0 += t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 58 * Delta;
		t2 = 48 * S;
		t0 -= t1;
		t2 *= epsilon;
		t0 += t2;
		t1 = 264 * S;
		t0 -= t1;
		t1 = 8 * P;
		t0 -= t1;
		t1 = spin * spin;
		t2 = 2 * spin;
		t0 -= t1;
		t2 *= epsilon;
		t1 = 4 * epsilon;
		t0 -= t2;
		t1 *= epsilon;
		t0 += t1;
		t1 = 62 * epsilon;
		t0 -= t1;
		b[2] = t0 + 149;
		t0 = 4 * Delta;
		t1 = 16 * S;
		t0 += t1;
		t1 = 4 * epsilon;
		t0 += t1;
		b[3] = t0 - 20;
		b[4] = 1;
		ps[6] = {b[0], b[1], b[2], b[3], b[4]};
		// \sum_i a[7, i] n ^ i
		//   = (n + 2 Delta - 7) (n + 2 epsilon - 5) (n + Delta - spin - 6) (n + Delta + 2 epsilon + spin - 6)
		ps[7] = {2 * Delta - 7, real(1)};
		ps[7]._mul_linear(real(2 * epsilon - 5));
		ps[7]._mul_linear(Delta - (spin + 6));
		ps[7]._mul_linear(Delta + (2 * epsilon + spin - 6));
		return ps;
	}
}  // namespace qboot
