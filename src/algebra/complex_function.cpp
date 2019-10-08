#include "qboot/algebra/complex_function.hpp"

#include "qboot/mp/real.hpp"  // for real, pow

using qboot::algebra::ComplexFunction, qboot::algebra::FunctionSymmetry;
using qboot::mp::real;
using std::move;

template <class Ring>
ComplexFunction<Ring> _proj(ComplexFunction<Ring>&& f, FunctionSymmetry sym)
{
	if (f.symmetry() == sym) return move(f);
	ComplexFunction<Ring> z(f.lambda(), sym);
	if (f.symmetry() == FunctionSymmetry::Mixed)
	{
		for (uint32_t dy = 0; dy <= f.lambda() / 2; ++dy)
			for (uint32_t dx = 0; dx + 2 * dy <= f.lambda(); ++dx)
				if (matches(sym, dx)) z.at(dx, dy).swap(f.at(dx, dy));
	}
	if (sym == FunctionSymmetry::Mixed)
	{
		for (uint32_t dy = 0; dy <= f.lambda() / 2; ++dy)
			for (uint32_t dx = 0; dx + 2 * dy <= f.lambda(); ++dx)
				if (matches(f.symmetry(), dx)) z.at(dx, dy).swap(f.at(dx, dy));
	}
	move(f)._reset();
	// otherwise (even to odd or odd to even), proj is vanishing
	return z;
}

namespace qboot::algebra
{
	template class ComplexFunction<real>;
	template class ComplexFunction<Polynomial>;

	template <>
	ComplexFunction<mp::real> ComplexFunction<mp::real>::proj(FunctionSymmetry sym) &&
	{
		return _proj(move(*this), sym);
	}
	template <>
	ComplexFunction<Polynomial> ComplexFunction<Polynomial>::proj(FunctionSymmetry sym) &&
	{
		return _proj(move(*this), sym);
	}

	ComplexFunction<real> v_to_d(const real& d, uint32_t lambda)
	{
		// f.at(m, n) = (der x) ^ m (der y) ^ n ((x - 1) ^ 2 - y) ^ d / (n! m!)
		//            = (-1) ^ {n + m} 2 ^ {2 n + m} lf(d, n) lf(2 (d - n), m) / (n! m! 4 ^ d)
		// where lf(x, n) = x (x - 1) ... (x - (n - 1)) (falling factorial)
		ComplexFunction<real> f(lambda);
		f.at(0, 0) = mp::pow(4, -d);
		for (uint32_t n = 0; n <= lambda / 2; ++n)
		{
			if (n > 0) f.at(0u, n) = f.at(0u, n - 1) * 4 * (-d + (n - 1)) / n;
			for (uint32_t m = 1; m + 2 * n <= lambda; ++m)
				f.at(m, n) = f.at(m - 1, n) * (-4 * d + 2 * (m + 2 * n - 1)) / m;
		}
		return f;
	}
}  // namespace qboot::algebra
