#include "qboot/algebra/complex_function.hpp"

#include "qboot/mp/real.hpp"  // for real, pow

using qboot::algebra::ComplexFunction, qboot::algebra::FunctionSymmetry;
using qboot::mp::real;
using std::move;

namespace
{
	template <class Ring>
	ComplexFunction<Ring> _proj(ComplexFunction<Ring>&& f, FunctionSymmetry sym)
	{
		if (f.symmetry() == sym) return move(f);
		ComplexFunction<Ring> z(f.lambda(), sym);
		if (f.symmetry() == FunctionSymmetry::Mixed)
		{
			for (uint32_t dy = 0; dy <= f.lambda() / 2; ++dy)
				for (uint32_t dx = 0; dx + 2 * dy <= f.lambda(); ++dx)
					if (_matches(sym, dx)) z.at(dx, dy).swap(f.at(dx, dy));
		}
		if (sym == FunctionSymmetry::Mixed)
		{
			for (uint32_t dy = 0; dy <= f.lambda() / 2; ++dy)
				for (uint32_t dx = 0; dx + 2 * dy <= f.lambda(); ++dx)
					if (_matches(f.symmetry(), dx)) z.at(dx, dy).swap(f.at(dx, dy));
		}
		move(f)._reset();
		// otherwise (even to odd or odd to even), proj is vanishing
		return z;
	}
}  // namespace

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
}  // namespace qboot::algebra
