#include "qboot/algebra/matrix.hpp"

#include "qboot/algebra/complex_function.hpp"  // for ComplexFunction
#include "qboot/algebra/polynomial.hpp"        // for Polynomial
#include "qboot/algebra/real_function.hpp"     // for ComplexFunction

using qboot::mp::real;
using std::move;

namespace qboot::algebra
{
	// these instantiation causes compilation error on MSVC
	// template class Vector<Vector<real>>;
	// template class Vector<Matrix<real>>;
	[[maybe_unused]] static void _instantiate()
	{
		Vector<real> _;
		Matrix<real> _0;
		// but you can use them... why?
		Vector<Vector<real>> _1;
		Vector<Matrix<real>> _2;
		Vector<Polynomial> _3;
		Matrix<Polynomial> _4;
		Vector<RealFunction<real>> _5;
		Matrix<RealFunction<real>> _6;
		Vector<ComplexFunction<real>> _7;
		Matrix<ComplexFunction<real>> _8;
		Vector<RealFunction<Polynomial>> _9;
		Vector<ComplexFunction<Polynomial>> _10;
		Matrix<RealFunction<Polynomial>> _11;
		Matrix<ComplexFunction<Polynomial>> _12;
	}

	static void _add_row(Matrix<real>* mat, uint32_t f, uint32_t t, const real& x, uint32_t c0 = 0)
	{
		for (uint32_t c = c0; c < mat->column(); ++c) mat->at(t, c) += x * mat->at(f, c);
	}
	static void _swap_row(Matrix<real>* mat, uint32_t r1, uint32_t r2, uint32_t c0 = 0)
	{
		for (uint32_t c = c0; c < mat->column(); ++c) mat->at(r1, c).swap(mat->at(r2, c));
	}
	static void _multiply_row(Matrix<real>* mat, uint32_t r, const real& x, uint32_t c0 = 0)
	{
		for (uint32_t c = c0; c < mat->column(); ++c) mat->at(r, c) *= x;
	}
	real determinant(Matrix<real>&& mat)
	{
		assert(mat.is_square());
		if (mat.row() == 0) return real(1);
		if (mat.row() == 1) return mat.at(0, 0);
		if (mat.row() == 2) return mat.at(0, 0) * mat.at(1, 1) - mat.at(0, 1) * mat.at(1, 0);
		real det(1);
		for (uint32_t j = 0; j < mat.row(); ++j)
		{
			uint32_t p = j;
			for (uint32_t i = j + 1; i < mat.row(); ++i)
				if (mp::cmpabs(mat.at(p, j), mat.at(i, j)) < 0) p = i;
			if (p != j)
			{
				det.negate();
				_swap_row(&mat, p, j, j);
			}
			det *= mat.at(j, j);
			if (det.iszero()) break;
			_multiply_row(&mat, j, 1 / mat.at(j, j), j);
			for (uint32_t i = j + 1; i < mat.row(); ++i) _add_row(&mat, j, i, -mat.at(i, j), j);
		}
		move(mat)._reset();
		return det;
	}
	Matrix<real> inverse(Matrix<real>&& mat)
	{
		assert(mat.is_square());
		auto inv = Matrix<real>::constant(real(1), mat.row());
		for (uint32_t j = 0; j < mat.row(); ++j)
		{
			uint32_t p = j;
			for (uint32_t i = j + 1; i < mat.row(); ++i)
				if (mp::cmpabs(mat.at(p, j), mat.at(i, j)) < 0) p = i;
			_swap_row(&inv, p, j);
			_swap_row(&mat, p, j, j);
			auto t = 1 / mat.at(j, j);
			_multiply_row(&inv, j, t);
			_multiply_row(&mat, j, t, j);
			for (uint32_t i = j + 1; i < mat.row(); ++i)
			{
				_add_row(&inv, j, i, -mat.at(i, j));
				_add_row(&mat, j, i, -mat.at(i, j), j);
			}
		}
		for (uint32_t j = mat.row() - 1; j < mat.row(); --j)
			for (uint32_t r = j - 1; r < j; --r) _add_row(&inv, j, r, -mat.at(r, j));
		move(mat)._reset();
		return inv;
	}
	Matrix<real> cholesky_decomposition(const Matrix<real>& mat)
	{
		assert(mat.is_square());
		Matrix<real> L(mat.row(), mat.row());
		real s;
		for (uint32_t i = 0; i < mat.row(); ++i)
		{
			for (uint32_t j = 0; j <= i; ++j)
			{
				s = {};
				for (uint32_t k = 0; k < j; ++k) s += L.at(i, k) * L.at(j, k);
				s = mat.at(i, j) - s;
				L.at(i, j) = i == j ? mp::sqrt(s) : s / L.at(j, j);
			}
		}
		return L;
	}
	Matrix<real> lower_triangular_inverse(const Matrix<real>& mat)
	{
		// this must be lower triangular, i.e., at(i, j) = 0 for i < j
		assert(mat.is_square());
		Matrix<real> res(mat.row(), mat.row());
		real s;
		for (uint32_t i = 0; i < mat.row(); ++i)
		{
			res.at(i, i) = 1 / mat.at(i, i);
			for (uint32_t j = 0; j < i; ++j)
			{
				s = {};
				for (uint32_t k = j; k < i; ++k) s += mat.at(i, k) * res.at(k, j);
				res.at(i, j) = -s / mat.at(i, i);
			}
		}
		return res;
	}
}  // namespace qboot::algebra
