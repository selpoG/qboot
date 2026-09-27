#include <cstdint>    // for uint32_t
#include <cstdlib>    // for getenv
#include <exception>  // for exception
#include <iostream>   // for cerr
#include <locale>     // for locale, numpunct
#include <stdexcept>  // for runtime_error
#include <utility>    // for move

#include "mpfr.h"  // for mpfr_free_cache

#include "qboot/qboot.hpp"  // for JSONInput, PVM, PolynomialProgram, numerical types

namespace
{
	using qboot::algebra::Matrix, qboot::algebra::Polynomial, qboot::algebra::Vector;
	using qboot::mp::real;
	namespace fs = qboot::fs;

	template <class Exception, class Function>
	void require_throws(Function action)
	{
		try
		{
			action();
		}
		catch (const Exception&)
		{
			return;
		}
		throw std::runtime_error("Expected exception was not thrown");
	}

	// Exercise locale independence without depending on installed OS locales.
	class DecimalComma : public std::numpunct<char>
	{
		char do_decimal_point() const override { return ','; }
	};

	qboot::PVM matrix(uint32_t degree)
	{
		Matrix<Vector<Polynomial>> polynomials(2, 2);
		for (uint32_t r = 0; r < 2; ++r)
			for (uint32_t c = 0; c < 2; ++c)
			{
				polynomials.at(r, c) = Vector<Polynomial>(2);
				if (r == c)
				{
					// Leading zeros deliberately test preservation of the sampling degree.
					polynomials.at(r, c)[0] = Polynomial(real("1.23456789012345678901234567890123456789"));
					polynomials.at(r, c)[1] = Polynomial(real(-1));
				}
			}
		Vector<real> points(degree + 1), scalings(degree + 1);
		for (uint32_t i = 0; i <= degree; ++i)
		{
			points[i] = real(i + 1);
			scalings[i] = real(1);
		}
		Vector<Polynomial> basis(degree / 2 + 1);
		basis[0] = Polynomial(real(1));
		if (degree >= 2) basis[1] = Polynomial{real(0), real(1)};
		return {std::move(polynomials), std::move(points), std::move(scalings), std::move(basis)};
	}

	qboot::PolynomialProgram bounded_program()
	{
		// Maximize 1/8 + y, with z = 2y + 1, y >= 0, z <= 5. Optimum: 17/8.
		qboot::PolynomialProgram program(2);
		program.objective_constant() = real("0.125");
		program.objectives(Vector<real>{real(1), real(0)});
		program.add_equation(Vector<real>{real(-2), real(1)}, real(1));
		program.add_inequality(qboot::PolynomialInequality(2, Vector<real>{real(1), real(0)}, real(0)));
		program.add_inequality(qboot::PolynomialInequality(2, Vector<real>{real(0), real(-1)}, real(-5)));
		return program;
	}

	void test(const fs::path& output)
	{
		fs::create_directories(output);
		qboot::JSONInput input(real("0.125"), Vector<real>{real(1)}, 3);
		require_throws<std::logic_error>([&] { input.write(output / "missing.json"); });
		require_throws<std::invalid_argument>([&] { input.register_constraint(3, matrix(0)); });
		for (uint32_t degree = 0; degree <= 2; ++degree) input.register_constraint(degree, matrix(degree));
		require_throws<std::logic_error>([&] { input.register_constraint(0, matrix(0)); });
		const auto previous = std::locale::global(std::locale(std::locale::classic(), new DecimalComma));
		input.write(output / "matrices.json");
		std::locale::global(previous);
		require_throws<std::ios_base::failure>([&] { input.write(output / "absent" / "input.json"); });
		qboot::JSONInput nonfinite(qboot::mp::nan(), Vector<real>{real(1)}, 1);
		nonfinite.register_constraint(0, matrix(0));
		require_throws<std::invalid_argument>([&] { nonfinite.write(output / "nonfinite.json"); });
		bounded_program().create_json(1).write(output / "bounded.json");
		bounded_program().create_json(2).write(output / "bounded-parallel.json");
		bounded_program().create_xml().write(output / "bounded.xml");
	}
}  // namespace

int main()
{
	qboot::mp::global_prec = 256;
	qboot::mp::global_rnd = MPFR_RNDN;
	try
	{
		const auto* output = std::getenv("QBOOT_TEST_OUTPUT");
		test(output != nullptr ? output : "json-output");
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		mpfr_free_cache();
		return 1;
	}
	mpfr_free_cache();
	return 0;
}
