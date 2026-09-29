#include <cstdint>     // for uint32_t
#include <cstdlib>     // for getenv
#include <exception>   // for exception
#include <filesystem>  // for path, create_directories
#include <iostream>    // for cerr
#include <locale>      // for locale, numpunct
#include <memory>      // for make_unique
#include <stdexcept>   // for runtime_error
#include <string>      // for to_string
#include <utility>     // for move

#include "mpfr.h"  // for mpfr_free_cache

#include "qboot/qboot.hpp"  // for JSONInput, PVM, PolynomialProgram, numerical types

namespace
{
	using qboot::algebra::Matrix, qboot::algebra::Polynomial, qboot::algebra::Vector;
	using qboot::mp::real;
	namespace fs = std::filesystem;

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
		// Maximize 1/8 + y, with z = 2y + 1, 3y >= 3, z <= 5. Optimum: 17/8.
		qboot::PolynomialProgram program(2);
		program.objective_constant() = real("0.125");
		program.objectives(Vector<real>{real(1), real(0)});
		program.add_equation(Vector<real>{real(-2), real(1)}, real(1));
		program.add_inequality(qboot::PolynomialInequality(2, Vector<real>{real(3), real(0)}, real(3)));
		program.add_inequality(qboot::PolynomialInequality(2, Vector<real>{real(0), real(-1)}, real(-5)));
		return program;
	}

	class TestScale : public qboot::ScaleFactor
	{
		uint32_t degree_;

	public:
		explicit TestScale(uint32_t degree) : degree_(degree) {}
		uint32_t max_degree() const override { return degree_; }
		real eval(const real& x) const override { return x + 1; }
		real sample_point(uint32_t k) const override { return real(k + 1); }
		Vector<real> sample_points() const override
		{
			Vector<real> points(degree_ + 1);
			for (uint32_t k = 0; k <= degree_; ++k) points[k] = sample_point(k);
			return points;
		}
		Vector<real> sample_scalings() const override
		{
			auto points = sample_points();
			for (auto& x : points) x = eval(x);
			return points;
		}
		Vector<Polynomial> bilinear_bases() const override
		{
			Vector<Polynomial> basis(degree_ / 2 + 1);
			for (uint32_t i = 0; i < basis.size(); ++i) basis[i] = Polynomial(i);
			return basis;
		}
	};

	qboot::PolynomialProgram consistency_program(uint32_t degree)
	{
		// Equations imply (y0, y1, y2) = (1/2 + t, 2 - t, t).
		qboot::PolynomialProgram program(3);
		program.objective_constant() = real("0.125");
		program.objectives(Vector<real>{real(3), real(5), real(7)});
		program.add_equation(Vector<real>{real(2), real(1), real(-1)}, real(3));
		program.add_equation(Vector<real>{real(0), real(1), real(1)}, real(2));
		const auto recovered = program.recover(Vector<real>{real(3)});
		if (recovered != Vector<real>{real("3.5"), real(-1), real(3)})
			throw std::runtime_error("variable recovery after two equations");
		auto scale = std::make_unique<TestScale>(degree);
		Vector<Vector<Matrix<real>>> mat(3);
		Vector<Matrix<real>> target(degree + 1);
		for (uint32_t n = 0; n < 4; ++n)
		{
			Vector<Matrix<real>> samples(degree + 1);
			for (uint32_t k = 0; k <= degree; ++k)
			{
				auto x = scale->sample_point(k);
				real value;
				for (uint32_t j = 0; j <= degree; ++j) value += (n + j + 1) * qboot::mp::pow(x, j);
				value *= scale->eval(x);
				samples[k] = Matrix<real>(2, 2);
				for (uint32_t r = 0; r < 2; ++r)
					for (uint32_t c = 0; c < 2; ++c) samples[k].at(r, c) = value * (r == c ? r + 2 : 1);
			}
			if (n < 3)
				mat[n] = std::move(samples);
			else
				target = std::move(samples);
		}
		program.add_inequality(qboot::PolynomialInequality(3, 2, std::move(scale), std::move(mat), std::move(target)));
		return program;
	}

	void test(const fs::path& output)
	{
		fs::create_directories(output);
		qboot::JSONInput input(real("0.125"), Vector<real>{real(1)}, 3);
		for (uint32_t degree = 0; degree <= 2; ++degree) input.register_constraint(degree, matrix(degree));
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
		for (uint32_t degree = 0; degree <= 2; ++degree)
		{
			auto name = "consistency-" + std::to_string(degree);
			consistency_program(degree).create_json(2).write(output / (name + ".json"));
			consistency_program(degree).create_xml().write(output / (name + ".xml"));
			consistency_program(degree).create_input(2).write(output / (name + "-sdp"), 2);
		}
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
