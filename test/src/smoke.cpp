#include <array>      // for array
#include <cstdint>    // for uint32_t
#include <fstream>    // for ifstream
#include <iostream>   // for cerr
#include <iterator>   // for istreambuf_iterator
#include <stdexcept>  // for runtime_error
#include <string>     // for string
#include <utility>    // for move
#include <vector>     // for vector

#include "qboot/qboot.hpp"  // for numerical types, bootstrap equations, SDPB output

namespace
{
	using qboot::mp::integer, qboot::mp::rational, qboot::mp::real;
	namespace fs = qboot::fs;

	void require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	void arithmetic()
	{
		require(integer("12345678901234567890") * 9 == integer("111111110111111111010"), "integer product");
		require(rational("1/3") + rational("1/6") == rational("1/2"), "rational sum");
		require(qboot::mp::parse("1.25e-2").value() == rational("1/80"), "decimal parsing");
		const auto root = qboot::mp::sqrt(real(2));
		require(qboot::mp::abs(root * root - 2) < real("1e-60"), "multiprecision square root");
		qboot::algebra::Polynomial p{real(1), real(-2), real(1)};
		require(p.eval(real(3)) == 4, "polynomial evaluation");
		require(p.derivative().eval(real(3)) == 4, "polynomial derivative");
		require(qboot::algebra::mul(p, p).eval(real(3)) == 16, "polynomial product");
	}

	void bootstrap(const fs::path& output, uint32_t parallel)
	{
		qboot::Context context(30, 3, rational(3), parallel);
		const qboot::PrimaryOperator sigma(real("0.518"), 0, context);
		const std::array external{sigma, sigma, sigma, sigma};
		std::vector<qboot::Sector> sectors{{"unit", 1, {real(1)}}};
		qboot::Sector even("even", 1, qboot::SectorType::Continuous);
		even.add_op(0, real("1.4"), real("1.5"));
		even.add_op(0, real(3));
		sectors.push_back(std::move(even));
		qboot::BootstrapEquation equations(context, std::move(sectors), 2);
		qboot::Equation equation(equations, qboot::algebra::FunctionSymmetry::Odd);
		equation.add("unit", qboot::PrimaryOperator(context), external);
		equation.add("even", external);
		equations.add_equation(std::move(equation));
		equations.finish();
		auto program = equations.convert(qboot::FindContradiction("unit"), parallel);
		auto input = std::move(program).create_input(parallel);
		require(input.num_constraints() == 2, "finite and infinite spectral constraints");
		std::move(input).write(output, parallel);
		auto xml_program = equations.convert(qboot::FindContradiction("unit"), parallel);
		auto xml = std::move(xml_program).create_xml(parallel);
		require(xml.num_constraints() == 2, "XML spectral constraints");
		xml.write(output / "input.xml");
	}

	std::string read(const fs::path& path)
	{
		std::ifstream input(path);
		require(input.good(), "output file missing");
		return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
	}

	void output_test()
	{
		const fs::path serial("smoke-serial"), parallel("smoke-parallel");
		fs::remove_all(serial);
		fs::remove_all(parallel);
		bootstrap(serial, 1);
		bootstrap(parallel, 2);
		require(fs::exists(serial / "blocks.0"), "block metadata missing");
		require(fs::exists(serial / "objectives"), "objectives missing");
		uint32_t count = 0;
		for (const auto& entry : fs::directory_iterator(serial))
		{
			const auto contents = read(entry.path());
			require(!contents.empty(), "empty SDPB file");
			require(contents == read(parallel / entry.path().filename()), "parallel output differs");
			require(contents.find("nan") == std::string::npos, "NaN in SDPB output");
			require(contents.find("inf") == std::string::npos, "infinity in SDPB output");
			++count;
		}
		require(count == 8, "unexpected SDPB output file count");
		fs::remove_all(serial);
		fs::remove_all(parallel);
	}
}  // namespace

int main()
{
	qboot::mp::global_prec = 256;
	qboot::mp::global_rnd = MPFR_RNDN;
	try
	{
		arithmetic();
		output_test();
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
