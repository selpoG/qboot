#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "qboot/qboot.hpp"

using qboot::Context, qboot::Sector, qboot::BootstrapEquation;
using qboot::mp::real, qboot::mp::rational, qboot::mp::parse;
using std::array, std::unique_ptr, std::cout, std::endl, std::map, std::move, std::vector, std::string;
namespace fs = qboot::fs;

template <class T>
using dict = map<string, T, std::less<>>;
using Op = qboot::PrimaryOperator;
using Eq = qboot::Equation;
constexpr auto ContinuousType = qboot::SectorType::Continuous;
constexpr auto Odd = qboot::algebra::FunctionSymmetry::Odd;
constexpr auto Even = qboot::algebra::FunctionSymmetry::Even;

static string name(const dict<rational>& deltas);
static BootstrapEquation create(const Context& c, const dict<rational>& deltas, uint32_t numax, uint32_t maxspin);

string name(const dict<rational>& deltas)
{
	return string("sdp-") + deltas.at("s").str('#') + "-" + deltas.at("e").str('#');
}

BootstrapEquation create(const Context& c, const dict<rational>& deltas, uint32_t numax = 20, uint32_t maxspin = 24)
{
	dict<Op> ops;
	ops.emplace("s", Op(real(deltas.at("s")), 0, c));
	ops.emplace("e", Op(real(deltas.at("e")), 0, c));
	auto ext = [&ops](auto o1, auto o2, auto o3, auto o4) {
		return array{ops.at(o1), ops.at(o2), ops.at(o3), ops.at(o4)};
	};
	Op u{c}, T{2, c};
	vector<Sector> secs{{"unit", 1, {real(1)}}, {"e", 1}, {"T", 1}};
	{
		Sector even("even", 1, ContinuousType);
		even.add_op(0, real("1.409"), real("1.4135"));
		even.add_op(0, real(6));
		even.add_op(2, real(5));
		for (uint32_t spin = 4; spin <= maxspin; spin += 2) even.add_op(spin);
		secs.push_back(move(even));
	}
	BootstrapEquation boot(c, move(secs), numax);
	Eq eq(boot, Odd);
	eq.add("unit", u, ext("s", "s", "s", "s"));
	eq.add("e", ops.at("e"), ext("s", "s", "s", "s"));
	eq.add("T", T, ext("s", "s", "s", "s"));
	eq.add("even", ext("s", "s", "s", "s"));
	boot.add_equation(move(eq));
	boot.finish();
	return boot;
}

int main(int argc, char* argv[])
{
	qboot::mp::global_prec = 1000;
	qboot::mp::global_rnd = MPFR_RNDN;
	constexpr uint32_t n_Max = 400, lambda = 16, maxspin = 28, numax = 12, parallel = 8;
	const rational dim("3");
	if (argc != 1 && argc != 3)
	{
		cout << "usage: ./program [delta_s delta_e]" << endl;
		return 1;
	}
	dict<rational> deltas;
	deltas["s"] = parse("0.5181489").value();
	deltas["e"] = parse("1.412625").value();
	if (argc == 3)
	{
		unique_ptr<char*[]> args(argv);
		deltas["s"] = parse(args[1]).value();
		deltas["e"] = parse(args[2]).value();
		args.release();
	}
	Context c(n_Max, lambda, dim, parallel);
	auto eqn = create(c, deltas, numax, maxspin);
	auto pmp = eqn.convert(qboot::FindContradiction("unit"), parallel);
	auto root = fs::current_path() / name(deltas);
	cout << root << endl;
	move(pmp).create_input(parallel).write(root, parallel);
	return 0;
}
