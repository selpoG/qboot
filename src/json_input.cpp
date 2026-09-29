#include "qboot/json_input.hpp"

#include <cassert>     // for assert
#include <cstdint>     // for int32_t, uint32_t
#include <filesystem>  // for path
#include <fstream>     // for ofstream
#include <iomanip>     // for setprecision
#include <ios>         // for defaultfloat, ios
#include <locale>      // for locale
#include <memory>      // for make_unique
#include <optional>    // for optional
#include <ostream>     // for ostream
#include <stdexcept>   // for invalid_argument
#include <utility>     // for move

namespace fs = std::filesystem;

using qboot::algebra::Vector, qboot::algebra::Polynomial;
using qboot::mp::real;
using std::make_unique, std::optional, fs::path, std::ostream, std::ofstream;

namespace
{
	void write_number(ostream& out, const real& x)
	{
		if (x.isnan() || x.isinf()) throw std::invalid_argument("Nonfinite number in JSON input");
		// SDPB reads decimal strings directly into multiprecision numbers.
		out << '"' << x << '"';
	}

	void write_vec(ostream& out, const Vector<real>& v)
	{
		out << '[';
		for (uint32_t i = 0; i < v.size(); ++i)
		{
			if (i != 0) out << ',';
			write_number(out, v[i]);
		}
		out << ']';
	}

	void write_pol(ostream& out, const Polynomial& v, uint32_t size)
	{
		out << '[';
		for (uint32_t i = 0; i < size; ++i)
		{
			if (i != 0) out << ',';
			if (int32_t(i) <= v.degree())
				write_number(out, v[i]);
			else
				out << "\"0\"";
		}
		out << ']';
	}

	void write_basis(ostream& out, const Vector<Polynomial>& v, uint32_t size)
	{
		out << '[';
		for (uint32_t i = 0; i < size; ++i)
		{
			if (i != 0) out << ',';
			write_pol(out, v[i], i + 1);
		}
		out << ']';
	}
}  // namespace

namespace qboot
{
	JSONInput::JSONInput(real&& constant, Vector<real>&& obj, uint32_t num_constraints)
	    : objectives_(obj.size() + 1),
	      constraints_(make_unique<optional<PVM>[]>(num_constraints)),
	      num_constraints_(num_constraints)
	{
		objectives_[0] = std::move(constant);
		for (uint32_t i = 0; i < obj.size(); ++i) objectives_[i + 1] = std::move(obj[i]);
		std::move(obj)._reset();
	}
	void JSONInput::register_constraint(uint32_t index, PVM&& c) &
	{
		assert(index < num_constraints_);
		assert(!constraints_[index].has_value());
		assert(c.dim() > 0 && objectives_.size() == c.num_of_vars() + 1);
		constraints_[index] = std::move(c);
	}
	void JSONInput::write(const path& path) const
	{
		for (uint32_t i = 0; i < num_constraints_; ++i) assert(constraints_[i].has_value());
		ofstream file;
		file.exceptions(std::ios::failbit | std::ios::badbit);
		file.imbue(std::locale::classic());
		file.open(path);
		file << std::defaultfloat << std::setprecision(3 + int32_t(double(mp::global_prec) * 0.302));
		file << "{\n\"objective\":";
		write_vec(file, objectives_);
		// The first coefficient is the constant term: its variable is fixed to one.
		file << ",\n\"normalization\":[\"1\"";
		for (uint32_t i = 1; i < objectives_.size(); ++i) file << ",\"0\"";
		file << "],\n\"PositiveMatrixWithPrefactorArray\":[\n";
		for (uint32_t j = 0; j < num_constraints_; ++j)
		{
			if (j != 0) file << ",\n";
			const auto& c = constraints_[j].value();
			file << "{\"polynomials\":[";
			for (uint32_t r = 0; r < c.dim(); ++r)
			{
				if (r != 0) file << ',';
				file << '[';
				for (uint32_t s = 0; s < c.dim(); ++s)
				{
					if (s != 0) file << ',';
					// Preserve the sampling degree even when leading coefficients cancel.
					file << '[';
					for (uint32_t n = 0; n <= c.num_of_vars(); ++n)
					{
						if (n != 0) file << ',';
						write_pol(file, c.matrices().at(r, s).at(n), c.deg() + 1);
					}
					file << ']';
				}
				file << ']';
			}
			file << "],\n\"samplePoints\":";
			write_vec(file, c.sample_points());
			file << ",\n\"sampleScalings\":";
			write_vec(file, c.sample_scalings());
			file << ",\n\"bilinearBasis_0\":";
			write_basis(file, c.bilinear(), c.deg() / 2 + 1);
			file << ",\n\"bilinearBasis_1\":";
			write_basis(file, c.bilinear(), (c.deg() + 1) / 2);
			file << '}';
		}
		file << "\n]}\n";
		file.close();
	}
}  // namespace qboot
