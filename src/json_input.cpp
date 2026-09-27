#include "qboot/json_input.hpp"

#include <cstdint>    // for int32_t, uint32_t
#include <fstream>    // for ofstream
#include <iomanip>    // for setprecision
#include <ios>        // for defaultfloat, ios
#include <locale>     // for locale
#include <memory>     // for make_unique
#include <optional>   // for optional
#include <ostream>    // for ostream
#include <stdexcept>  // for invalid_argument, logic_error
#include <utility>    // for move

using qboot::algebra::Polynomial, qboot::algebra::Vector;
using qboot::mp::real;

namespace
{
	void write_number(std::ostream& out, const real& value)
	{
		if (value.isnan() || value.isinf()) throw std::invalid_argument("Nonfinite number in JSON input");
		// SDPB reads decimal strings directly into multiprecision numbers.
		out << '"' << value << '"';
	}

	template <class T, class Writer>
	void write_array(std::ostream& out, const T& values, Writer write)
	{
		out << '[';
		bool first = true;
		for (const auto& value : values)
		{
			if (!first) out << ',';
			first = false;
			write(out, value);
		}
		out << ']';
	}

	void write_polynomial(std::ostream& out, const Polynomial& polynomial, uint32_t size)
	{
		out << '[';
		for (uint32_t i = 0; i < size; ++i)
		{
			if (i != 0) out << ',';
			if (int32_t(i) <= polynomial.degree())
				write_number(out, polynomial[i]);
			else
				out << "\"0\"";
		}
		out << ']';
	}

	void write_basis(std::ostream& out, const qboot::PVM& constraint, uint32_t size)
	{
		out << '[';
		for (uint32_t i = 0; i < size; ++i)
		{
			if (i != 0) out << ',';
			write_polynomial(out, constraint.bilinear()[i], i + 1);
		}
		out << ']';
	}
}  // namespace

namespace qboot
{
	JSONInput::JSONInput(real&& constant, Vector<real>&& obj, uint32_t num_constraints)
	    : objectives_(obj.size() + 1),
	      constraints_(std::make_unique<std::optional<PVM>[]>(num_constraints)),
	      num_constraints_(num_constraints)
	{
		objectives_[0] = std::move(constant);
		for (uint32_t i = 0; i < obj.size(); ++i) objectives_[i + 1] = std::move(obj[i]);
		std::move(obj)._reset();
	}

	void JSONInput::register_constraint(uint32_t index, PVM&& constraint) &
	{
		if (index >= num_constraints_ || constraint.dim() == 0 || objectives_.size() != constraint.num_of_vars() + 1)
			throw std::invalid_argument("Invalid JSON constraint");
		if (constraints_[index]) throw std::logic_error("JSON constraint already registered");
		constraints_[index] = std::move(constraint);
	}

	void JSONInput::write(const fs::path& path) const
	{
		for (uint32_t i = 0; i < num_constraints_; ++i)
			if (!constraints_[i]) throw std::logic_error("Missing JSON constraint");
		std::ofstream file;
		file.exceptions(std::ios::failbit | std::ios::badbit);
		file.imbue(std::locale::classic());
		file.open(path);
		file << std::defaultfloat << std::setprecision(3 + int32_t(double(mp::global_prec) * 0.302));
		file << "{\n\"objective\":";
		write_array(file, objectives_, write_number);
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
					write_array(file, c.matrices().at(r, s),
					            [&c](std::ostream& out, const Polynomial& p)
					            { write_polynomial(out, p, c.deg() + 1); });
				}
				file << ']';
			}
			file << "],\n\"samplePoints\":";
			write_array(file, c.sample_points(), write_number);
			file << ",\n\"sampleScalings\":";
			write_array(file, c.sample_scalings(), write_number);
			file << ",\n\"bilinearBasis_0\":";
			write_basis(file, c, c.deg() / 2 + 1);
			file << ",\n\"bilinearBasis_1\":";
			write_basis(file, c, (c.deg() + 1) / 2);
			file << '}';
		}
		file << "\n]}\n";
		file.close();
	}
}  // namespace qboot
