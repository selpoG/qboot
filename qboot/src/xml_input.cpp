#include "xml_input.hpp"

#include <fstream>      // for ofstream
#include <memory>       // for make_unique
#include <optional>     // for optional
#include <ostream>      // for ostream
#include <string_view>  // for string_view
#include <utility>      // for move

namespace fs = qboot::fs;

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using mp::real;
using std::move, std::make_unique, std::optional, fs::path, std::ostream, std::ofstream, std::string_view;

class ScopedTag
{
	ostream& os_;
	string_view tag_;

public:
	ScopedTag(ostream& os, string_view tag) : os_(os), tag_(tag) { os_ << "<" << tag_ << ">\n"; }
	ScopedTag(const ScopedTag&) = delete;
	ScopedTag(ScopedTag&&) = delete;
	ScopedTag& operator=(const ScopedTag&) = delete;
	ScopedTag& operator=(ScopedTag&&) = delete;
	~ScopedTag() { os_ << "</" << tag_ << ">\n"; }
};

static void write_vec(ostream& out, const Vector<real>& v)
{
	for (const auto& x : v)
	{
		ScopedTag elt(out, "elt");
		out << x;
	}
}

static void write_vec(ostream& out, const Vector<real>& v, string_view tag)
{
	ScopedTag tmp(out, tag);
	write_vec(out, v);
}

static void write_pol(ostream& out, const Polynomial& v)
{
	ScopedTag pol(out, "polynomial");
	if (v.iszero())
	{
		ScopedTag elt(out, "coeff");
		out << 0;
	}
	else
		for (const auto& x : v)
		{
			ScopedTag elt(out, "coeff");
			out << x;
		}
}

namespace qboot
{
	PVM::PVM(Matrix<Vector<Polynomial>>&& M, Vector<real>&& x, Vector<real>&& scale, Vector<Polynomial>&& bilinear)
	    : M_(move(M)), sample_points_(move(x)), sample_scalings_(move(scale)), bilinear_(move(bilinear))
	{
		dim_ = M_.row();
		assert(sample_points_.size() > 0);
		deg_ = sample_points_.size() - 1;
		assert(M_.is_square());
		assert(dim_ > 0 && M_.at(0, 0).size() > 0);
		num_of_vars_ = M_.at(0, 0).size() - 1;
		for (uint32_t r = 0; r < dim_; ++r)
			for (uint32_t c = 0; c < dim_; ++c)
			{
				assert(M_.at(r, c).size() == num_of_vars_ + 1);
				for (uint32_t i = 0; i <= num_of_vars_; ++i) assert(M_.at(r, c).at(i).degree() <= int32_t(deg_));
			}
		assert(sample_scalings_.size() == deg_ + 1);
		assert(bilinear_.size() == deg_ / 2 + 1);
		for (uint32_t m = 0; m <= deg_ / 2; ++m) assert(bilinear_[m].degree() == int32_t(m));
	}

	XMLInput::XMLInput(real&& constant, Vector<real>&& obj, uint32_t num_constraints)
	    : objectives_(obj.size() + 1),
	      constraints_(make_unique<optional<PVM>[]>(num_constraints)),
	      num_constraints_(num_constraints)
	{
		objectives_[0] = move(constant);
		for (uint32_t i = 0; i < obj.size(); ++i) objectives_[i + 1] = move(obj[i]);
	}
	void XMLInput::register_constraint(uint32_t index, PVM&& c) &
	{
		assert(index < num_constraints_);
		assert(!constraints_[index].has_value());
		assert(objectives_.size() == c.num_of_vars() + 1);
		constraints_[index] = move(c);
	}
	void XMLInput::write(const path& path) const
	{
		for (uint32_t i = 0; i < num_constraints_; ++i) assert(constraints_[i].has_value());
		ofstream file(path);
		file << std::defaultfloat << std::setprecision(3 + int32_t(double(mp::global_prec) * 0.302));
		ScopedTag sdp(file, "sdp");
		write_vec(file, objectives_, "objective");
		ScopedTag pvms(file, "polynomialVectorMatrices");
		for (uint32_t j = 0; j < num_constraints_; ++j)
		{
			const auto& c = *constraints_[j];
			ScopedTag pvm(file, "polynomialVectorMatrix");
			{
				ScopedTag tmp(file, "rows");
				file << c.dim();
			}
			{
				ScopedTag tmp(file, "cols");
				file << c.dim();
			}
			{
				ScopedTag elms(file, "elements");
				for (uint32_t r = 0; r < c.dim(); ++r)
					for (uint32_t s = 0; s < c.dim(); ++s)
					{
						ScopedTag pv(file, "polynomialVector");
						for (uint32_t n = 0; n <= c.num_of_vars(); ++n) write_pol(file, c.matrices().at(r, s).at(n));
					}
			}
			write_vec(file, c.sample_points(), "samplePoints");
			write_vec(file, c.sample_scalings(), "sampleScalings");
			{
				ScopedTag bil(file, "bilinearBasis");
				for (const auto& q : c.bilinear()) write_pol(file, q);
			}
		}
	}
}  // namespace qboot
