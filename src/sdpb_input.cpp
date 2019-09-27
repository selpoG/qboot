#include "qboot/sdpb_input.hpp"

#include <future>  // for future
#include <vector>  // for vector

#include "qboot/task_queue.hpp"  // for QBOOT_scope, TaskQueue, _event_base

namespace fs = qboot::fs;

using qboot::algebra::Vector, qboot::algebra::Matrix;
using qboot::mp::real;
using std::array, std::move, std::optional, std::make_unique, std::vector;
using std::ostream, std::ofstream, std::string, std::to_string, fs::path, fs::create_directory;

static ostream& write_mat(ostream& out, const Matrix<real>& v)
{
	out << v.row() << " " << v.column() << "\n";
	for (uint32_t r = 0; r < v.row(); ++r)
		for (uint32_t c = 0; c < v.column(); ++c) out << v.at(r, c) << "\n";
	return out;
}

static ostream& write_vec(ostream& out, const Vector<real>& v)
{
	out << v.size() << "\n";
	for (const auto& x : v) out << x << "\n";
	return out;
}

static void set_manip(ostream& out)
{
	out << std::defaultfloat << std::setprecision(3 + int32_t(double(qboot::mp::global_prec) * 0.302));
}

namespace qboot
{
	DualConstraint::DualConstraint(uint32_t dim, uint32_t deg, Matrix<real>&& B, Vector<real>&& c,
	                               array<Matrix<real>, 2>&& bilinear)
	    : dim_(dim), deg_(deg), constraint_B_(move(B)), constraint_c_(move(c)), bilinear_(move(bilinear))
	{
		[[maybe_unused]] auto schur = schur_size();
		assert(constraint_B_.row() == schur);
		assert(constraint_c_.size() == schur);
		assert(bilinear_[0].row() == (deg_ / 2) + 1);
		assert(bilinear_[0].column() == deg_ + 1);
		assert(bilinear_[1].row() == (deg_ == 0 ? 1 : ((deg_ - 1) / 2) + 1));
		assert(bilinear_[1].column() == deg_ + 1);
	}

	SDPBInput::SDPBInput(real&& constant, Vector<real>&& obj, uint32_t num_constraints)
	    : constant_term_(move(constant)),
	      objectives_(move(obj)),
	      constraints_(make_unique<optional<DualConstraint>[]>(num_constraints)),
	      num_constraints_(num_constraints)
	{
	}
	void SDPBInput::register_constraint(uint32_t index, DualConstraint&& c) &
	{
		assert(index < num_constraints_);
		assert(!constraints_[index].has_value());
		assert(objectives_.size() == c.obj_B().column());
		constraints_[index] = move(c);
	}
	void SDPBInput::write_objectives(const path& root) const
	{
		ofstream file(root / "objectives");
		set_manip(file);
		write_objectives(file);
	}
	void SDPBInput::write_objectives(ostream& out) const
	{
		out << constant_term_ << "\n";
		write_vec(out, objectives_);
	}
	void SDPBInput::write_bilinear_bases(const path& root) const
	{
		ofstream file(root / "bilinear_bases.0");
		set_manip(file);
		write_bilinear_bases(file);
	}
	void SDPBInput::write_bilinear_bases(ostream& out) const
	{
		out << num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i)
			for (const auto& bas : constraints_[i].value().bilinear()) write_mat(out, bas);
	}
	void SDPBInput::write_free_var_matrix(const path& root, uint32_t i) const
	{
		string filename = "free_var_matrix." + to_string(i);
		ofstream file(root / filename);
		set_manip(file);
		write_free_var_matrix(file, constraints_[i].value());
	}
	void SDPBInput::write_free_var_matrix(ostream& out, const DualConstraint& cons) const
	{
		write_mat(out, cons.obj_B());
	}
	void SDPBInput::write_primal_objective_c(const path& root, uint32_t i) const
	{
		string filename = "primal_objective_c." + to_string(i);
		ofstream file(root / filename);
		set_manip(file);
		write_primal_objective_c(file, constraints_[i].value());
	}
	void SDPBInput::write_primal_objective_c(ostream& out, const DualConstraint& cons) const
	{
		write_vec(out, cons.obj_c());
	}
	void SDPBInput::write_blocks(const path& root) const
	{
		ofstream file(root / "blocks.0");
		set_manip(file);
		write_blocks(file);
	}
	void SDPBInput::write_blocks(ostream& out) const
	{
		out << 1 << "\n";
		out << num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i) out << i << "\n";
		out << num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i].value().dim() << "\n";
		out << num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i].value().degree() << "\n";
		out << num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i].value().schur_size() << "\n";
		out << 2 * num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i)
			for (const auto& m : constraints_[i].value().bilinear())
				out << m.row() * constraints_[i].value().dim() << "\n";
		out << 2 * num_constraints_ << "\n";
		for (uint32_t i = 0; i < num_constraints_; ++i)
			for (const auto& m : constraints_[i].value().bilinear())
				out << m.column() * constraints_[i].value().dim() << "\n";
	}
	void SDPBInput::write(const path& root_, uint32_t parallel, _event_base* event) const
	{
		for (uint32_t i = 0; i < num_constraints_; ++i) assert(constraints_[i].has_value());
		// ensure root to be path to directory
		auto root = root_ / "";
		create_directory(root);
		TaskQueue q(parallel);
		vector<std::future<bool>> tasks;
		tasks.push_back(q.push([this, &root, event]() {
			_scoped_event scope("write_blocks", event);
			write_blocks(root);
			return true;
		}));
		tasks.push_back(q.push([this, &root, event]() {
			_scoped_event scope("write_objectives", event);
			write_objectives(root);
			return true;
		}));
		tasks.push_back(q.push([this, &root, event]() {
			_scoped_event scope("write_bilinear_bases", event);
			write_bilinear_bases(root);
			return true;
		}));
		for (uint32_t i = 0; i < num_constraints_; ++i)
			tasks.push_back(q.push([this, &root, i, event]() {
				_scoped_event scope("write constraint " + std::to_string(i), event);
				write_free_var_matrix(root, i);
				write_primal_objective_c(root, i);
				return true;
			}));
		for (auto&& x : tasks) x.get();
	}
}  // namespace qboot
