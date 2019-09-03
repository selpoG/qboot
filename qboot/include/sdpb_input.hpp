#ifndef QBOOT_SDPB_INPUT_HPP_
#define QBOOT_SDPB_INPUT_HPP_

#include <array>     // for array
#include <cassert>   // for assert
#include <cstdint>   // for uint32_t
#include <fstream>   // for ofstream
#include <iomanip>   // for setprecision
#include <ios>       // for fixed
#include <memory>    // for unique_ptr
#include <optional>  // for optional

#if __has_include(<filesystem>)
#include <filesystem>  // for path, create_directory
#else
#include <experimental/filesystem>
namespace std  // NOLINT
{
	// contamination of namespace std might be an undefined behavior
	namespace filesystem = experimental::filesystem;
}  // namespace std
#endif

#include "matrix.hpp"  // for Matrix, Vector
#include "real.hpp"    // for real

namespace qboot
{
	template <class OStr, class Real>
	OStr& write_mat(OStr& out, const algebra::Matrix<Real>& v)
	{
		out << v.row() << " " << v.column() << "\n";
		for (uint32_t r = 0; r < v.row(); ++r)
			for (uint32_t c = 0; c < v.column(); ++c) out << v.at(r, c) << "\n";
		return out;
	}

	template <class OStr, class Real>
	OStr& write_vec(OStr& out, const algebra::Vector<Real>& v)
	{
		out << v.size() << "\n";
		for (uint32_t i = 0; i < v.size(); ++i) out << v.at(i) << "\n";
		return out;
	}

	// represents a semidefiniteness of a polynomial matrix
	//   M_0(x) + \sum_{n=1}^{N} y_n M_n(x) >= 0, (M_n.row() = M_n.column() = dim)
	// which is equivalent to a set of constraints
	//   Tr(A_p Y) + (B y)_p = c_p (B.row() = P, B.column() == N).
	// p <-> (r, s, k) where 0 <= r <= s < dim, 0 <= k <= deg,
	// and 0 <= p < P = schur_size().
	template <class Real>
	class DualConstraint
	{
		uint32_t dim_{}, deg_{};
		algebra::Matrix<Real> constraint_B_{};
		algebra::Vector<Real> constraint_c_{};
		std::array<algebra::Matrix<Real>, 2> bilinear_{};

	public:
		DualConstraint() = default;
		DualConstraint(uint32_t dim, uint32_t deg, algebra::Matrix<Real>&& B, algebra::Vector<Real>&& c,
		               std::array<algebra::Matrix<Real>, 2>&& bilinear)
		    : dim_(dim),
		      deg_(deg),
		      constraint_B_(std::move(B)),
		      constraint_c_(std::move(c)),
		      bilinear_(std::move(bilinear))
		{
			uint32_t schur = schur_size();
			assert(constraint_B_.row() == schur);
			assert(constraint_c_.size() == schur);
			assert(bilinear_[0].row() == (deg_ / 2) + 1);
			assert(bilinear_[0].column() == deg_ + 1);
			assert(bilinear_[1].row() == (deg_ == 0 ? 1 : ((deg_ - 1) / 2) + 1));
			assert(bilinear_[1].column() == deg_ + 1);
		}
		[[nodiscard]] uint32_t dim() const { return dim_; }
		[[nodiscard]] uint32_t degree() const { return deg_; }
		[[nodiscard]] uint32_t schur_size() const { return (deg_ + 1) * dim_ * (dim_ + 1) / 2; }
		[[nodiscard]] const std::array<algebra::Matrix<Real>, 2>& bilinear() const { return bilinear_; }
		[[nodiscard]] const algebra::Matrix<Real>& obj_B() const { return constraint_B_; }
		[[nodiscard]] const algebra::Vector<Real>& obj_c() const { return constraint_c_; }
	};

	// maximize b_0 + \sum_{n=1}^{N} b_n y_n over y_n such that all constraints are satisfied.
	// constant_term_ = b_0, objectives_[n] = b_n
	template <class Real>
	class SDPBInput
	{
		Real constant_term_;
		algebra::Vector<Real> objectives_;
		std::unique_ptr<std::optional<DualConstraint<Real>>[]> constraints_;
		uint32_t num_constraints_;
		template <class OStr>
		void set_manip(OStr& out) const
		{
			out << std::defaultfloat << std::setprecision(3 + int32_t(Real::prec * 0.302));
		}

	public:
		SDPBInput(const Real& constant, algebra::Vector<Real>&& obj, uint32_t num_constraints)
		    : constant_term_(constant),
		      objectives_(std::move(obj)),
		      constraints_(std::make_unique<std::optional<DualConstraint<Real>>[]>(num_constraints)),
		      num_constraints_(num_constraints)
		{
		}
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		void register_constraint(uint32_t index, DualConstraint<Real>&& c) &
		{
			assert(!constraints_[index].has_value());
			assert(objectives_.size() == c.obj_B().column());
			constraints_[index] = std::move(c);
		}
		void write_objectives(const std::filesystem::path& root) const
		{
			std::ofstream file(root / "objectives");
			set_manip(file);
			write_objectives(file);
		}
		template <class OStr, class = std::enable_if_t<!std::is_same_v<OStr, std::filesystem::path>>>
		void write_objectives(OStr& out) const
		{
			out << constant_term_ << "\n";
			write_vec(out, objectives_);
		}
		void write_bilinear_bases(const std::filesystem::path& root) const
		{
			std::ofstream file(root / "bilinear_bases.0");
			set_manip(file);
			write_bilinear_bases(file);
		}
		template <class OStr, class = std::enable_if_t<!std::is_same_v<OStr, std::filesystem::path>>>
		void write_bilinear_bases(OStr& out) const
		{
			out << num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i)
				for (const auto& bas : constraints_[i]->bilinear()) write_mat(out, bas);
		}
		void write_free_var_matrix(const std::filesystem::path& root, uint32_t i) const
		{
			std::string filename = "free_var_matrix." + std::to_string(i);
			std::ofstream file(root / filename);
			set_manip(file);
			write_free_var_matrix(file, i);
		}
		template <class OStr, class = std::enable_if_t<!std::is_same_v<OStr, std::filesystem::path>>>
		void write_free_var_matrix(OStr& out, uint32_t i) const
		{
			write_free_var_matrix(out, *constraints_[i]);
		}
		template <class OStr>
		void write_free_var_matrix(OStr& out, const DualConstraint<Real>& cons) const
		{
			write_mat(out, cons.obj_B());
		}
		void write_primal_objective_c(const std::filesystem::path& root, uint32_t i) const
		{
			std::string filename = "primal_objective_c." + std::to_string(i);
			std::ofstream file(root / filename);
			set_manip(file);
			write_primal_objective_c(file, i);
		}
		template <class OStr, class = std::enable_if_t<!std::is_same_v<OStr, std::filesystem::path>>>
		void write_primal_objective_c(OStr& out, uint32_t i) const
		{
			write_primal_objective_c(out, *constraints_[i]);
		}
		template <class OStr>
		void write_primal_objective_c(OStr& out, const DualConstraint<Real>& cons) const
		{
			write_vec(out, cons.obj_c());
		}
		void write_blocks(const std::filesystem::path& root) const
		{
			std::ofstream file(root / "blocks.0");
			set_manip(file);
			write_blocks(file);
		}
		template <class OStr, class = std::enable_if_t<!std::is_same_v<OStr, std::filesystem::path>>>
		void write_blocks(OStr& out) const
		{
			out << 1 << "\n";
			out << num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i) out << i << "\n";
			out << num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i]->dim() << "\n";
			out << num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i]->degree() << "\n";
			out << num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i) out << constraints_[i]->schur_size() << "\n";
			out << 2 * num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i)
				for (const auto& m : constraints_[i]->bilinear()) out << m.row() * constraints_[i]->dim() << "\n";
			out << 2 * num_constraints_ << "\n";
			for (uint32_t i = 0; i < num_constraints_; ++i)
				for (const auto& m : constraints_[i]->bilinear()) out << m.column() * constraints_[i]->dim() << "\n";
		}
		void write_all(const std::filesystem::path& root_) const
		{
			// ensure root to be path to directory
			auto root = root_ / "";
			std::filesystem::create_directory(root);
			write_blocks(root);
			write_objectives(root);
			write_bilinear_bases(root);
			for (uint32_t i = 0; i < num_constraints_; ++i)
			{
				write_free_var_matrix(root, i);
				write_primal_objective_c(root, i);
			}
		}
	};
}  // namespace qboot

#endif  // QBOOT_SDPB_INPUT_HPP_
