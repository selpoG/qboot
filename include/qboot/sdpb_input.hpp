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

#include "qboot/algebra/matrix.hpp"  // for Matrix, Vector
#include "qboot/mp/real.hpp"         // for real
#include "qboot/my_filesystem.hpp"   // for path, create_directory
#include "qboot/task_queue.hpp"      // for _event_base

namespace qboot
{
	// represents a semidefiniteness of a polynomial matrix
	//   M_0(x) + \sum_{n=1}^{N} y_n M_n(x) >= 0, (M_n.row() = M_n.column() = dim)
	// which is equivalent to a set of constraints
	//   Tr(A_p Y) + (B y)_p = c_p (B.row() = P, B.column() == N).
	// p <-> (r, s, k) where 0 <= r <= s < dim, 0 <= k <= deg,
	// and 0 <= p < P = schur_size().
	class DualConstraint
	{
		uint32_t dim_{}, deg_{};
		algebra::Matrix<mp::real> constraint_B_{};
		algebra::Vector<mp::real> constraint_c_{};
		std::array<algebra::Matrix<mp::real>, 2> bilinear_{};

	public:
		void _reset() &&
		{
			dim_ = deg_ = 0;
			std::move(constraint_B_)._reset();
			std::move(constraint_c_)._reset();
			for (auto&& x : bilinear_) std::move(x)._reset();
		}
		DualConstraint() = default;
		DualConstraint(uint32_t dim, uint32_t deg, algebra::Matrix<mp::real>&& B, algebra::Vector<mp::real>&& c,
		               std::array<algebra::Matrix<mp::real>, 2>&& bilinear);
		[[nodiscard]] uint32_t dim() const { return dim_; }
		[[nodiscard]] uint32_t degree() const { return deg_; }
		[[nodiscard]] uint32_t schur_size() const { return (deg_ + 1) * dim_ * (dim_ + 1) / 2; }
		[[nodiscard]] const std::array<algebra::Matrix<mp::real>, 2>& bilinear() const { return bilinear_; }
		[[nodiscard]] uint32_t _total_memory() const
		{
			uint32_t sum = constraint_B_.row() * constraint_B_.column() + constraint_c_.size();
			for (const auto& q : bilinear_) sum += q.row() * q.column();
			return sum;
		}
		[[nodiscard]] const algebra::Matrix<mp::real>& obj_B() const { return constraint_B_; }
		[[nodiscard]] const algebra::Vector<mp::real>& obj_c() const { return constraint_c_; }
	};

	// maximize b_0 + \sum_{n=1}^{N} b_n y_n over y_n such that all constraints are satisfied.
	// constant_term_ = b_0, objectives_[n] = b_n
	class SDPBInput
	{
		mp::real constant_term_;
		algebra::Vector<mp::real> objectives_;
		std::unique_ptr<std::optional<DualConstraint>[]> constraints_;
		uint32_t num_constraints_;
		void write_objectives(const fs::path& root) const;
		void write_objectives(std::ostream& out) const;
		void write_bilinear_bases(const fs::path& root) const;
		void write_bilinear_bases(std::ostream& out) const;
		void write_free_var_matrix(const fs::path& root, uint32_t i) const;
		void write_free_var_matrix(std::ostream& out, const DualConstraint& cons) const;
		void write_primal_objective_c(const fs::path& root, uint32_t i) const;
		void write_primal_objective_c(std::ostream& out, const DualConstraint& cons) const;
		void write_blocks(const fs::path& root) const;
		void write_blocks(std::ostream& out) const;

	public:
		void _reset() &&
		{
			std::move(constant_term_)._reset();
			std::move(objectives_)._reset();
			constraints_.reset();
			num_constraints_ = 0;
		}
		SDPBInput(mp::real&& constant, algebra::Vector<mp::real>&& obj, uint32_t num_constraints);
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		[[nodiscard]] uint32_t _total_memory() const
		{
			uint32_t sum = 1 + objectives_.size();
			for (uint32_t i = 0; i < num_constraints_; ++i)
				if (constraints_[i].has_value()) sum += constraints_[i].value()._total_memory();
			return sum;
		}
		// call this for index = 0, ..., num_constraints() - 1 before call of write
		void register_constraint(uint32_t index, DualConstraint&& c) &;
		void write(const fs::path& root_, uint32_t parallel = 1, _event_base* event = nullptr) const&;
		void write(const fs::path& root_, uint32_t parallel = 1, _event_base* event = nullptr) &&;
	};
}  // namespace qboot

#endif  // QBOOT_SDPB_INPUT_HPP_
