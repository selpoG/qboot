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
	// represents a semidefiniteness of a polynomial matrix
	//   M_0(x) + \sum_{n=1}^{N} y_n M_n(x) >= 0, (M_n.row() = M_n.column() = dim)
	// which is equivalent to a set of constraints
	//   Tr(A_p Y) + (B y)_p = c_p (B.row() = P, B.column() == N).
	// p <-> (r, s, k) where 0 <= r <= s < dim, 0 <= k <= deg,
	// and 0 <= p < P = schur_size().
	class DualConstraint
	{
		uint32_t dim_{}, deg_{};
		algebra::Matrix<mpfr::real> constraint_B_{};
		algebra::Vector<mpfr::real> constraint_c_{};
		std::array<algebra::Matrix<mpfr::real>, 2> bilinear_{};

	public:
		DualConstraint() = default;
		DualConstraint(uint32_t dim, uint32_t deg, algebra::Matrix<mpfr::real>&& B, algebra::Vector<mpfr::real>&& c,
		               std::array<algebra::Matrix<mpfr::real>, 2>&& bilinear);
		[[nodiscard]] uint32_t dim() const { return dim_; }
		[[nodiscard]] uint32_t degree() const { return deg_; }
		[[nodiscard]] uint32_t schur_size() const { return (deg_ + 1) * dim_ * (dim_ + 1) / 2; }
		[[nodiscard]] const std::array<algebra::Matrix<mpfr::real>, 2>& bilinear() const { return bilinear_; }
		[[nodiscard]] const algebra::Matrix<mpfr::real>& obj_B() const { return constraint_B_; }
		[[nodiscard]] const algebra::Vector<mpfr::real>& obj_c() const { return constraint_c_; }
	};

	// maximize b_0 + \sum_{n=1}^{N} b_n y_n over y_n such that all constraints are satisfied.
	// constant_term_ = b_0, objectives_[n] = b_n
	class SDPBInput
	{
		mpfr::real constant_term_;
		algebra::Vector<mpfr::real> objectives_;
		std::unique_ptr<std::optional<DualConstraint>[]> constraints_;
		uint32_t num_constraints_;
		void write_objectives(const std::filesystem::path& root) const;
		void write_objectives(std::ostream& out) const;
		void write_bilinear_bases(const std::filesystem::path& root) const;
		void write_bilinear_bases(std::ostream& out) const;
		void write_free_var_matrix(const std::filesystem::path& root, uint32_t i) const;
		void write_free_var_matrix(std::ostream& out, const DualConstraint& cons) const;
		void write_primal_objective_c(const std::filesystem::path& root, uint32_t i) const;
		void write_primal_objective_c(std::ostream& out, const DualConstraint& cons) const;
		void write_blocks(const std::filesystem::path& root) const;
		void write_blocks(std::ostream& out) const;

	public:
		SDPBInput(mpfr::real&& constant, algebra::Vector<mpfr::real>&& obj, uint32_t num_constraints);
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		// call this for index = 0, ..., num_constraints() - 1 before call of write
		void register_constraint(uint32_t index, DualConstraint&& c) &;
		void write(const std::filesystem::path& root_) const;
	};
}  // namespace qboot

#endif  // QBOOT_SDPB_INPUT_HPP_
