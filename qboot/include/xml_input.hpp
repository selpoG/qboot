#ifndef QBOOT_XML_INPUT_HPP_
#define QBOOT_XML_INPUT_HPP_

#include <memory>    // for unique_ptr
#include <optional>  // for optional

#include "matrix.hpp"      // for Vector, Matrix
#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for real

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

namespace qboot
{
	class PVM
	{
		uint32_t dim_{}, deg_{}, num_of_vars_{};
		algebra::Matrix<algebra::Vector<algebra::Polynomial>> M_{};
		algebra::Vector<mpfr::real> sample_points_{}, sample_scalings_{};
		algebra::Vector<algebra::Polynomial> bilinear_{};

	public:
		PVM() = default;
		PVM(algebra::Matrix<algebra::Vector<algebra::Polynomial>>&& M, algebra::Vector<mpfr::real>&& x,
		    algebra::Vector<mpfr::real>&& scale, algebra::Vector<algebra::Polynomial>&& bilinear);
		[[nodiscard]] uint32_t dim() const { return dim_; }
		[[nodiscard]] uint32_t deg() const { return deg_; }
		[[nodiscard]] uint32_t num_of_vars() const { return num_of_vars_; }
		[[nodiscard]] const algebra::Vector<mpfr::real>& sample_points() const { return sample_points_; }
		[[nodiscard]] const algebra::Vector<mpfr::real>& sample_scalings() const { return sample_scalings_; }
		[[nodiscard]] const algebra::Matrix<algebra::Vector<algebra::Polynomial>>& matrices() const { return M_; }
		[[nodiscard]] const algebra::Vector<algebra::Polynomial>& bilinear() const { return bilinear_; }
	};
	class XMLInput
	{
		algebra::Vector<mpfr::real> objectives_;
		std::unique_ptr<std::optional<PVM>[]> constraints_;
		uint32_t num_constraints_;

	public:
		XMLInput(mpfr::real&& constant, algebra::Vector<mpfr::real>&& obj, uint32_t num_constraints);
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		void register_constraint(uint32_t index, PVM&& c) &;
		void write(const std::filesystem::path& path) const;
	};
}  // namespace qboot

#endif  // QBOOT_XML_INPUT_HPP_
