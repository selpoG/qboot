#ifndef QBOOT_XML_INPUT_HPP_
#define QBOOT_XML_INPUT_HPP_

#include <filesystem>  // for path
#include <memory>      // for unique_ptr
#include <optional>    // for optional

#include "qboot/algebra/matrix.hpp"      // for Vector, Matrix
#include "qboot/algebra/polynomial.hpp"  // for Polynomial
#include "qboot/mp/real.hpp"             // for real

namespace qboot
{
	class PVM
	{
		uint32_t dim_{}, deg_{}, num_of_vars_{};
		algebra::Matrix<algebra::Vector<algebra::Polynomial>> M_{};
		algebra::Vector<mp::real> sample_points_{}, sample_scalings_{};
		algebra::Vector<algebra::Polynomial> bilinear_{};

	public:
		void _reset() &&
		{
			dim_ = deg_ = num_of_vars_ = 0;
			std::move(M_)._reset();
			std::move(sample_points_)._reset();
			std::move(sample_scalings_)._reset();
			std::move(bilinear_)._reset();
		}
		PVM() = default;
		PVM(algebra::Matrix<algebra::Vector<algebra::Polynomial>>&& M, algebra::Vector<mp::real>&& x,
		    algebra::Vector<mp::real>&& scale, algebra::Vector<algebra::Polynomial>&& bilinear);
		[[nodiscard]] uint32_t dim() const { return dim_; }
		[[nodiscard]] uint32_t deg() const { return deg_; }
		[[nodiscard]] uint32_t num_of_vars() const { return num_of_vars_; }
		[[nodiscard]] const algebra::Vector<mp::real>& sample_points() const { return sample_points_; }
		[[nodiscard]] const algebra::Vector<mp::real>& sample_scalings() const { return sample_scalings_; }
		[[nodiscard]] const algebra::Matrix<algebra::Vector<algebra::Polynomial>>& matrices() const { return M_; }
		[[nodiscard]] const algebra::Vector<algebra::Polynomial>& bilinear() const { return bilinear_; }
	};
	class XMLInput
	{
		algebra::Vector<mp::real> objectives_;
		std::unique_ptr<std::optional<PVM>[]> constraints_;
		uint32_t num_constraints_;

	public:
		void _reset() &&
		{
			std::move(objectives_)._reset();
			constraints_.reset();
			num_constraints_ = 0;
		}
		XMLInput(mp::real&& constant, algebra::Vector<mp::real>&& obj, uint32_t num_constraints);
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		void register_constraint(uint32_t index, PVM&& c) &;
		void write(const std::filesystem::path& path) const;
	};
}  // namespace qboot

#endif  // QBOOT_XML_INPUT_HPP_
