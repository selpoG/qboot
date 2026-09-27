#ifndef QBOOT_JSON_INPUT_HPP_
#define QBOOT_JSON_INPUT_HPP_

#include <cstdint>   // for uint32_t
#include <memory>    // for unique_ptr
#include <optional>  // for optional
#include <utility>   // for move

#include "qboot/algebra/matrix.hpp"  // for Vector
#include "qboot/mp/real.hpp"         // for real
#include "qboot/my_filesystem.hpp"   // for path
#include "qboot/xml_input.hpp"       // for PVM

namespace qboot
{
	// Polynomial matrix program for SDPB 3.1+ pmp2sdp.
	class JSONInput
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
		JSONInput(mp::real&& constant, algebra::Vector<mp::real>&& obj, uint32_t num_constraints);
		[[nodiscard]] uint32_t num_constraints() const { return num_constraints_; }
		void register_constraint(uint32_t index, PVM&& c) &;
		void write(const fs::path& path) const;
	};
}  // namespace qboot

#endif  // QBOOT_JSON_INPUT_HPP_
