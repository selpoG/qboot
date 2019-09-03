#ifndef QBOOT_BOOTSTRAP_EQUATION_HPP_
#define QBOOT_BOOTSTRAP_EQUATION_HPP_

#include <cstdint>      // for uint32_t
#include <iterator>     // for begin, end
#include <map>          // for map
#include <string>       // for string
#include <type_traits>  // for is_convertible_v
#include <variant>      // for visit
#include <vector>       // for vector

#include "complex_function.hpp"    // for FunctionSymmetry
#include "context_variables.hpp"   // for Context
#include "matrix.hpp"              // for Vector
#include "polynomial_program.hpp"  // for PolynomialProgramming
#include "primary_op.hpp"          // for Operator
#include "real.hpp"                // for real

namespace qboot
{
	template <class Real>
	class Entry
	{
		uint32_t r_, c_;
		Real coeff_;
		Block<Real> block_;

	public:
		Entry(uint32_t r, uint32_t c, const Real& coeff, const Block<Real> block)
		    : r_(r), c_(c), coeff_(coeff), block_(block)
		{
		}
		Entry(const Real& coeff, const Block<Real> block) : r_(0), c_(0), coeff_(coeff), block_(block) {}
		Entry(uint32_t r, uint32_t c, const Block<Real> block) : r_(r), c_(c), coeff_(1), block_(block) {}
		explicit Entry(const Block<Real> block) : r_(0), c_(0), coeff_(1), block_(block) {}
		[[nodiscard]] const Block<Real>& block() const { return block_; }
		[[nodiscard]] uint32_t row() const { return r_; }
		[[nodiscard]] uint32_t column() const { return c_; }
		[[nodiscard]] const Real& coeff() const { return coeff_; }
	};
	// one bootrap equation contains
	// - constant terms (from the unit operator)
	// - two OPE coefficients times conformal block (from discrete or continuous spectrum)
	// and their symmetries (Odd or Even) are common (Even and Odd do not mix in one equation).
	// constant terms contributes as a normalization for linear functional alpha.
	// others form up matrices for each sector.
	// we need labels for sectors.
	// for example, in the mixed ising bootstrap, {"unit", "epsilon or sigma", "T", "even", "odd+", "odd-"}
	// (you do not have to include spin in sector name).
	// each term must tell its sector name and its position in the matrix.
	template <class Real>
	class BootstrapEquation
	{
		const Context<Real>& cont_;
		// maps from sector name to its unique id
		std::map<std::string, uint32_t> sectors_{};
		// sz_[id] is the size of matrix of id-th sector
		std::vector<uint32_t> sz_{};
		// ops_[id] is empty if operators in id-th sector are fixed
		// otherwise, ops_[id] is a list of operators
		std::vector<std::vector<GeneralPrimaryOperator<Real>>> ops_{};
		uint32_t N_ = 0, numax_;

	public:
		// seq is a sequence of tuples (operator, sector name, size)
		template <class Container>
		explicit BootstrapEquation(const Context<Real>& cont, uint32_t numax, const Container& seq)
		    : BootstrapEquation(cont, numax, std::begin(seq), std::end(seq))
		{
		}
		template <class InputIterator>
		BootstrapEquation(const Context<Real>& cont, uint32_t numax, InputIterator first, InputIterator last)
		    : cont_(cont), numax_(numax)
		{
			for (; first != last; ++first)
			{
				const auto& [sec, sz] = *first;
				static_assert(std::is_convertible_v<decltype(sec), std::string>);
				static_assert(std::is_convertible_v<decltype(sz), uint32_t>);
				auto id = uint32_t(sectors_.size());
				sectors_[sec] = id;
				sz_.push_back(uint32_t(sz));
			}
		}
		void register_operator(const std::string& sec, const GeneralPrimaryOperator<Real>& op) &
		{
			ops_[sectors_[sec]].push_back(op);
		}
		// sequence of tuples (sector name, entry)
		template <class Container>
		void add_equation(algebra::FunctionSymmetry sym, const Container& cont) &
		{
			add_equation(sym, std::begin(cont), std::end(cont));
		}
		template <class InputIterator>
		void add_equation(algebra::FunctionSymmetry sym, InputIterator first, InputIterator last) &
		{
			N_ += algebra::function_dimension(cont_.lambda(), sym);
			for (; first != last; ++first)
			{
				const auto& [sec, entry] = *first;
				static_assert(std::is_convertible_v<decltype(sec), std::string>);
				static_assert(std::is_convertible_v<decltype(entry), Entry<Real>>);
				assert(std::visit([](auto b) { return b.symmetry(); }, entry.block()) == sym);
				auto id = sectors_[sec];
				assert(entry.row() < sz_[id] && entry.column() < sz_[id]);
			}
		}
		[[nodiscard]] PolynomialProgramming<Real> create_pmp() const
		{
			PolynomialProgramming<Real> prg(N_);
			prg.objective_constant() = Real(0);
			prg.objectives(algebra::Vector(N_, Real(0)));
			// prg.add_equation
			// prg.add_inequality
			return prg;
		}
	};
}  // namespace qboot

#endif  // QBOOT_BOOTSTRAP_EQUATION_HPP_
