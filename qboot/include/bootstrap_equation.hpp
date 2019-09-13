#ifndef QBOOT_BOOTSTRAP_EQUATION_HPP_
#define QBOOT_BOOTSTRAP_EQUATION_HPP_

#include <cstdint>      // for uint32_t
#include <map>          // for map
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "complex_function.hpp"    // for FunctionSymmetry
#include "conformal_scale.hpp"     // for ConformalScale
#include "context_variables.hpp"   // for Context
#include "matrix.hpp"              // for Vector, Matrix
#include "polynomial_program.hpp"  // for PolynomialProgram
#include "primary_op.hpp"          // for GeneralPrimaryOperator
#include "real.hpp"                // for real

namespace qboot
{
	// an element in bootstrap equation
	// coeff ope[r] ope[c] block
	class Entry
	{
		uint32_t r_, c_;
		mpfr::real coeff_;
		Block block_;

	public:
		Entry(uint32_t r, uint32_t c, const mpfr::real& coeff, const Block& block)
		    : r_(r), c_(c), coeff_(coeff), block_(block)
		{
		}
		Entry(const mpfr::real& coeff, const Block& block) : Entry(0, 0, coeff, block) {}
		Entry(uint32_t r, uint32_t c, const Block& block) : Entry(r, c, mpfr::real(1), block) {}
		explicit Entry(const Block& block) : Entry(0, 0, mpfr::real(1), block) {}
		[[nodiscard]] const Block& block() const { return block_; }
		[[nodiscard]] uint32_t row() const { return r_; }
		[[nodiscard]] uint32_t column() const { return c_; }
		[[nodiscard]] const mpfr::real& coeff() const { return coeff_; }
	};
	class Equation;
	// key: sector name, value: size of matrices (or the number of OPEs in the sector)
	using SectorInformation = std::map<std::string, uint32_t, std::less<>>;
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
	// each section contains a discrete or continuous sprectrum.
	// for discrete spectrum,
	class BootstrapEquation
	{
		friend class Equation;
		const Context& cont_;
		// maps from sector name to its unique id
		std::map<std::string, uint32_t, std::less<>> sectors_{};
		// sz_[id] is the size of matrix of id-th sector
		std::vector<uint32_t> sz_{};
		// ops_[id] is empty if operators in id-th sector are fixed
		// otherwise, ops_[id] is a list of operators
		std::vector<std::vector<GeneralPrimaryOperator>> ops_{};
		// eqs_[i]: i-th equation
		// eqs_[i][id]: terms in id-th sector in i-th equation
		// std::vector<std::vector<std::vector<Entry>>> eqs_{};
		std::vector<Equation> eqs_{};
		// syms_[i]: symmetry of the i-th equation
		std::vector<algebra::FunctionSymmetry> syms_{};
		std::vector<uint32_t> offsets_{};
		std::vector<std::optional<algebra::Vector<mpfr::real>>> ope_{};
		uint32_t N_ = 0;

		[[nodiscard]] mpfr::real take_element(uint32_t id, const algebra::Matrix<mpfr::real>& m) const
		{
			if (ope_[id].has_value()) return m.inner_product(*ope_[id]);
			return m.at(0, 0);
		}
		[[nodiscard]] algebra::Vector<algebra::Matrix<mpfr::real>> make_disc_mat(uint32_t id) const;
		[[nodiscard]] std::unique_ptr<ConformalScale> common_scale(uint32_t id, const GeneralPrimaryOperator& op) const;
		[[nodiscard]] algebra::Vector<algebra::Vector<algebra::Matrix<mpfr::real>>> make_cont_mat(
		    uint32_t id, const GeneralPrimaryOperator& op, std::unique_ptr<ConformalScale>* ag) const;

		// alpha maximizes alpha(norm)
		// and satisfies alpha(target) = N and alpha(sec) >= 0 for each sector sec (!= target, norm)
		[[nodiscard]] PolynomialProgram ope_maximize(std::string_view target, std::string_view norm, mpfr::real&& N,
		                                             bool verbose = false) const;

	public:
		BootstrapEquation(const Context& cont, const SectorInformation& sectors) : cont_(cont)
		{
			for (const auto& [sec, sz] : sectors)
			{
				auto id = uint32_t(sectors_.size());
				sectors_[sec] = id;
				sz_.push_back(uint32_t(sz));
				ops_.emplace_back();
				ope_.emplace_back();
			}
		}
		void register_operator(std::string_view sec, const GeneralPrimaryOperator& op) &
		{
			ops_[sectors_.find(sec)->second].push_back(op);
		}
		void register_ope(std::string_view sec, algebra::Vector<mpfr::real>&& ope) &
		{
			auto id = sectors_.find(sec)->second;
			assert(ope.size() == sz_[id]);
			assert(!ope_[id].has_value());
			ope_[id] = std::move(ope);
		}
		void add_equation(Equation&& eq) { eqs_.push_back(std::move(eq)); }

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha(norm) = 1 and alpha(sec) >= 0 for each sector sec (!= norm)
		// the size of matrices in norm sector must be 1
		[[nodiscard]] PolynomialProgram find_contradiction(std::string_view norm, bool verbose = false) const;

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = 1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an upper bound on lambda, lambda ^ 2 <= -alpha(norm)
		[[nodiscard]] PolynomialProgram ope_maximize(std::string_view target, std::string_view norm,
		                                             bool verbose = false) const
		{
			return ope_maximize(target, norm, mpfr::real(1), verbose);
		}

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = -1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an lower bound on lambda, lambda ^ 2 >= alpha(norm)
		[[nodiscard]] PolynomialProgram ope_minimize(std::string_view target, std::string_view norm,
		                                             bool verbose = false) const
		{
			return ope_maximize(target, norm, mpfr::real(-1), verbose);
		}
	};
	class Equation
	{
		const BootstrapEquation& eq_;
		algebra::FunctionSymmetry sym_;
		std::vector<std::vector<Entry>> terms_{};

	public:
		Equation(const BootstrapEquation& eq, algebra::FunctionSymmetry sym)
		    : eq_(eq), sym_(sym), terms_(eq_.sectors_.size())
		{
		}
		const std::vector<Entry>& operator[](uint32_t id) const { return terms_[id]; }
		void add_term(std::string_view sec, uint32_t r, uint32_t c, const mpfr::real& coeff, const PrimaryOperator& o,
		              const PrimaryOperator& o1, const PrimaryOperator& o2, const PrimaryOperator& o3,
		              const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			assert(r < eq_.sz_[id] && c < eq_.sz_[id]);
			terms_[id].emplace_back(r, c, coeff, ConformalBlock<PrimaryOperator>(o, o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, const mpfr::real& coeff, const PrimaryOperator& o,
		              const PrimaryOperator& o1, const PrimaryOperator& o2, const PrimaryOperator& o3,
		              const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			terms_[id].emplace_back(coeff, ConformalBlock<PrimaryOperator>(o, o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, uint32_t r, uint32_t c, const PrimaryOperator& o, const PrimaryOperator& o1,
		              const PrimaryOperator& o2, const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			assert(r < eq_.sz_[id] && c < eq_.sz_[id]);
			terms_[id].emplace_back(r, c, ConformalBlock<PrimaryOperator>(o, o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, const PrimaryOperator& o, const PrimaryOperator& o1,
		              const PrimaryOperator& o2, const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			terms_[id].emplace_back(ConformalBlock<PrimaryOperator>(o, o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, uint32_t r, uint32_t c, const mpfr::real& coeff, const PrimaryOperator& o1,
		              const PrimaryOperator& o2, const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			assert(r < eq_.sz_[id] && c < eq_.sz_[id]);
			terms_[id].emplace_back(r, c, coeff, GeneralConformalBlock(o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, const mpfr::real& coeff, const PrimaryOperator& o1,
		              const PrimaryOperator& o2, const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			terms_[id].emplace_back(coeff, GeneralConformalBlock(o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, uint32_t r, uint32_t c, const PrimaryOperator& o1,
		              const PrimaryOperator& o2, const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			assert(r < eq_.sz_[id] && c < eq_.sz_[id]);
			terms_[id].emplace_back(r, c, GeneralConformalBlock(o1, o2, o3, o4, sym_));
		}
		void add_term(std::string_view sec, const PrimaryOperator& o1, const PrimaryOperator& o2,
		              const PrimaryOperator& o3, const PrimaryOperator& o4)
		{
			auto id = eq_.sectors_.find(sec)->second;
			terms_[id].emplace_back(GeneralConformalBlock(o1, o2, o3, o4, sym_));
		}
	};
}  // namespace qboot

#endif  // QBOOT_BOOTSTRAP_EQUATION_HPP_
