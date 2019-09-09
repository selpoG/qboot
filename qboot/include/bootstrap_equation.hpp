#ifndef QBOOT_BOOTSTRAP_EQUATION_HPP_
#define QBOOT_BOOTSTRAP_EQUATION_HPP_

#include <cstdint>  // for uint32_t
#include <functional>
#include <iterator>     // for begin, end
#include <map>          // for map
#include <string>       // for string
#include <type_traits>  // for is_convertible_v
#include <variant>      // for visit
#include <vector>       // for vector

#include "complex_function.hpp"    // for FunctionSymmetry
#include "context_variables.hpp"   // for Context
#include "matrix.hpp"              // for Vector
#include "pole_data.hpp"           // for ConformalScale
#include "polynomial_program.hpp"  // for PolynomialProgram
#include "primary_op.hpp"          // for Operator
#include "real.hpp"                // for real

namespace qboot
{
	// an element in bootstrap equation
	// coeff ope[r] ope[c] block
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
		Entry(const Real& coeff, const Block<Real> block) : Entry(0, 0, coeff, block) {}
		Entry(uint32_t r, uint32_t c, const Block<Real> block) : Entry(r, c, Real(1), block) {}
		explicit Entry(const Block<Real> block) : Entry(0, 0, Real(1), block) {}
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
	// each section contains a discrete or continuous sprectrum.
	// for discrete spectrum,
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
		// eqs_[i]: i-th equation
		// eqs_[i][id]: terms in id-th sector in i-th equation
		std::vector<std::vector<std::vector<Entry<Real>>>> eqs_{};
		// syms_[i]: symmetry of the i-th equation
		std::vector<algebra::FunctionSymmetry> syms_{};
		std::vector<uint32_t> offsets_{};
		std::vector<std::optional<algebra::Vector<Real>>> ope_{};
		uint32_t N_ = 0, numax_;
		[[nodiscard]] Real take_element(uint32_t id, const algebra::Matrix<Real>& m) const
		{
			if (ope_[id].has_value()) return m.inner_product(*ope_[id]);
			return m.at(0, 0);
		}
		[[nodiscard]] algebra::Vector<algebra::Matrix<Real>> make_disc_mat(uint32_t id) const
		{
			auto sz = sz_[id];
			algebra::Vector<algebra::Matrix<Real>> mat(N_);
			for (uint32_t i = 0; i < N_; ++i) mat[i] = {sz, sz};
			for (uint32_t i = 0; i < eqs_.size(); ++i)
			{
				if (eqs_[i][id].empty()) continue;
				auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
				algebra::Matrix<algebra::Vector<Real>> tmp(sz, sz);
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = algebra::Vector<Real>{n};
				for (const auto& term : eqs_[i][id])
				{
					const auto& block = std::get<FixedBlock>(term.block());
					uint32_t r = term.row(), c = term.column();
					auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2), cont_.evaluate(block).flatten());
					tmp.at(r, c) += val;
					if (c != r) tmp.at(c, r) += val;
				}
				for (uint32_t j = 0; j < n; ++j)
					for (uint32_t r = 0; r < sz; ++r)
						for (uint32_t c = 0; c < sz; ++c) mat[j + offsets_[i]].at(r, c) = std::move(tmp.at(r, c)[j]);
			}
			return mat;
		}
		[[nodiscard]] std::unique_ptr<ConformalScale<Real>> common_scale(uint32_t id,
		                                                                 const GeneralPrimaryOperator<Real>& op) const
		{
			auto include_odd = false;
			for (uint32_t i = 0; !include_odd && i < eqs_.size(); ++i)
				for (const auto& term : eqs_[i][id])
					if (std::get<GeneralBlock>(term.block()).include_odd())
					{
						include_odd = true;
						break;
					}
			auto ag = std::make_unique<ConformalScale<Real>>(op, cont_, include_odd);
			ag->set_gap(op.lower_bound(), op.upper_bound_safe());
			return ag;
		}
		[[nodiscard]] algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> make_cont_mat(
		    uint32_t id, const GeneralPrimaryOperator<Real>& op, std::unique_ptr<ConformalScale<Real>>& ag) const
		{
			auto sz = sz_[id];
			auto sp = ag->sample_points();
			algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> mat(N_);
			for (uint32_t i = 0; i < N_; ++i) mat[i] = algebra::Vector<algebra::Matrix<Real>>(sp.size(), {sz, sz});
			for (uint32_t i = 0; i < eqs_.size(); ++i)
			{
				if (eqs_[i][id].empty()) continue;
				auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
				for (uint32_t k = 0; k < sp.size(); ++k)
				{
					algebra::Matrix<algebra::Vector<Real>> tmp(sz, sz);
					for (uint32_t r = 0; r < sz; ++r)
						for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = algebra::Vector<Real>(n);
					for (const auto& term : eqs_[i][id])
					{
						uint32_t r = term.row(), c = term.column();
						const auto& block = std::get<GeneralBlock>(term.block()).fix_op(op);
						auto delta = ag->get_delta(sp[k]);
						auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2), cont_.evaluate(block, delta).flatten());
						tmp.at(r, c) += val;
						if (c != r) tmp.at(c, r) += val;
					}
					for (uint32_t j = 0; j < n; ++j)
						for (uint32_t r = 0; r < sz; ++r)
							for (uint32_t c = 0; c < sz; ++c)
								mat[j + offsets_[i]][k].at(r, c) = std::move(tmp.at(r, c)[j]);
				}
			}
			return mat;
		}
		// alpha maximizes alpha(norm)
		// and satisfies alpha(target) = N and alpha(sec) >= 0 for each sector sec (!= target, norm)
		[[nodiscard]] PolynomialProgram<Real> ope_maximize(const std::string& target, const std::string& norm, Real&& N,
		                                                   bool verbose = false) const
		{
			PolynomialProgram<Real> prg(N_);
			{
				if (verbose) std::cout << "[" << norm << "]" << std::endl;
				auto id = sectors_.at(norm);
				assert(sz_[id] == 1 || ope_[id].has_value());
				assert(ops_[id].empty());
				auto mat = make_disc_mat(id);
				algebra::Vector<Real> v(N_);
				for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
				prg.objective_constant() = Real(0);
				prg.objectives(std::move(v));
			}
			{
				if (verbose) std::cout << "[" << target << "]" << std::endl;
				auto id = sectors_.at(target);
				assert(sz_[id] == 1 || ope_[id].has_value());
				assert(ops_[id].empty());
				auto mat = make_disc_mat(id);
				algebra::Vector<Real> v(N_);
				for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
				prg.add_equation(std::move(v), std::move(N));
			}
			for (const auto& [sec, id] : sectors_)
			{
				if (sec == norm || sec == target) continue;
				if (verbose) std::cout << "[" << sec << "]" << std::endl;
				auto sz = sz_[id];
				if (ops_[id].empty())
					if (ope_[id].has_value())
					{
						auto m = make_disc_mat(id);
						algebra::Vector<Real> v(N_);
						for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
						prg.add_inequality(std::make_unique<Ineq>(N_, std::move(v), Real(0)));
					}
					else
						prg.add_inequality(
						    std::make_unique<Ineq>(N_, sz, make_disc_mat(id), algebra::Matrix<Real>(sz, sz)));
				else
					for (const auto& op : ops_[id])
					{
						if (verbose) std::cout << op.str() << std::endl;
						auto ag = common_scale(id, op);
						auto mat = make_cont_mat(id, op, ag);
						prg.add_inequality(std::make_unique<Ineq>(
						    N_, sz, std::move(ag), std::move(mat),
						    algebra::Vector<algebra::Matrix<Real>>(ag->max_degree() + 1, {sz, sz})));
					}
			}
			return prg;
		}

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
				ops_.emplace_back();
				ope_.emplace_back();
			}
		}
		void register_operator(const std::string& sec, const GeneralPrimaryOperator<Real>& op) &
		{
			ops_[sectors_[sec]].push_back(op);
		}
		void register_ope(const std::string& sec, algebra::Vector<Real>&& ope) &
		{
			auto id = sectors_[sec];
			assert(ope.size() == sz_[id]);
			assert(!ope_[id].has_value());
			ope_[id] = std::move(ope);
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
			syms_.push_back(sym);
			std::vector<std::vector<Entry<Real>>> eq(sectors_.size());
			offsets_.push_back(N_);
			N_ += algebra::function_dimension(cont_.lambda(), sym);
			for (; first != last; ++first)
			{
				const auto& [sec, entry] = *first;
				static_assert(std::is_convertible_v<decltype(sec), std::string>);
				static_assert(std::is_convertible_v<decltype(entry), Entry<Real>>);
				assert(std::visit([](auto b) { return b.symmetry(); }, entry.block()) == sym);
				auto id = sectors_.at(sec);
				assert(entry.row() < sz_[id] && entry.column() < sz_[id]);
				eq[id].push_back(entry);
			}
			eqs_.push_back(std::move(eq));
		}
		using FixedBlock = ConformalBlock<Real, PrimaryOperator<Real>>;
		using GeneralBlock = GeneralConformalBlock<Real>;
		using Ineq = PolynomialInequalityEvaluated<Real>;
		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha(norm) = 1 and alpha(sec) >= 0 for each sector sec (!= norm)
		// the size of matrices in norm sector must be 1
		[[nodiscard]] PolynomialProgram<Real> find_contradiction(const std::string& norm, bool verbose = false) const
		{
			PolynomialProgram<Real> prg(N_);
			prg.objective_constant() = Real(0);
			prg.objectives(algebra::Vector(N_, Real(0)));
			{
				if (verbose) std::cout << "[" << norm << "]" << std::endl;
				auto id = sectors_.at(norm);
				assert(sz_[id] == 1 || ope_[id].has_value());
				assert(ops_[id].empty());
				auto mat = make_disc_mat(id);
				algebra::Vector<Real> v(N_);
				for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
				prg.add_equation(std::move(v), Real(1));
			}
			for (const auto& [sec, id] : sectors_)
			{
				if (sec == norm) continue;
				if (verbose) std::cout << "[" << sec << "]" << std::endl;
				auto sz = sz_[id];
				if (ops_[id].empty())
					if (ope_[id].has_value())
					{
						auto m = make_disc_mat(id);
						algebra::Vector<Real> v(N_);
						for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
						prg.add_inequality(std::make_unique<Ineq>(N_, std::move(v), Real(0)));
					}
					else
						prg.add_inequality(
						    std::make_unique<Ineq>(N_, sz, make_disc_mat(id), algebra::Matrix<Real>(sz, sz)));
				else
					for (const auto& op : ops_[id])
					{
						if (verbose) std::cout << op.str() << std::endl;
						auto ag = common_scale(id, op);
						auto mat = make_cont_mat(id, op, ag);
						prg.add_inequality(std::make_unique<Ineq>(
						    N_, sz, std::move(ag), std::move(mat),
						    algebra::Vector<algebra::Matrix<Real>>(ag->max_degree() + 1, {sz, sz})));
					}
			}
			return prg;
		}
		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = 1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an upper bound on lambda, lambda ^ 2 <= -alpha(norm)
		[[nodiscard]] PolynomialProgram<Real> ope_maximize(const std::string& target, const std::string& norm,
		                                                   bool verbose = false) const
		{
			return ope_maximize(target, norm, Real(1), verbose);
		}
		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = -1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an lower bound on lambda, lambda ^ 2 >= alpha(norm)
		[[nodiscard]] PolynomialProgram<Real> ope_minimize(const std::string& target, const std::string& norm,
		                                                   bool verbose = false) const
		{
			return ope_maximize(target, norm, Real(-1), verbose);
		}
	};
}  // namespace qboot

#endif  // QBOOT_BOOTSTRAP_EQUATION_HPP_
