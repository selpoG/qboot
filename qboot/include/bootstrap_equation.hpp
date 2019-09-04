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
#include "polynomial_program.hpp"  // for PolynomialProgramming
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
				ops_.emplace_back();
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
			syms_.push_back(sym);
			std::vector<std::vector<Entry<Real>>> eq(sectors_.size());
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
		// Func: spin |-> num of poles
		// uint32_t -> uint32_t
		template <class Func>
		[[nodiscard]] PolynomialProgramming<Real> create_pmp(const std::string& norm, Func num_poles) const
		{
			PolynomialProgramming<Real> prg(N_);
			prg.objective_constant() = Real(0);
			prg.objectives(algebra::Vector(N_, Real(0)));
			{
				auto id = sectors_.at(norm);
				assert(sz_[id] == 1);
				assert(ops_[id].empty());
				algebra::Vector<Real> v(N_);
				uint32_t p = 0;
				for (uint32_t i = 0; i < eqs_.size(); ++i)
				{
					auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
					if (!eqs_[i][id].empty())
					{
						algebra::Vector<Real> tmp(n);
						for (const auto& term : eqs_[i][id])
							tmp +=
							    mul_scalar(term.coeff(), cont_.evaluate(std::get<FixedBlock>(term.block())).flatten());
						for (uint32_t j = 0; j < n; ++j) v[j + p] = std::move(tmp[j]);
					}
					p += n;
				}
				prg.add_equation(std::move(v), Real(1));
			}
			for (const auto& [sec, id] : sectors_)
			{
				if (sec == norm) continue;
				auto sz = sz_[id];
				if (ops_[id].empty())
				{
					// point-like
					algebra::Vector<algebra::Matrix<Real>> mat(N_);
					for (uint32_t i = 0; i < N_; ++i) mat[i] = {sz, sz};
					uint32_t p = 0;
					for (uint32_t i = 0; i < eqs_.size(); ++i)
					{
						auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
						if (!eqs_[i][id].empty())
						{
							algebra::Matrix<algebra::Vector<Real>> tmp(sz, sz);
							for (uint32_t r = 0; r < sz; ++r)
								for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = algebra::Vector<Real>(n);
							for (const auto& term : eqs_[i][id])
							{
								uint32_t r = term.row(), c = term.column();
								auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2),
								                      cont_.evaluate(std::get<FixedBlock>(term.block())).flatten());
								tmp.at(r, c) += val;
								if (c != r) tmp.at(c, r) += val;
							}
							for (uint32_t j = 0; j < n; ++j)
								for (uint32_t r = 0; r < sz; ++r)
									for (uint32_t c = 0; c < sz; ++c) mat[j + p].at(r, c) = std::move(tmp.at(r, c)[j]);
						}
						p += n;
					}
					prg.add_inequality(std::make_unique<Ineq>(N_, sz, std::move(mat), algebra::Matrix<Real>(sz, sz)));
				}
				else
					for (const auto& op : ops_[id])
					{
						auto ag = std::make_unique<ConformalScale<Real>>(num_poles(op.spin()), op.spin(), cont_, false);
						ag->set_gap(op.lower_bound(), op.upper_bound_safe());
						for (uint32_t i = 0; !ag->odd_included() && i < eqs_.size(); ++i)
							for (const auto& term : eqs_[i][id])
								if (std::get<GeneralBlock>(term.block()).include_odd())
								{
									ag = std::make_unique<ConformalScale<Real>>(num_poles(op.spin()), op.spin(), cont_,
									                                            true);
									ag->set_gap(op.lower_bound(), op.upper_bound_safe());
									break;
								}
						auto sp = ag->sample_points();
						algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> mat(N_);
						for (uint32_t i = 0; i < N_; ++i)
							mat[i] = algebra::Vector<algebra::Matrix<Real>>(sp.size(), {sz, sz});
						uint32_t p = 0;
						for (uint32_t i = 0; i < eqs_.size(); ++i)
						{
							auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
							if (!eqs_[i][id].empty())
							{
								for (uint32_t k = 0; k < sp.size(); ++k)
								{
									algebra::Matrix<algebra::Vector<Real>> tmp(sz, sz);
									for (uint32_t r = 0; r < sz; ++r)
										for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = algebra::Vector<Real>(n);
									for (const auto& term : eqs_[i][id])
									{
										uint32_t r = term.row(), c = term.column();
										auto block = std::get<GeneralBlock>(term.block()).fix_op(op);
										auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2),
										                      cont_.evaluate(block, ag->get_delta(sp[k])).flatten());
										tmp.at(r, c) += val;
										if (c != r) tmp.at(c, r) += val;
									}
									for (uint32_t j = 0; j < n; ++j)
										for (uint32_t r = 0; r < sz; ++r)
											for (uint32_t c = 0; c < sz; ++c)
												mat[j + p][k].at(r, c) = std::move(tmp.at(r, c)[j]);
								}
							}
							p += n;
						}
						prg.add_inequality(
						    std::make_unique<Ineq>(N_, sz, std::move(ag), std::move(mat),
						                           algebra::Vector<algebra::Matrix<Real>>(sp.size(), {sz, sz})));
					}
			}
			return prg;
		}
	};
}  // namespace qboot

#endif  // QBOOT_BOOTSTRAP_EQUATION_HPP_
