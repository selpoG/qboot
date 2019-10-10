#ifndef QBOOT_BOOTSTRAP_EQUATION_HPP_
#define QBOOT_BOOTSTRAP_EQUATION_HPP_

#include <algorithm>    // for min
#include <array>        // for array
#include <cstdint>      // for uint32_t
#include <functional>   // for function
#include <map>          // for map
#include <memory>       // for unique_ptr
#include <optional>     // for optional
#include <string>       // for string
#include <string_view>  // for string_view
#include <tuple>        // for tuple
#include <utility>      // for move
#include <vector>       // for vector

#include "qboot/algebra/complex_function.hpp"  // for FunctionSymmetry, function_dimension
#include "qboot/algebra/matrix.hpp"            // for Vector, Matrix
#include "qboot/block.hpp"                     // for ConformalBlock
#include "qboot/conformal_scale.hpp"           // for ConformalScale
#include "qboot/context.hpp"                   // for Context
#include "qboot/mp/rational.hpp"               // for rational
#include "qboot/mp/real.hpp"                   // for real
#include "qboot/polynomial_program.hpp"        // for PolynomialProgram
#include "qboot/primary_op.hpp"                // for GeneralPrimaryOperator, PrimaryOperator
#include "qboot/task_queue.hpp"                // for _event_base

namespace qboot
{
	// function object
	class _pole_selector
	{
		uint32_t numax_;

	public:
		explicit _pole_selector(uint32_t numax) : numax_(numax) {}
		// num of poles to pick
		uint32_t operator()(uint32_t spin) const { return numax_ + std::min(numax_, spin) / 2; }
	};
	// an element in bootstrap equation
	// coeff ope[r] ope[c] block
	class Entry
	{
		uint32_t r_, c_;
		mp::real coeff_;
		Block block_;

	public:
		void _reset() &&
		{
			r_ = c_ = 0;
			std::move(coeff_)._reset();
			std::visit([](auto&& v) { std::move(v)._reset(); }, std::move(block_));
		}
		Entry(uint32_t r, uint32_t c, const mp::real& coeff, const Block& block)
		    : r_(r), c_(c), coeff_(coeff), block_(block)
		{
		}
		Entry(const mp::real& coeff, const Block& block) : Entry(0, 0, coeff, block) {}
		Entry(uint32_t r, uint32_t c, const Block& block) : Entry(r, c, mp::real(1), block) {}
		explicit Entry(const Block& block) : Entry(0, 0, mp::real(1), block) {}
		[[nodiscard]] const Block& block() const { return block_; }
		[[nodiscard]] uint32_t row() const { return r_; }
		[[nodiscard]] uint32_t column() const { return c_; }
		[[nodiscard]] const mp::real& coeff() const { return coeff_; }
	};
	enum class SectorType : bool
	{
		Continuous = false,
		Discrete = true
	};
	class Sector
	{
		friend class BootstrapEquation;
		std::string name_;
		uint32_t sz_;
		SectorType type_;
		// tuple(spin, lower bound of delta, upper bound of delta)
		std::vector<std::tuple<uint32_t, std::optional<mp::real>, std::optional<mp::real>>> op_args_{};
		std::vector<GeneralPrimaryOperator> ops_{};
		std::optional<algebra::Vector<mp::real>> ope_{};
		static std::optional<algebra::Vector<mp::real>> clone(const std::optional<algebra::Vector<mp::real>>& ope)
		{
			if (ope.has_value()) return ope.value().clone();
			return {};
		}
		void set_operators(const mp::rational& epsilon, const std::function<uint32_t(uint32_t)>& num_poles) &
		{
			assert(type_ == SectorType::Continuous);
			ops_ = {};
			for (const auto& [sp, lb, ub] : op_args_)
			{
				if (ub.has_value())
					ops_.emplace_back(sp, num_poles(sp), epsilon, lb.value(), ub.value());
				else if (lb.has_value())
					ops_.emplace_back(sp, num_poles(sp), epsilon, lb.value());
				else
					ops_.emplace_back(sp, num_poles(sp), epsilon);
			}
		}

	public:
		void _reset() &&
		{
			std::string{}.swap(name_);
			sz_ = 0;
			type_ = SectorType::Continuous;
			std::vector<std::tuple<uint32_t, std::optional<mp::real>, std::optional<mp::real>>>{}.swap(op_args_);
			std::vector<GeneralPrimaryOperator>{}.swap(ops_);
			ope_.reset();
		}
		Sector(std::string_view name, uint32_t size, SectorType type = SectorType::Discrete)
		    : name_(name), sz_(size), type_(type)
		{
			assert(sz_ > 0);
		}
		Sector(std::string_view name, uint32_t size, algebra::Vector<mp::real>&& ope)
		    : name_(name), sz_(size), type_(SectorType::Discrete), ope_{std::move(ope)}
		{
			assert(ope_->size() == sz_);
			assert(sz_ > 0);
		}
		~Sector() = default;
		Sector(const Sector& s)
		    : name_(s.name_), sz_(s.sz_), type_(s.type_), op_args_(s.op_args_), ops_(s.ops_), ope_(clone(s.ope_))
		{
		}
		Sector(Sector&&) = default;
		Sector& operator=(const Sector& s)
		{
			if (this != &s) *this = Sector(s);
			return *this;
		}
		Sector& operator=(Sector&&) = default;
		[[nodiscard]] const std::string& name() const { return name_; }
		[[nodiscard]] uint32_t size() const { return sz_; }
		[[nodiscard]] SectorType type() const { return type_; }
		[[nodiscard]] bool is_matrix() const { return sz_ > 1 && !ope_.has_value(); }
		[[nodiscard]] const std::optional<algebra::Vector<mp::real>>& ope() const { return ope_; }
		// delta in [unitarity bound, inf)
		void add_op(uint32_t spin) &
		{
			assert(type_ == SectorType::Continuous);
			op_args_.emplace_back(spin, std::optional<mp::real>{}, std::optional<mp::real>{});
		}
		// delta in [lb, inf)
		void add_op(uint32_t spin, const mp::real& lb) &
		{
			assert(type_ == SectorType::Continuous);
			op_args_.emplace_back(spin, lb, std::optional<mp::real>{});
		}
		// delta in [lb, ub)
		void add_op(uint32_t spin, const mp::real& lb, const mp::real& ub) &
		{
			assert(type_ == SectorType::Continuous);
			ub.isinf() ? add_op(spin, lb) : void(op_args_.emplace_back(spin, lb, ub));
		}
	};
	class Equation;
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
	// an instance of BootstrapEquation refers to Context,
	// so you must gurantee that Context survives longer than BootstrapEquation
	class BootstrapEquation
	{
		const Context& cont_;
		std::vector<Sector> sectors_{};
		// maps from sector name to its unique id
		std::map<std::string, uint32_t, std::less<>> sector_id_{};
		std::vector<Equation> eqs_{};
		// total dimension (index of equations, index of derivatives)
		uint32_t N_ = 0;

		[[nodiscard]] algebra::Vector<mp::real> make_disc_mat_v(uint32_t id) const
		{
			auto m = make_disc_mat(id);
			algebra::Vector<mp::real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, std::move(m[i]));
			return v;
		}
		[[nodiscard]] mp::real take_element(uint32_t id, algebra::Matrix<mp::real>&& m) const
		{
			const auto& ope = sector(id).ope();
			if (ope) return m.inner_product(ope.value());
			return std::move(m.at(0, 0));
		}
		[[nodiscard]] algebra::Vector<algebra::Matrix<mp::real>> make_disc_mat(uint32_t id) const;
		[[nodiscard]] std::unique_ptr<ConformalScale> common_scale(uint32_t id, const GeneralPrimaryOperator& op) const;
		[[nodiscard]] algebra::Vector<algebra::Vector<algebra::Matrix<mp::real>>> make_cont_mat(
		    uint32_t id, const GeneralPrimaryOperator& op, const std::unique_ptr<ConformalScale>& ag) const;

		// alpha maximizes alpha(norm)
		// and satisfies alpha(target) = N and alpha(sec) >= 0 for each sector sec (!= target, norm)
		[[nodiscard]] PolynomialProgram ope_maximize(std::string_view target, std::string_view norm, mp::real&& N,
		                                             uint32_t parallel,
		                                             const std::unique_ptr<_event_base>& event) const;
		void add_ineqs(PolynomialProgram* prg, const std::function<bool(const std::string&)>& filter, uint32_t parallel,
		               const std::unique_ptr<_event_base>& event) const;

	public:
		void _reset() &&
		{
			std::vector<Sector>{}.swap(sectors_);
			std::map<std::string, uint32_t, std::less<>>{}.swap(sector_id_);
			std::vector<Equation>{}.swap(eqs_);
			N_ = 0;
		}
		BootstrapEquation(const Context& cont, const std::vector<Sector>& sectors, uint32_t numax)
		    : BootstrapEquation(cont, sectors, _pole_selector(numax))
		{
		}
		BootstrapEquation(const Context& cont, std::vector<Sector>&& sectors, uint32_t numax)
		    : BootstrapEquation(cont, std::move(sectors), _pole_selector(numax))
		{
		}
		BootstrapEquation(const Context& cont, const std::vector<Sector>& sectors,
		                  const std::function<uint32_t(uint32_t)>& num_poles)
		    : BootstrapEquation(cont, std::vector(sectors), num_poles)
		{
		}
		// num_poles: spin -> num of poles
		BootstrapEquation(const Context& cont, std::vector<Sector>&& sectors,
		                  const std::function<uint32_t(uint32_t)>& num_poles)
		    : cont_(cont), sectors_(std::move(sectors))
		{
			for (uint32_t id = 0; id < sectors_.size(); ++id)
			{
				sector_id_[sectors_[id].name()] = id;
				if (sectors_[id].type() == SectorType::Continuous)
					sectors_[id].set_operators(cont.epsilon(), num_poles);
			}
		}
		void add_equation(const Equation& eq) &
		{
			assert(N_ == 0);
			eqs_.push_back(eq);
		}
		void add_equation(Equation&& eq) &
		{
			assert(N_ == 0);
			eqs_.push_back(std::move(eq));
		}
		// call this to finish add_equation
		void finish() &;
		[[nodiscard]] uint32_t lambda() const { return cont_.lambda(); }
		[[nodiscard]] uint32_t num_sectors() const { return uint32_t(sectors_.size()); }
		[[nodiscard]] uint32_t get_id(std::string_view sector_name) const
		{
			return sector_id_.find(sector_name)->second;
		}
		[[nodiscard]] const Sector& sector(uint32_t id) const { return sectors_[id]; }

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha(norm) = 1 and alpha(sec) >= 0 for each sector sec (!= norm)
		// the size of matrices in norm sector must be 1
		[[nodiscard]] PolynomialProgram find_contradiction(std::string_view norm, uint32_t parallel = 1,
		                                                   const std::unique_ptr<_event_base>& event = {}) const;

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = 1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an upper bound on lambda, lambda ^ 2 <= -alpha(norm)
		[[nodiscard]] PolynomialProgram ope_maximize(std::string_view target, std::string_view norm,
		                                             uint32_t parallel = 1,
		                                             const std::unique_ptr<_event_base>& event = {}) const
		{
			return ope_maximize(target, norm, mp::real(1), parallel, event);
		}

		// create a PolynomialProgram which finds a linear functional alpha
		// s.t. alpha maximizes alpha(norm)
		// and satisfies alpha(target) = -1 and alpha(sec) >= 0 for each sector sec (!= target, norm)
		// the size of matrices in target, norm sector must be 1
		// such a alpha gives an lower bound on lambda, lambda ^ 2 >= alpha(norm)
		[[nodiscard]] PolynomialProgram ope_minimize(std::string_view target, std::string_view norm,
		                                             uint32_t parallel = 1,
		                                             const std::unique_ptr<_event_base>& event = {}) const
		{
			return ope_maximize(target, norm, mp::real(-1), parallel, event);
		}
	};
	using Externals = std::array<PrimaryOperator, 4>;
	class Equation
	{
		const BootstrapEquation& boot_;
		algebra::FunctionSymmetry sym_;
		uint32_t dim_;
		// terms_[id]: terms in id-th sector
		std::vector<std::vector<Entry>> terms_{};

	public:
		void _reset() &&
		{
			sym_ = algebra::FunctionSymmetry::Mixed;
			dim_ = 0;
			std::vector<std::vector<Entry>>{}.swap(terms_);
		}
		Equation(const BootstrapEquation& boot, algebra::FunctionSymmetry sym)
		    : boot_(boot), sym_(sym), dim_(algebra::function_dimension(boot.lambda(), sym)), terms_(boot_.num_sectors())
		{
		}
		[[nodiscard]] algebra::FunctionSymmetry symmetry() const { return sym_; }
		[[nodiscard]] uint32_t dimension() const { return dim_; }
		const std::vector<Entry>& operator[](uint32_t id) const { return terms_[id]; }
		void add(std::string_view sec, uint32_t r, uint32_t c, const mp::real& coeff, const PrimaryOperator& o,
		         const Externals& os)
		{
			auto id = boot_.get_id(sec);
			assert(r < boot_.sector(id).size() && c < boot_.sector(id).size());
			terms_[id].emplace_back(r, c, coeff, ConformalBlock<PrimaryOperator>(o, os[0], os[1], os[2], os[3], sym_));
		}
		void add(std::string_view sec, const mp::real& coeff, const PrimaryOperator& o, const Externals& os)
		{
			add(sec, 0, 0, coeff, o, os);
		}
		void add(std::string_view sec, uint32_t r, uint32_t c, const PrimaryOperator& o, const Externals& os)
		{
			add(sec, r, c, mp::real(1), o, os);
		}
		void add(std::string_view sec, const PrimaryOperator& o, const Externals& os) { add(sec, 0, 0, o, os); }
		void add(std::string_view sec, uint32_t r, uint32_t c, const mp::real& coeff, const Externals& os)
		{
			auto id = boot_.get_id(sec);
			assert(r < boot_.sector(id).size() && c < boot_.sector(id).size());
			terms_[id].emplace_back(r, c, coeff, GeneralConformalBlock(os[0], os[1], os[2], os[3], sym_));
		}
		void add(std::string_view sec, const mp::real& coeff, const Externals& os) { add(sec, 0, 0, coeff, os); }
		void add(std::string_view sec, uint32_t r, uint32_t c, const Externals& os) { add(sec, r, c, mp::real(1), os); }
		void add(std::string_view sec, const Externals& os) { add(sec, 0, 0, os); }
	};
}  // namespace qboot

#endif  // QBOOT_BOOTSTRAP_EQUATION_HPP_
