#include "bootstrap_equation.hpp"

#include <iostream>  // for cout, endl
#include <memory>    // for unique_ptr
#include <string>    // for string
#include <utility>   // for move

using algebra::Vector, algebra::Matrix;
using mpfr::real;
using std::move, std::unique_ptr, std::make_unique, std::cout, std::endl, std::string;

namespace qboot
{
	using FixedBlock = ConformalBlock<PrimaryOperator>;
	using GeneralBlock = GeneralConformalBlock;
	using Ineq = PolynomialInequalityEvaluated;
	[[nodiscard]] Vector<Matrix<real>> BootstrapEquation::make_disc_mat(uint32_t id) const
	{
		auto sz = sz_[id];
		Vector<Matrix<real>> mat(N_);
		for (uint32_t i = 0; i < N_; ++i) mat[i] = {sz, sz};
		for (uint32_t i = 0; i < eqs_.size(); ++i)
		{
			if (eqs_[i][id].empty()) continue;
			auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
			Matrix<Vector<real>> tmp(sz, sz);
			for (uint32_t r = 0; r < sz; ++r)
				for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = Vector<real>{n};
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
					for (uint32_t c = 0; c < sz; ++c) mat[j + offsets_[i]].at(r, c) = move(tmp.at(r, c)[j]);
		}
		return mat;
	}
	[[nodiscard]] unique_ptr<ConformalScale> BootstrapEquation::common_scale(uint32_t id,
	                                                                         const GeneralPrimaryOperator& op) const
	{
		auto include_odd = false;
		for (uint32_t i = 0; !include_odd && i < eqs_.size(); ++i)
			for (const auto& term : eqs_[i][id])
				if (std::get<GeneralBlock>(term.block()).include_odd())
				{
					include_odd = true;
					break;
				}
		auto ag = make_unique<ConformalScale>(op, cont_, include_odd);
		ag->set_gap(op.lower_bound(), op.upper_bound_safe());
		return ag;
	}
	[[nodiscard]] Vector<Vector<Matrix<real>>> BootstrapEquation::make_cont_mat(uint32_t id,
	                                                                            const GeneralPrimaryOperator& op,
	                                                                            unique_ptr<ConformalScale>* ag) const
	{
		auto sz = sz_[id];
		auto sp = (*ag)->sample_points();
		Vector<Vector<Matrix<real>>> mat(N_);
		for (uint32_t i = 0; i < N_; ++i) mat[i] = Vector<Matrix<real>>(sp.size(), {sz, sz});
		for (uint32_t i = 0; i < eqs_.size(); ++i)
		{
			if (eqs_[i][id].empty()) continue;
			auto n = algebra::function_dimension(cont_.lambda(), syms_[i]);
			for (uint32_t k = 0; k < sp.size(); ++k)
			{
				Matrix<Vector<real>> tmp(sz, sz);
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = Vector<real>(n);
				for (const auto& term : eqs_[i][id])
				{
					uint32_t r = term.row(), c = term.column();
					const auto& block = std::get<GeneralBlock>(term.block()).fix_op(op);
					auto delta = (*ag)->get_delta(sp[k]);
					auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2), cont_.evaluate(block, delta).flatten());
					tmp.at(r, c) += val;
					if (c != r) tmp.at(c, r) += val;
				}
				for (uint32_t j = 0; j < n; ++j)
					for (uint32_t r = 0; r < sz; ++r)
						for (uint32_t c = 0; c < sz; ++c) mat[j + offsets_[i]][k].at(r, c) = move(tmp.at(r, c)[j]);
			}
		}
		return mat;
	}

	[[nodiscard]] PolynomialProgram BootstrapEquation::ope_maximize(const string& target, const string& norm, real&& N,
	                                                                bool verbose) const
	{
		PolynomialProgram prg(N_);
		{
			if (verbose) cout << "[" << norm << "]" << endl;
			auto id = sectors_.at(norm);
			assert(sz_[id] == 1 || ope_[id].has_value());
			assert(ops_[id].empty());
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.objective_constant() = real(0);
			prg.objectives(move(v));
		}
		{
			if (verbose) cout << "[" << target << "]" << endl;
			auto id = sectors_.at(target);
			assert(sz_[id] == 1 || ope_[id].has_value());
			assert(ops_[id].empty());
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.add_equation(move(v), move(N));
		}
		for (const auto& [sec, id] : sectors_)
		{
			if (sec == norm || sec == target) continue;
			if (verbose) cout << "[" << sec << "]" << endl;
			auto sz = sz_[id];
			if (ops_[id].empty())
				if (ope_[id].has_value())
				{
					auto m = make_disc_mat(id);
					Vector<real> v(N_);
					for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
					prg.add_inequality(make_unique<Ineq>(N_, move(v), real(0)));
				}
				else
					prg.add_inequality(make_unique<Ineq>(N_, sz, make_disc_mat(id), Matrix<real>(sz, sz)));
			else
				for (const auto& op : ops_[id])
				{
					if (verbose) cout << op.str() << endl;
					auto ag = common_scale(id, op);
					auto mat = make_cont_mat(id, op, &ag);
					prg.add_inequality(make_unique<Ineq>(N_, sz, move(ag), move(mat),
					                                     Vector<Matrix<real>>(ag->max_degree() + 1, {sz, sz})));
				}
		}
		return prg;
	}

	[[nodiscard]] PolynomialProgram BootstrapEquation::find_contradiction(const string& norm, bool verbose) const
	{
		PolynomialProgram prg(N_);
		prg.objective_constant() = real(0);
		prg.objectives(Vector(N_, real(0)));
		{
			if (verbose) cout << "[" << norm << "]" << endl;
			auto id = sectors_.at(norm);
			assert(sz_[id] == 1 || ope_[id].has_value());
			assert(ops_[id].empty());
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.add_equation(move(v), real(1));
		}
		for (const auto& [sec, id] : sectors_)
		{
			if (sec == norm) continue;
			if (verbose) cout << "[" << sec << "]" << endl;
			auto sz = sz_[id];
			if (ops_[id].empty())
				if (ope_[id].has_value())
				{
					auto m = make_disc_mat(id);
					Vector<real> v(N_);
					for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
					prg.add_inequality(make_unique<Ineq>(N_, move(v), real(0)));
				}
				else
					prg.add_inequality(make_unique<Ineq>(N_, sz, make_disc_mat(id), Matrix<real>(sz, sz)));
			else
				for (const auto& op : ops_[id])
				{
					if (verbose) cout << op.str() << endl;
					auto ag = common_scale(id, op);
					auto mat = make_cont_mat(id, op, &ag);
					prg.add_inequality(make_unique<Ineq>(N_, sz, move(ag), move(mat),
					                                     Vector<Matrix<real>>(ag->max_degree() + 1, {sz, sz})));
				}
		}
		return prg;
	}
}  // namespace qboot
