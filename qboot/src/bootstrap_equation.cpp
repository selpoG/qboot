#include "bootstrap_equation.hpp"

#include <future>    // for future
#include <iostream>  // for cout, endl
#include <memory>    // for unique_ptr
#include <string>    // for string
#include <utility>   // for move

#include "task_queue.hpp"  // for TaskQueue

using algebra::Vector, algebra::Matrix;
using mp::real;
using std::move, std::unique_ptr, std::make_unique, std::cout, std::endl, std::string, std::string_view, std::vector;

namespace qboot
{
	using FixedBlock = ConformalBlock<PrimaryOperator>;
	using GeneralBlock = GeneralConformalBlock;
	using Ineq = PolynomialInequalityEvaluated;
	[[nodiscard]] Vector<Matrix<real>> BootstrapEquation::make_disc_mat(uint32_t id) const
	{
		auto sz = sector(id).size();
		Vector<Matrix<real>> mat(N_);
		for (uint32_t i = 0; i < N_; ++i) mat[i] = {sz, sz};
		uint32_t p = 0;
		for (const auto& eq : eqs_)
		{
			auto n = eq.dimension();
			if (eq[id].empty())
			{
				p += n;
				continue;
			}
			Matrix<Vector<real>> tmp(sz, sz);
			for (uint32_t r = 0; r < sz; ++r)
				for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = Vector<real>{n};
			for (const auto& term : eq[id])
			{
				const auto& block = std::get<FixedBlock>(term.block());
				uint32_t r = term.row(), c = term.column();
				auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2), cont_.evaluate(block).flatten());
				if (c != r) tmp.at(c, r) += val;
				tmp.at(r, c) += val;
			}
			for (uint32_t j = 0; j < n; ++j)
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) mat[j + p].at(r, c) = move(tmp.at(r, c)[j]);
			p += n;
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
		auto sz = sector(id).size();
		auto sp = (*ag)->sample_points();
		Vector<Vector<Matrix<real>>> mat(N_);
		for (uint32_t i = 0; i < N_; ++i) mat[i] = Vector<Matrix<real>>(sp.size(), {sz, sz});
		uint32_t p = 0;
		for (const auto& eq : eqs_)
		{
			auto n = eq.dimension();
			if (eq[id].empty())
			{
				p += n;
				continue;
			}
			for (uint32_t k = 0; k < sp.size(); ++k)
			{
				Matrix<Vector<real>> tmp(sz, sz);
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) tmp.at(r, c) = Vector<real>(n);
				for (const auto& term : eq[id])
				{
					uint32_t r = term.row(), c = term.column();
					const auto& block = std::get<GeneralBlock>(term.block()).fix_op(op);
					auto delta = (*ag)->get_delta(sp[k]);
					auto val = mul_scalar(term.coeff() / (r == c ? 1 : 2), cont_.evaluate(block, delta).flatten());
					if (c != r) tmp.at(c, r) += val;
					tmp.at(r, c) += val;
				}
				for (uint32_t j = 0; j < n; ++j)
					for (uint32_t r = 0; r < sz; ++r)
						for (uint32_t c = 0; c < sz; ++c) mat[j + p][k].at(r, c) = move(tmp.at(r, c)[j]);
			}
			p += n;
		}
		return mat;
	}

	[[nodiscard]] PolynomialProgram BootstrapEquation::ope_maximize(string_view target, string_view norm, real&& N,
	                                                                uint32_t parallel, bool verbose) const
	{
		assert(N_ > 0);
		PolynomialProgram prg(N_);
		{
			if (verbose) cout << "[" << norm << "]" << endl;
			auto id = sector_id_.find(norm)->second;
			assert(!sector(id).is_matrix());
			assert(sector(id).type() == SectorType::Discrete);
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.objective_constant() = real(0);
			prg.objectives(move(v));
		}
		{
			if (verbose) cout << "[" << target << "]" << endl;
			auto id = sector_id_.find(target)->second;
			assert(!sector(id).is_matrix());
			assert(sector(id).type() == SectorType::Discrete);
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.add_equation(move(v), move(N));
		}
		TaskQueue q(parallel);
		vector<std::future<unique_ptr<Ineq>>> ineqs;
		for (const auto& [sec, id] : sector_id_)
		{
			if (sec == norm || sec == target) continue;
			if (verbose) cout << "[" << sec << "]" << endl;
			auto sz = sector(id).size();
			if (sector(id).type() == SectorType::Discrete)
				if (sector(id).is_matrix())
					ineqs.push_back(q.push([this, id = id, sz]() {
						return make_unique<Ineq>(N_, sz, make_disc_mat(id), Matrix<real>(sz, sz));
					}));
				else
					ineqs.push_back(q.push([this, id = id]() {
						auto m = make_disc_mat(id);
						Vector<real> v(N_);
						for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
						return make_unique<Ineq>(N_, move(v), real(0));
					}));
			else
				for (const auto& op : sector(id).ops_)
					ineqs.push_back(q.push([this, id = id, sz, &op]() {
						auto ag = common_scale(id, op);
						auto mat = make_cont_mat(id, op, &ag);
						return make_unique<Ineq>(N_, sz, move(ag), move(mat),
						                         Vector<Matrix<real>>(ag->max_degree() + 1, {sz, sz}));
					}));
		}
		for (auto&& x : ineqs) prg.add_inequality(x.get());
		return prg;
	}

	void BootstrapEquation::finish() &
	{
		for (const auto& eq : eqs_) N_ += eq.dimension();
	}

	[[nodiscard]] PolynomialProgram BootstrapEquation::find_contradiction(string_view norm, uint32_t parallel,
	                                                                      bool verbose) const
	{
		assert(N_ > 0);
		PolynomialProgram prg(N_);
		prg.objective_constant() = real(0);
		prg.objectives(Vector(N_, real(0)));
		{
			if (verbose) cout << "[" << norm << "]" << endl;
			auto id = sector_id_.find(norm)->second;
			assert(!sector(id).is_matrix());
			assert(sector(id).type() == SectorType::Discrete);
			auto mat = make_disc_mat(id);
			Vector<real> v(N_);
			for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, mat[i]);
			prg.add_equation(move(v), real(1));
		}
		TaskQueue q(parallel);
		vector<std::future<unique_ptr<Ineq>>> ineqs;
		for (const auto& [sec, id] : sector_id_)
		{
			if (sec == norm) continue;
			if (verbose) cout << "[" << sec << "]" << endl;
			auto sz = sector(id).size();
			if (sector(id).type() == SectorType::Discrete)
				if (sector(id).is_matrix())
					ineqs.push_back(q.push([this, id = id, sz]() {
						return make_unique<Ineq>(N_, sz, make_disc_mat(id), Matrix<real>(sz, sz));
					}));
				else
					ineqs.push_back(q.push([this, id = id]() {
						auto m = make_disc_mat(id);
						Vector<real> v(N_);
						for (uint32_t i = 0; i < N_; ++i) v[i] = take_element(id, m[i]);
						return make_unique<Ineq>(N_, move(v), real(0));
					}));
			else
				for (const auto& op : sector(id).ops_)
					ineqs.push_back(q.push([this, id = id, sz, &op]() {
						auto ag = common_scale(id, op);
						auto mat = make_cont_mat(id, op, &ag);
						return make_unique<Ineq>(N_, sz, move(ag), move(mat),
						                         Vector<Matrix<real>>(ag->max_degree() + 1, {sz, sz}));
					}));
		}
		for (auto&& x : ineqs) prg.add_inequality(x.get());
		return prg;
	}
}  // namespace qboot
