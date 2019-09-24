#include "polynomial_program.hpp"

#include <future>   // for future, promise
#include <memory>   // for unique_ptr
#include <utility>  // for move
#include <vector>   // for vector

#include "task_queue.hpp"  // for TaskQueue

using algebra::Vector, algebra::Matrix, algebra::Polynomial;
using mp::real;
using std::move, std::unique_ptr, std::make_unique, std::vector;

namespace qboot
{
	PolynomialInequality::~PolynomialInequality() = default;
	PolynomialInequalityEvaluated::~PolynomialInequalityEvaluated() = default;
	PolynomialInequalityWithCoeffs::~PolynomialInequalityWithCoeffs() = default;

	[[nodiscard]] Matrix<Polynomial> PolynomialInequalityEvaluated::as_polynomial(const Vector<Matrix<real>>& vals)
	{
		uint32_t deg = vals.size() - 1;
		Vector<Matrix<real>> ev(deg + 1);
		const auto& chi = PolynomialInequality::get_scale();
		auto xs = chi->sample_points();
		for (uint32_t k = 0; k <= deg; ++k) ev[k] = vals[k] / chi->eval(xs[k]);
		return algebra::polynomial_interpolate(ev, xs);
	}

	PolynomialInequalityEvaluated::PolynomialInequalityEvaluated(uint32_t N, unique_ptr<ScaleFactor>&& scale,
	                                                             Vector<Vector<real>>&& mat, Vector<real>&& target)
	    : PolynomialInequality(N, 1, move(scale)), mat_(N)
	{
		uint32_t deg = PolynomialInequality::get_scale()->max_degree();
		for (uint32_t i = 0; i < N; ++i)
		{
			assert(mat[i].size() == deg + 1);
			mat_[i] = Vector<Matrix<real>>{deg + 1};
			for (uint32_t k = 0; k <= deg; ++k)
			{
				mat_[i][k] = {1, 1};
				mat_[i][k].at(0, 0) = move(mat[i][k]);
			}
		}
		assert(target.size() == deg + 1);
		target_ = Vector<Matrix<real>>{deg + 1};
		for (uint32_t k = 0; k <= deg; ++k)
		{
			target_[k] = {1, 1};
			target_[k].at(0, 0) = move(target[k]);
		}
	}

	PolynomialInequalityEvaluated::PolynomialInequalityEvaluated(uint32_t N, Vector<real>&& mat, real&& target)
	    : PolynomialInequality(N, 1, make_unique<TrivialScale>()), mat_(N), target_(1)
	{
		assert(mat.size() == N);
		for (uint32_t i = 0; i < N; ++i)
		{
			mat_[i] = Vector<Matrix<real>>(1);
			mat_[i][0] = {1, 1};
			mat_[i][0].at(0, 0) = move(mat[i]);
		}
		target_[0] = {1, 1};
		target_[0].at(0, 0) = move(target);
	}

	PolynomialInequalityEvaluated::PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, Vector<Matrix<real>>&& mat,
	                                                             Matrix<real>&& target)
	    : PolynomialInequality(N, sz, make_unique<TrivialScale>()), mat_(N), target_(1)
	{
		assert(mat.size() == N);
		for (uint32_t i = 0; i < N; ++i)
		{
			assert(mat[i].is_square() && mat[i].row() == sz);
			mat_[i] = Vector<Matrix<real>>(1);
			mat_[i][0] = {sz, sz};
			for (uint32_t r = 0; r < sz; ++r)
				for (uint32_t c = 0; c < sz; ++c) mat_[i][0].at(r, c) = move(mat[i].at(r, c));
		}
		assert(target.is_square() && target.row() == sz);
		target_[0] = {sz, sz};
		for (uint32_t r = 0; r < sz; ++r)
			for (uint32_t c = 0; c < sz; ++c) target_[0].at(r, c) = move(target.at(r, c));
	}

	PolynomialInequalityEvaluated::PolynomialInequalityEvaluated(uint32_t N, uint32_t sz,
	                                                             std::unique_ptr<ScaleFactor>&& scale,
	                                                             Vector<Vector<Matrix<real>>>&& mat,
	                                                             Vector<Matrix<real>>&& target)
	    : PolynomialInequality(N, sz, move(scale)), mat_(move(mat)), target_(move(target))
	{
		uint32_t deg = PolynomialInequality::get_scale()->max_degree();
		assert(mat_.size() == N);
		for (uint32_t i = 0; i < N; ++i)
		{
			assert(mat_[i].size() == deg + 1);
			for (uint32_t k = 0; k <= deg; ++k) assert(mat_[i][k].is_square() && mat_[i][k].row() == sz);
		}
		assert(target_.size() == deg + 1);
		for (uint32_t k = 0; k <= deg; ++k) assert(target_[k].is_square() && target_[k].row() == sz);
	}

	PolynomialInequalityWithCoeffs::PolynomialInequalityWithCoeffs(uint32_t N, unique_ptr<ScaleFactor>&& scale,
	                                                               Vector<Polynomial>&& mat, Polynomial&& target)
	    : PolynomialInequality(N, 1, move(scale))
	{
		[[maybe_unused]] auto deg = int32_t(PolynomialInequality::get_scale()->max_degree());
		assert(mat.size() == N);
		mat_ = Vector<Matrix<Polynomial>>(N);
		for (uint32_t i = 0; i < N; ++i)
		{
			mat_[i] = {1, 1};
			assert(mat[i].degree() <= deg);
			mat_[i].at(0, 0) = move(mat[i]);
		}
		target_ = {1, 1};
		assert(target.degree() <= deg);
		target_.at(0, 0) = move(target);
	}

	PolynomialInequalityWithCoeffs::PolynomialInequalityWithCoeffs(uint32_t N, uint32_t sz,
	                                                               unique_ptr<ScaleFactor>&& scale,
	                                                               Vector<Matrix<Polynomial>>&& mat,
	                                                               Matrix<Polynomial>&& target)
	    : PolynomialInequality(N, sz, move(scale)), mat_(move(mat)), target_(move(target))
	{
		[[maybe_unused]] auto deg = int32_t(PolynomialInequality::get_scale()->max_degree());
		assert(mat_.size() == N);
		for (uint32_t i = 0; i < N; ++i)
		{
			assert(mat_[i].is_square() && mat_[i].row() == sz);
			for (uint32_t r = 0; r < sz; ++r)
				for (uint32_t c = 0; c < sz; ++c) assert(mat_[i].at(r, c).degree() <= deg);
		}
		assert(target_.is_square() && target_.row() == sz);
		for (uint32_t r = 0; r < sz; ++r)
			for (uint32_t c = 0; c < sz; ++c) assert(target_.at(r, c).degree() <= deg);
	}

	void PolynomialProgram::add_equation(Vector<real>&& vec, real&& target) &
	{
		assert(vec.size() == N_);
		// apply previous equations to new equation
		for (uint32_t e = 0; e < equation_.size(); ++e)
		{
			const auto& t = vec[leading_indices_[e]];
			target -= t * equation_targets_[e];
			vec -= mul_scalar(t, equation_[e]);
		}
		uint32_t argmax = 0;
		for (uint32_t n = 1; n < N_; ++n)
			if (mp::cmpabs(vec[argmax], vec[n]) < 0) argmax = n;
		auto max = vec[argmax];
		target /= max;
		vec /= max;
		vec[argmax] = 1;
		// apply new equation to previous equations
		for (uint32_t e = 0; e < equation_.size(); ++e)
		{
			const auto& t = equation_[e][argmax];
			equation_targets_[e] -= t * target;
			equation_[e] -= mul_scalar(t, vec);
		}
		equation_.push_back(move(vec));
		equation_targets_.push_back(move(target));
		leading_indices_.push_back(argmax);
		for (uint32_t n = 0;; ++n)
			if (free_indices_[n] == argmax)
			{
				free_indices_.erase(free_indices_.begin() + n);
				break;
			}
	}
	SDPBInput PolynomialProgram::create_input(uint32_t parallel, _event_base* event) &&
	{
		uint32_t eq_sz = uint32_t(equation_.size()), M = N_ - eq_sz;
		// y[leading_indices_[0]], ..., y[leading_indices_[eq_sz - 1]] are eliminated
		// we rearrange y[0], ..., y[N - 1] as
		//   w[0] = y[leading_indices_[0]], ..., w[eq_sz - 1] = y[leading_indices_[eq_sz - 1]], z[0], ..., z[M - 1]
		//   where M = N_ - eq_sz and z[n] = y[free_indices_[n]]
		// for e = 0, ..., eq_sz - 1,
		//   w[e] + \sum_{m = 0}^{M - 1} equation_[e][free_indices_[m]] y[free_indices_[m]] = equation_targets_[e]
		// terms from eliminated variables w[e]
		Vector<real> obj_new(M);
		for (uint32_t m = 0; m < M; ++m) obj_new[m] = move(obj_[free_indices_[m]]);
		for (uint32_t e = 0; e < eq_sz; ++e)
		{
			const auto& t = obj_[leading_indices_[e]];
			for (uint32_t m = 0; m < M; ++m) obj_new[m] -= equation_[e][free_indices_[m]] * t;
			obj_const_ += equation_targets_[e] * t;
		}
		SDPBInput sdpb(move(obj_const_), move(obj_new), uint32_t(inequality_.size()));
		TaskQueue q(parallel);
		std::vector<std::future<bool>> tasks;
		for (uint32_t j = 0; j < inequality_.size(); ++j)
			tasks.push_back(q.push([this, &sdpb, j, M, eq_sz, event]() {
				QBOOT_scope(scope, std::to_string(j), event);
				auto&& ineq = inequality_[j];
				// convert ineq to DualConstraint
				uint32_t sz = ineq->size(), deg = ineq->max_degree(), schur_sz = (deg + 1) * sz * (sz + 1) / 2,
				         d0 = deg / 2, d1 = deg == 0 ? 0 : (deg - 1) / 2;
				Vector<real> d_c(schur_sz);
				Matrix<real> d_B(schur_sz, M);
				Matrix<real> q0(d0 + 1, deg + 1);
				Matrix<real> q1(d1 + 1, deg + 1);
				const auto& xs = ineq->sample_points();
				const auto& scs = ineq->sample_scalings();
				Vector<Matrix<real>> e_c(deg + 1);
				Vector<Vector<Matrix<real>>> e_B(N_);
				for (uint32_t k = 0; k <= deg; ++k) e_c[k] = -ineq->target_eval_with_scale(k);
				for (uint32_t n = 0; n < N_; ++n)
				{
					e_B[n] = Vector<Matrix<real>>{deg + 1};
					for (uint32_t k = 0; k <= deg; ++k) e_B[n][k] = -ineq->matrix_eval_with_scale(n, k);
				}
				// Tr(A_p Y) + (e_B y)_p = (e_c)_p
				// convert to Tr(A_p Z) + (d_B z)_p = (d_c)_p
				uint32_t p = 0;
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c <= r; ++c)
						for (uint32_t k = 0; k <= deg; ++k)
						{
							d_c.at(p) = move(e_c[k].at(r, c));
							for (uint32_t m = 0; m < M; ++m) d_B.at(p, m) = move(e_B[free_indices_[m]][k].at(r, c));
							// terms from eliminated variables w[e]
							for (uint32_t e = 0; e < eq_sz; ++e)
							{
								const auto& t = e_B[leading_indices_[e]][k].at(r, c);
								for (uint32_t m = 0; m < M; ++m) d_B.at(p, m) -= equation_[e][free_indices_[m]] * t;
								d_c.at(p) -= equation_targets_[e] * t;
							}
							++p;
						}
				for (uint32_t m = 0; m <= d0; ++m)
					for (uint32_t k = 0; k <= deg; ++k)
						q0.at(m, k) = ineq->bilinear_bases()[m].eval(xs[k]) * mp::sqrt(scs[k]);
				for (uint32_t m = 0; m <= d1; ++m)
					for (uint32_t k = 0; k <= deg; ++k)
						q1.at(m, k) = ineq->bilinear_bases()[m].eval(xs[k]) * mp::sqrt(scs[k] * xs[k]);
				sdpb.register_constraint(
				    j, DualConstraint(ineq->size(), ineq->max_degree(), move(d_B), move(d_c), {move(q0), move(q1)}));
				return true;
			}));
		for (auto&& x : tasks) x.get();
		return sdpb;
	}
	XMLInput PolynomialProgram::create_xml(uint32_t parallel, _event_base* event) &&
	{
		uint32_t eq_sz = uint32_t(equation_.size()), M = N_ - eq_sz;
		Vector<real> obj_new(M);
		for (uint32_t m = 0; m < M; ++m) obj_new[m] = move(obj_[free_indices_[m]]);
		for (uint32_t e = 0; e < eq_sz; ++e)
		{
			const auto& t = obj_[leading_indices_[e]];
			for (uint32_t m = 0; m < M; ++m) obj_new[m] -= equation_[e][free_indices_[m]] * t;
			obj_const_ += equation_targets_[e] * t;
		}
		XMLInput sdpb(move(obj_const_), move(obj_new), uint32_t(inequality_.size()));
		TaskQueue q(parallel);
		std::vector<std::future<bool>> tasks;
		for (uint32_t j = 0; j < inequality_.size(); ++j)
			tasks.push_back(q.push([this, &sdpb, j, M, eq_sz, event]() {
				QBOOT_scope(scope, std::to_string(j), event);
				auto&& ineq = inequality_[j];
				uint32_t sz = ineq->size();
				Matrix<Vector<Polynomial>> mat(sz, sz);
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) mat.at(r, c) = Vector<Polynomial>(M + 1);
				auto target = -ineq->target_polynomial();
				Vector<Matrix<Polynomial>> new_mat(M);
				for (uint32_t m = 0; m < M; ++m) new_mat[m] = ineq->matrix_polynomial(free_indices_[m]);
				for (uint32_t e = 0; e < eq_sz; ++e)
				{
					auto t = ineq->matrix_polynomial(leading_indices_[e]);
					for (uint32_t m = 0; m < M; ++m) new_mat[m] -= mul_scalar(equation_[e][free_indices_[m]], t);
					target += mul_scalar(equation_targets_[e], t);
				}
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c)
					{
						mat.at(r, c).at(0) = move(target.at(r, c));
						for (uint32_t m = 0; m < M; ++m) mat.at(r, c).at(m + 1) = move(new_mat[m].at(r, c));
					}
				sdpb.register_constraint(
				    j, PVM(move(mat), ineq->sample_points(), ineq->sample_scalings(), ineq->bilinear_bases()));
				return true;
			}));
		for (auto&& x : tasks) x.get();
		return sdpb;
	}
}  // namespace qboot