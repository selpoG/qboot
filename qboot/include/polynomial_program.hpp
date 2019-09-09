#ifndef QBOOT_POLYNOMIAL_PROGRAM_HPP_
#define QBOOT_POLYNOMIAL_PROGRAM_HPP_

#include <cstdint>  // for uint32_t
#include <memory>   // for unique_ptr
#include <vector>   // for vector

#include "matrix.hpp"        // for Matrix, Vector
#include "polynomial.hpp"    // for Polynomial
#include "real.hpp"          // for real
#include "scale_factor.hpp"  // for ScaleFactor
#include "sdpb_input.hpp"    // for SDPBInput

namespace qboot
{
	// represents an inequalities:
	//   \sum_{n = 0}^{N - 1} y[n] M[n](x) >= M[N](x) for all x >= 0
	// where M[n] are symmetric matrices of x whose elements are products of ScaleFactor chi and polynomial.
	// Let D be the max degree of M[n] / chi.
	// to solve this by SDPB, we need
	//   - bilinear bases q[m](x) (for m = 0, ..., D / 2)
	//     deg(q[m](x)) = m
	//   - sample points x[k] (for k = 0, ..., D)
	//   - sample scalings s[k] (for k = 0, ..., D)
	// each element of M[n] / chi is a polynomial of x, but SDPB v.2 requires only there evaluation at x = x[k]
	// while SDPB v.1 requires explicit coefficients of polynomials in xml file
	// class PolynomialInequalityEvaluated holds only evaluated matrices
	// class PolynomialInequalityWithCoeffs holds full coefficients of matrices
	// class PolynomialInequality is the base abstract class of the two
	template <class Real>
	class PolynomialInequality
	{
		uint32_t N_, sz_;
		std::unique_ptr<ScaleFactor<Real>> chi_;

	public:
		PolynomialInequality(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor<Real>>&& chi)
		    : N_(N), sz_(sz), chi_(std::move(chi))
		{
		}
		PolynomialInequality(const PolynomialInequality&) = delete;
		PolynomialInequality& operator=(const PolynomialInequality&) = delete;
		PolynomialInequality(PolynomialInequality&&) noexcept = default;
		PolynomialInequality& operator=(PolynomialInequality&&) noexcept = default;
		virtual ~PolynomialInequality() = default;
		[[nodiscard]] const std::unique_ptr<ScaleFactor<Real>>& get_scale() const& { return chi_; }
		[[nodiscard]] std::unique_ptr<ScaleFactor<Real>> get_scale() && { return std::move(chi_); }
		// number of free variables N
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }

		// size of matrix
		[[nodiscard]] uint32_t size() const { return sz_; }

		// max degree D among all elements in M[n] / chi
		[[nodiscard]] uint32_t max_degree() const { return chi_->max_degree(); }

		// {q[0](x), ..., q[D / 2](x)}
		[[nodiscard]] algebra::Vector<algebra::Polynomial<Real>> bilinear_bases() { return chi_->bilinear_bases(); }

		// {x_0, ..., x_D}
		[[nodiscard]] algebra::Vector<Real> sample_points() { return chi_->sample_points(); }

		// {s_0, ..., s_D}
		[[nodiscard]] algebra::Vector<Real> sample_scalings() { return chi_->sample_scalings(); }

		// M[n] / chi (0 <= n < N)
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) = 0;

		// M[N] / chi
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() = 0;

		// evaluate M[n] at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] virtual algebra::Matrix<Real> matrix_eval_with_scale(uint32_t n, uint32_t k) = 0;

		// evaluate M[N] at x = x_k
		[[nodiscard]] virtual algebra::Matrix<Real> target_eval_with_scale(uint32_t k) = 0;

		// evaluate M[n] / chi at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] virtual algebra::Matrix<Real> matrix_eval_without_scale(uint32_t n, uint32_t k) = 0;

		// evaluate M[N] / chi at x = x_k
		[[nodiscard]] virtual algebra::Matrix<Real> target_eval_without_scale(uint32_t k) = 0;
	};

	template <class Real>
	class PolynomialInequalityEvaluated : public PolynomialInequality<Real>
	{
		// mat[n]: evaluated values of M[n]
		// mat[n][k]: M[n] evaluated at x = x[k]
		algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> mat_;
		// M[N]
		algebra::Vector<algebra::Matrix<Real>> target_;

		// get M[n] / chi or M[N] / chi as a polynomial matrix
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> as_polynomial(
		    const algebra::Vector<algebra::Matrix<Real>>& vals)
		{
			uint32_t deg = vals.size() - 1;
			algebra::Vector<algebra::Matrix<Real>> ev(deg + 1);
			const auto& chi = PolynomialInequality<Real>::get_scale();
			auto xs = chi->sample_points();
			for (uint32_t k = 0; k <= deg; ++k) ev[k] = vals[k] / chi->eval(xs[k]);
			return algebra::polynomial_interpolate(ev, xs);
		}

	public:
		PolynomialInequalityEvaluated(uint32_t N, std::unique_ptr<ScaleFactor<Real>>&& scale,
		                              algebra::Vector<algebra::Vector<Real>>&& mat, algebra::Vector<Real>&& target)
		    : PolynomialInequality<Real>(N, 1, std::move(scale)), mat_(N), target_{}
		{
			uint32_t deg = PolynomialInequality<Real>::get_scale()->max_degree();
			for (uint32_t i = 0; i < N; ++i)
			{
				assert(mat[i].size() == deg + 1);
				mat_[i] = algebra::Vector<algebra::Matrix<Real>>{deg + 1};
				for (uint32_t k = 0; k <= deg; ++k)
				{
					mat_[i][k] = {1, 1};
					mat_[i][k].at(0, 0) = std::move(mat[i][k]);
				}
			}
			assert(target.size() == deg + 1);
			target_ = algebra::Vector<algebra::Matrix<Real>>{deg + 1};
			for (uint32_t k = 0; k <= deg; ++k)
			{
				target_[k] = {1, 1};
				target_[k].at(0, 0) = std::move(target[k]);
			}
		}
		PolynomialInequalityEvaluated(uint32_t N, algebra::Vector<Real>&& mat, Real&& target)
		    : PolynomialInequality<Real>(N, 1, std::make_unique<TrivialScale<Real>>()), mat_(N), target_(1)
		{
			assert(mat.size() == N);
			for (uint32_t i = 0; i < N; ++i)
			{
				mat_[i] = algebra::Vector<algebra::Matrix<Real>>(1);
				mat_[i][0] = {1, 1};
				mat_[i][0].at(0, 0) = std::move(mat[i]);
			}
			target_[0] = {1, 1};
			target_[0].at(0, 0) = std::move(target);
		}
		PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, algebra::Vector<algebra::Matrix<Real>>&& mat,
		                              algebra::Matrix<Real>&& target)
		    : PolynomialInequality<Real>(N, sz, std::make_unique<TrivialScale<Real>>()), mat_(N), target_(1)
		{
			assert(mat.size() == N);
			for (uint32_t i = 0; i < N; ++i)
			{
				assert(mat[i].is_square() && mat[i].row() == sz);
				mat_[i] = algebra::Vector<algebra::Matrix<Real>>(1);
				mat_[i][0] = {sz, sz};
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c < sz; ++c) mat_[i][0].at(r, c) = std::move(mat[i].at(r, c));
			}
			assert(target.is_square() && target.row() == sz);
			target_[0] = {sz, sz};
			for (uint32_t r = 0; r < sz; ++r)
				for (uint32_t c = 0; c < sz; ++c) target_[0].at(r, c) = std::move(target.at(r, c));
		}
		PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor<Real>>&& scale,
		                              algebra::Vector<algebra::Vector<algebra::Matrix<Real>>>&& mat,
		                              algebra::Vector<algebra::Matrix<Real>>&& target)
		    : PolynomialInequality<Real>(N, sz, std::move(scale)), mat_(std::move(mat)), target_(std::move(target))
		{
			uint32_t deg = PolynomialInequality<Real>::get_scale()->max_degree();
			assert(mat_.size() == N);
			for (uint32_t i = 0; i < N; ++i)
			{
				assert(mat_[i].size() == deg + 1);
				for (uint32_t k = 0; k <= deg; ++k) assert(mat_[i][k].is_square() && mat_[i][k].row() == sz);
			}
			assert(target_.size() == deg + 1);
			for (uint32_t k = 0; k <= deg; ++k) assert(target_[k].is_square() && target_[k].row() == sz);
		}
		PolynomialInequalityEvaluated(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated& operator=(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated(PolynomialInequalityEvaluated&&) noexcept = default;
		PolynomialInequalityEvaluated& operator=(PolynomialInequalityEvaluated&&) noexcept = default;
		~PolynomialInequalityEvaluated() override = default;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) override
		{
			return as_polynomial(mat_[n]);
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() override
		{
			return as_polynomial(target_);
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval_with_scale(uint32_t n, uint32_t k) override
		{
			return mat_[n][k].clone();
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval_with_scale(uint32_t k) override { return target_[k].clone(); }
		[[nodiscard]] algebra::Matrix<Real> matrix_eval_without_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			return mat_[n][k] / chi->eval(chi->sample_point(k));
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval_without_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			return target_[k] / chi->eval(chi->sample_point(k));
		}
	};

	template <class Real>
	class PolynomialInequalityWithCoeffs : public PolynomialInequality<Real>
	{
		// M[n] / chi
		algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>> mat_;
		// M[N] / chi
		algebra::Matrix<algebra::Polynomial<Real>> target_;

	public:
		PolynomialInequalityWithCoeffs(uint32_t N, std::unique_ptr<ScaleFactor<Real>>&& scale,
		                               algebra::Vector<algebra::Polynomial<Real>>&& mat,
		                               algebra::Polynomial<Real>&& target)
		    : PolynomialInequality<Real>(N, 1, std::move(scale)), mat_{}, target_{}
		{
			auto deg = int32_t(PolynomialInequality<Real>::get_scale()->max_degree());
			assert(mat.size() == N);
			mat_ = algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>>(N);
			for (uint32_t i = 0; i < N; ++i)
			{
				mat_[i] = {1, 1};
				assert(mat[i].degree() <= deg);
				mat_[i].at(0, 0) = std::move(mat[i]);
			}
			target_ = {1, 1};
			assert(target.degree() <= deg);
			target_.at(0, 0) = std::move(target);
		}
		PolynomialInequalityWithCoeffs(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor<Real>>&& scale,
		                               algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>>&& mat,
		                               algebra::Matrix<algebra::Polynomial<Real>>&& target)
		    : PolynomialInequality<Real>(N, sz, std::move(scale)), mat_(std::move(mat)), target_(std::move(target))
		{
			auto deg = int32_t(PolynomialInequality<Real>::get_scale()->max_degree());
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
		PolynomialInequalityWithCoeffs(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs& operator=(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs(PolynomialInequalityWithCoeffs&&) noexcept = default;
		PolynomialInequalityWithCoeffs& operator=(PolynomialInequalityWithCoeffs&&) noexcept = default;
		~PolynomialInequalityWithCoeffs() override = default;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) override
		{
			return mat_[n].clone();
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() override
		{
			return target_.clone();
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval_with_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			auto x = chi->sample_point(k);
			return mul_scalar(chi->eval(x), mat_[n].eval(x));
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval_with_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			auto x = chi->sample_point(k);
			return mul_scalar(chi->eval(x), target_.eval(x));
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval_without_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			return mat_[n].eval(chi->sample_point(k));
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval_without_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality<Real>::get_scale();
			return target_.eval(chi->sample_point(k));
		}
	};

	// represents a polynomial matrix programming (maximize some linear quantity subject to semidefinite positivity)
	// maximize \sum_{n=0}^{N-1} b[n] y[n] + b[N]
	// over free (real) variables y[0], ..., y[N - 1]
	// e-th equations:
	//   \sum_{n = 0}^{N - 1} y[n] M_e[n] = M_e[N]
	//   (each of M_e[n] must be a real constant, which is independent from x)
	// j-th inequalities:
	//   \sum_{n = 0}^{N - 1} y[n] M_j[n](x) >= M_j[N](x) for all x >= 0
	//   (each of M_j[n] can be a polynomial matrix of x)
	template <class Real>
	class PolynomialProgram
	{
		// num of free variables
		uint32_t N_;
		// b[N]
		Real obj_const_{};
		// b[n]
		algebra::Vector<Real> obj_;
		// M_e[n]
		std::vector<algebra::Vector<Real>> equation_{};
		// M_e[N]
		std::vector<Real> equation_targets_{};
		// pivots for gaussian elimination
		// y[leading_indices[j]] are not free
		std::vector<uint32_t> leading_indices_{};
		// indices which are not eliminated
		// union of leading_indices_ and free_indices_ is always {0, ..., N - 1}
		std::vector<uint32_t> free_indices_{};
		// TODO(selpo): unique_ptr is the best choice?
		std::vector<std::unique_ptr<PolynomialInequality<Real>>> inequality_{};

	public:
		explicit PolynomialProgram(uint32_t num_of_vars) : N_(num_of_vars), obj_(num_of_vars)
		{
			for (uint32_t i = 0; i < N_; ++i) free_indices_.push_back(i);
		}
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }
		[[nodiscard]] Real& objective_constant() { return obj_const_; }
		[[nodiscard]] const Real& objective_constant() const { return obj_const_; }
		[[nodiscard]] const algebra::Vector<Real>& objectives() const { return obj_; }
		void objectives(algebra::Vector<Real>&& obj) &
		{
			assert(obj.size() == N_);
			obj_ = std::move(obj);
		}
		// add a constraint \sum_{n = 0}^{N - 1} y[n] mat[n] = target
		// we assume that all equations are linear independent
		// the order of call of this function may affects the resulting SDPB input
		// to guarantee the reproducibility, call this function in some fixed order
		void add_equation(algebra::Vector<Real>&& vec, Real&& target) &
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
				if (mpfr::cmpabs(vec[argmax], vec[n]) < 0) argmax = n;
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
			equation_.push_back(std::move(vec));
			equation_targets_.push_back(std::move(target));
			leading_indices_.push_back(argmax);
			for (uint32_t n = 0;; ++n)
				if (free_indices_[n] == argmax)
				{
					free_indices_.erase(free_indices_.begin() + n);
					break;
				}
		}
		void add_inequality(std::unique_ptr<PolynomialInequality<Real>>&& ineq) &
		{
			assert(ineq->num_of_variables() == N_);
			inequality_.push_back(std::move(ineq));
		}
		SDPBInput<Real> create_input() &&
		{
			uint32_t eq_sz = uint32_t(equation_.size()), M = N_ - eq_sz;
			// y[leading_indices_[0]], ..., y[leading_indices_[eq_sz - 1]] are eliminated
			// we rearrange y[0], ..., y[N - 1] as
			//   w[0] = y[leading_indices_[0]], ..., w[eq_sz - 1] = y[leading_indices_[eq_sz - 1]], z[0], ..., z[M - 1]
			//   where M = N_ - eq_sz and z[n] = y[free_indices_[n]]
			// for e = 0, ..., eq_sz - 1,
			//   w[e] + \sum_{m = 0}^{M - 1} equation_[e][free_indices_[m]] y[free_indices_[m]] = equation_targets_[e]
			// terms from eliminated variables w[e]
			algebra::Vector<Real> obj_new(M);
			for (uint32_t m = 0; m < M; ++m) obj_new[m] = std::move(obj_[free_indices_[m]]);
			for (uint32_t e = 0; e < eq_sz; ++e)
			{
				const auto& t = obj_[leading_indices_[e]];
				for (uint32_t m = 0; m < M; ++m) obj_new[m] -= equation_[e][free_indices_[m]] * t;
				obj_const_ += equation_targets_[e] * t;
			}
			SDPBInput<Real> sdpb(std::move(obj_const_), std::move(obj_new), uint32_t(inequality_.size()));
			for (uint32_t j = 0; j < inequality_.size(); ++j)
			{
				auto&& ineq = inequality_[j];
				// convert ineq to DualConstraint
				uint32_t sz = ineq->size(), deg = ineq->max_degree(), schur_sz = (deg + 1) * sz * (sz + 1) / 2,
				         d0 = deg / 2, d1 = deg == 0 ? 0 : (deg - 1) / 2;
				algebra::Vector<Real> d_c(schur_sz);
				algebra::Matrix<Real> d_B(schur_sz, M);
				algebra::Matrix<Real> q0(d0 + 1, deg + 1);
				algebra::Matrix<Real> q1(d1 + 1, deg + 1);
				const auto& xs = ineq->sample_points();
				const auto& scs = ineq->sample_scalings();
				algebra::Vector<algebra::Matrix<Real>> e_c(deg + 1);
				algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> e_B(N_);
				for (uint32_t k = 0; k <= deg; ++k) e_c[k] = mul_scalar(-scs[k], ineq->target_eval_without_scale(k));
				for (uint32_t n = 0; n < N_; ++n)
				{
					e_B[n] = algebra::Vector<algebra::Matrix<Real>>{deg + 1};
					for (uint32_t k = 0; k <= deg; ++k)
						e_B[n][k] = mul_scalar(-scs[k], ineq->matrix_eval_without_scale(n, k));
				}
				// Tr(A_p Y) + (e_B y)_p = (e_c)_p
				// convert to Tr(A_p Z) + (d_B z)_p = (d_c)_p
				uint32_t p = 0;
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c <= r; ++c)
						for (uint32_t k = 0; k <= deg; ++k)
						{
							d_c.at(p) = std::move(e_c[k].at(r, c));
							for (uint32_t m = 0; m < M; ++m)
								d_B.at(p, m) = std::move(e_B[free_indices_[m]][k].at(r, c));
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
						q0.at(m, k) = ineq->bilinear_bases()[m].eval(xs[k]) * mpfr::sqrt(scs[k]);
				for (uint32_t m = 0; m <= d1; ++m)
					for (uint32_t k = 0; k <= deg; ++k)
						q1.at(m, k) = ineq->bilinear_bases()[m].eval(xs[k]) * mpfr::sqrt(scs[k] * xs[k]);
				sdpb.register_constraint(j, DualConstraint<Real>(ineq->size(), ineq->max_degree(), std::move(d_B),
				                                                 std::move(d_c), {std::move(q0), std::move(q1)}));
			}
			return sdpb;
		}
	};
}  // namespace qboot

#endif  // QBOOT_POLYNOMIAL_PROGRAM_HPP_
