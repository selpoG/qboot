#ifndef QBOOT_POLYNOMIAL_PROGRAM_HPP_
#define QBOOT_POLYNOMIAL_PROGRAM_HPP_

#include <cstdint>  // for uint32_t
#include <memory>   // for unique_ptr
#include <vector>   // for vector

#include "matrix.hpp"      // for Matrix, Vector
#include "polynomial.hpp"  // for Polynomial
#include "real.hpp"        // for real
#include "sdpb_input.hpp"  // for SDPBInput

namespace qboot
{
	// represents an inequalities:
	//   \sum_{n = 0}^{N - 1} y[n] M[n](x) >= C(x) for all x >= 0
	// where M[n] and C are polynomial (square) matrices of x
	// Let d be the max degree of M[n] and C.
	// to solve this by SDPB, we need
	//   - bilinear bases q[m](x) (for m = 0, ..., d / 2)
	//     deg(q[m](x)) = m
	//   - sample points x[k] (for k = 0, ..., d)
	//   - sample scalings s[k] (for k = 0, ..., d)
	// each element of M[n] and C are polynomial of x, but SDPB v.2 requires only there evaluation at x = x[k]
	// while SDPB v.1 requires explicit coefficients of polynomials in xml file
	// class PolynomialInequalityEvaluated holds only evaluated matrices
	// class PolynomialInequalityWithCoeffs holds full coefficients of matrices
	// class PolynomialInequality is the base abstract class of the two
	template <class Real>
	class PolynomialInequality
	{
		uint32_t N_, sz_, max_deg_;
		algebra::Vector<algebra::Polynomial<Real>> bilinear_;
		algebra::Vector<Real> sample_x, sample_sc;

	public:
		PolynomialInequality(uint32_t N, uint32_t sz, uint32_t max_deg,
		                     algebra::Vector<algebra::Polynomial<Real>>&& bilinear_bases,
		                     algebra::Vector<Real>&& sample_points, algebra::Vector<Real>&& sample_scalings)
		    : N_(N),
		      sz_(sz),
		      max_deg_(max_deg),
		      bilinear_(std::move(bilinear_bases)),
		      sample_x(std::move(sample_points)),
		      sample_sc(std::move(sample_scalings))
		{
			assert(bilinear_.size() >= max_deg / 2 + 1);
			assert(sample_x.size() == max_deg + 1);
			assert(sample_sc.size() == max_deg + 1);
			for (uint32_t i = 0; i < bilinear_.size(); ++i) assert(bilinear_[i].degree() == int32_t(i));
		}
		PolynomialInequality(const PolynomialInequality&) = delete;
		PolynomialInequality& operator=(const PolynomialInequality&) = delete;
		PolynomialInequality(PolynomialInequality&&) noexcept = default;
		PolynomialInequality& operator=(PolynomialInequality&&) noexcept = default;
		virtual ~PolynomialInequality() = default;
		// size of matrix
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }
		// size of matrix
		[[nodiscard]] uint32_t size() const { return sz_; }
		// size of matrix
		[[nodiscard]] uint32_t max_degree() const { return max_deg_; }
		// q[m](x)
		[[nodiscard]] const algebra::Vector<algebra::Polynomial<Real>>& bilinear_bases() const { return bilinear_; }
		// x[k]
		[[nodiscard]] const algebra::Vector<Real>& sample_points() const { return sample_x; }
		// s[k]
		[[nodiscard]] const algebra::Vector<Real>& sample_scalings() const { return sample_sc; }
		// M[n]
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& = 0;
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && = 0;
		// C
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& = 0;
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && = 0;
		// evaluate M[n] at x = sample_points[k]
		[[nodiscard]] virtual algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& = 0;
		[[nodiscard]] virtual algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && = 0;
		// evaluate C at x = sample_points[k]
		[[nodiscard]] virtual algebra::Matrix<Real> target_eval(uint32_t k) const& = 0;
		[[nodiscard]] virtual algebra::Matrix<Real> target_eval(uint32_t k) && = 0;
	};

	template <class Real>
	class PolynomialInequalityEvaluated : public PolynomialInequality<Real>
	{
		// mat[n]: evaluated values of M[n]
		// mat[n][k]: M[n] evaluated at x = x[k]
		algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> mat_;
		algebra::Vector<algebra::Matrix<Real>> target_;

	public:
		PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, uint32_t max_deg,
		                              algebra::Vector<algebra::Vector<algebra::Matrix<Real>>>&& mat,
		                              algebra::Vector<algebra::Matrix<Real>>&& target,
		                              algebra::Vector<algebra::Polynomial<Real>>&& bilinear_bases,
		                              algebra::Vector<Real>&& sample_points, algebra::Vector<Real>&& sample_scalings)
		    : PolynomialInequality<Real>(N, sz, max_deg, std::move(bilinear_bases), std::move(sample_points),
		                                 std::move(sample_scalings)),
		      mat_(std::move(mat)),
		      target_(std::move(target))
		{
			assert(mat_.size() == N);
			for (uint32_t i = 0; i < N; i++)
			{
				assert(mat_[i].size() == max_deg + 1);
				for (uint32_t k = 0; k <= max_deg; k++) assert(mat_[i][k].is_square() && mat_[i][k].row() == sz);
			}
			assert(target_.size() == max_deg + 1);
			for (uint32_t k = 0; k <= max_deg; k++) assert(target_[k].is_square() && target_[k].row() == sz);
		}
		PolynomialInequalityEvaluated(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated& operator=(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated(PolynomialInequalityEvaluated&&) noexcept = default;
		PolynomialInequalityEvaluated& operator=(PolynomialInequalityEvaluated&&) noexcept = default;
		~PolynomialInequalityEvaluated() override = default;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& override
		{
			return algebra::polynomial_interpolate(mat_[n], PolynomialInequality<Real>::sample_points());
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && override
		{
			return matrix_polynomial(n);
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& override
		{
			return algebra::polynomial_interpolate(target_, PolynomialInequality<Real>::sample_points());
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && override
		{
			return target_polynomial();
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& override
		{
			return mat_[n][k].clone();
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && override
		{
			return std::move(mat_[n][k]);
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval(uint32_t k) const& override { return target_[k].clone(); }
		[[nodiscard]] algebra::Matrix<Real> target_eval(uint32_t k) && override { return std::move(target_[k]); }
	};

	template <class Real>
	class PolynomialInequalityWithCoeffs : public PolynomialInequality<Real>
	{
		algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>> mat_;
		algebra::Matrix<algebra::Polynomial<Real>> target_;

	public:
		PolynomialInequalityWithCoeffs(uint32_t N, uint32_t max_deg, algebra::Vector<algebra::Polynomial<Real>>&& mat,
		                               algebra::Polynomial<Real>&& target,
		                               algebra::Vector<algebra::Polynomial<Real>>&& bilinear_bases,
		                               algebra::Vector<Real>&& sample_points, algebra::Vector<Real>&& sample_scalings)
		    : PolynomialInequality<Real>(N, 1, max_deg, std::move(bilinear_bases), std::move(sample_points),
		                                 std::move(sample_scalings)),
		      mat_{},
		      target_{}
		{
			assert(mat.size() == N);
			mat_ = algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>>(N);
			for (uint32_t i = 0; i < N; i++)
			{
				mat_[i] = {1, 1};
				assert(mat[i].degree() <= int32_t(max_deg));
				mat_[i].at(0, 0) = std::move(mat[i]);
			}
			target_ = {1, 1};
			assert(target.degree() <= int32_t(max_deg));
			target_.at(0, 0) = std::move(target);
		}
		PolynomialInequalityWithCoeffs(uint32_t N, uint32_t sz, uint32_t max_deg,
		                               algebra::Vector<algebra::Matrix<algebra::Polynomial<Real>>>&& mat,
		                               algebra::Matrix<algebra::Polynomial<Real>>&& target,
		                               algebra::Vector<algebra::Polynomial<Real>>&& bilinear_bases,
		                               algebra::Vector<Real>&& sample_points, algebra::Vector<Real>&& sample_scalings)
		    : PolynomialInequality<Real>(N, sz, max_deg, std::move(bilinear_bases), std::move(sample_points),
		                                 std::move(sample_scalings)),
		      mat_(std::move(mat)),
		      target_(std::move(target))
		{
			assert(mat_.size() == N);
			for (uint32_t i = 0; i < N; i++)
			{
				assert(mat_[i].is_square() && mat_[i].row() == sz);
				for (uint32_t r = 0; r < sz; r++)
					for (uint32_t c = 0; c < sz; c++) assert(mat_[i].at(r, c).degree() <= int32_t(max_deg));
			}
			assert(target_.is_square() && target_.row() == sz);
			for (uint32_t r = 0; r < sz; r++)
				for (uint32_t c = 0; c < sz; c++) assert(target_.at(r, c).degree() <= int32_t(max_deg));
		}
		PolynomialInequalityWithCoeffs(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs& operator=(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs(PolynomialInequalityWithCoeffs&&) noexcept = default;
		PolynomialInequalityWithCoeffs& operator=(PolynomialInequalityWithCoeffs&&) noexcept = default;
		~PolynomialInequalityWithCoeffs() override = default;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& override
		{
			return mat_[n].clone();
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && override
		{
			return std::move(mat_[n]);
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& override
		{
			return target_.clone();
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && override
		{
			return std::move(target_);
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& override
		{
			return mat_[n].eval(PolynomialInequality<Real>::sample_points()[k]);
		}
		[[nodiscard]] algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && override
		{
			return matrix_eval(n, k);
		}
		[[nodiscard]] algebra::Matrix<Real> target_eval(uint32_t k) const& override
		{
			return target_.eval(PolynomialInequality<Real>::sample_points()[k]);
		}
		algebra::Matrix<Real> target_eval(uint32_t k) && override { return target_eval(k); }
	};

	// represents a polynomial matrix programming (maximize some linear quantity subject to semidefinite positivity)
	// maximize \sum_{n=0}^{N-1} b[n] y[n] + b[N]
	// over free (real) variables y[0], ..., y[N - 1]
	// e-th equations:
	//   \sum_{n = 0}^{N - 1} y[n] M_{e}^{n} = C_{e}
	//   (each of M_{e}^{n} or C_{e} must be a real constant, which is independent from x)
	// j-th inequalities:
	//   \sum_{n = 0}^{N - 1} y[n] M_{j}^{n}(x) >= C_{j}(x) for all x >= 0
	//   (each of M_{j}^{n} or C_{j} can be a polynomial matrix of x)
	template <class Real>
	class PolynomialProgramming
	{
		// num of free variables
		uint32_t N_;
		// b[N]
		Real obj_const_{};
		// b[n]
		algebra::Vector<Real> obj_;
		// M_{e}^{n}
		std::vector<algebra::Vector<Real>> equation_{};
		// C_{e}
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
		explicit PolynomialProgramming(uint32_t num_of_vars) : N_(num_of_vars), obj_(num_of_vars)
		{
			for (uint32_t i = 0; i < N_; i++) free_indices_.push_back(i);
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
			target /= vec[argmax];
			vec /= vec[argmax];
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
			SDPBInput<Real> sdpb(obj_const_, std::move(obj_), uint32_t(inequality_.size()));
			for (uint32_t j = 0; j < inequality_.size(); j++)
			{
				auto&& ineq = inequality_[j];
				// convert ineq to DualConstraint
				uint32_t sz = ineq->size(), deg = ineq->max_degree(), schur_sz = (deg + 1) * sz * (sz + 1) / 2,
				         d0 = deg / 2, d1 = deg == 0 ? 0 : (deg - 1) / 2;
				algebra::Matrix<Real> d_B(schur_sz, M);
				algebra::Vector<Real> d_c(schur_sz);
				algebra::Matrix<Real> q0(d0 + 1, deg + 1);
				algebra::Matrix<Real> q1(d1 + 1, deg + 1);
				const auto& xs = ineq->sample_points();
				const auto& scs = ineq->sample_scalings();
				algebra::Vector<algebra::Vector<algebra::Matrix<Real>>> e_B(N_);
				algebra::Vector<algebra::Matrix<Real>> e_c(deg + 1);
				for (uint32_t k = 0; k <= deg; ++k) e_c[k] = mul_scalar(-scs[k], std::move(*ineq).target_eval(k));
				for (uint32_t n = 0; n < N_; ++n)
				{
					e_B[n] = algebra::Vector<algebra::Matrix<Real>>{deg + 1};
					for (uint32_t k = 0; k <= deg; ++k)
						e_B[n][k] = mul_scalar(-scs[k], std::move(*ineq).matrix_eval(n, k));
				}
				uint32_t p = 0;
				for (uint32_t r = 0; r < sz; ++r)
					for (uint32_t c = 0; c <= r; ++c)
						for (uint32_t k = 0; k <= deg; ++k)
						{
							d_c.at(p) = std::move(e_c[k].at(r, c));
							for (uint32_t m = 0; m < M; ++m)
								d_B.at(p, m) = std::move(e_B[free_indices_[m]][k].at(r, c));
							for (uint32_t e = 0; e < eq_sz; ++e)
							{
								for (uint32_t m = 0; m < M; m++)
									d_B.at(p, free_indices_[m]) -=
									    equation_[e][free_indices_[m]] * d_B.at(p, leading_indices_[e]);
								d_c.at(p) -= equation_targets_[e] * d_B.at(p, leading_indices_[e]);
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
