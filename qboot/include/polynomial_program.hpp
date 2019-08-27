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
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
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
		PolynomialInequality(PolynomialInequality&&) = default;
		PolynomialInequality& operator=(PolynomialInequality&&) = default;
		virtual ~PolynomialInequality() = default;
		// size of matrix
		uint32_t num_of_variables() const { return N_; }
		// size of matrix
		uint32_t size() const { return sz_; }
		// size of matrix
		uint32_t max_degree() const { return max_deg_; }
		// q[m](x)
		const algebra::Vector<algebra::Polynomial<Real>>& bilinear_bases() const { return bilinear_; }
		// x[k]
		const algebra::Vector<Real>& sample_points() const { return sample_x; }
		// s[k]
		const algebra::Vector<Real>& sample_scalings() const { return sample_sc; }
		// M[n]
		virtual algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& = 0;
		virtual algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && = 0;
		// C
		virtual algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& = 0;
		virtual algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && = 0;
		// evaluate M[n] at x = sample_points[k]
		virtual algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& = 0;
		virtual algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && = 0;
		// evaluate C at x = sample_points[k]
		virtual algebra::Matrix<Real> target_eval(uint32_t k) const& = 0;
		virtual algebra::Matrix<Real> target_eval(uint32_t k) && = 0;
	};

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
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
		      target_(std::move(target_))
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
		PolynomialInequalityEvaluated(PolynomialInequalityEvaluated&&) = default;
		PolynomialInequalityEvaluated& operator=(PolynomialInequalityEvaluated&&) = default;
		~PolynomialInequalityEvaluated() = default;
		algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& override
		{
			return algebra::polynomial_interpolate(mat_[n], PolynomialInequality<Real>::sample_points());
		}
		algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && override
		{
			return matrix_polynomial(n);
		}
		algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& override
		{
			return algebra::polynomial_interpolate(target_, PolynomialInequality<Real>::sample_points());
		}
		algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && override { return target_polynomial(); }
		algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& override { return mat_[n][k].clone(); }
		algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && override { return std::move(mat_[n][k]); }
		algebra::Matrix<Real> target_eval(uint32_t k) const& override { return target_[k].clone(); }
		algebra::Matrix<Real> target_eval(uint32_t k) && override { return std::move(target_[k]); }
	};

	template <class Real = mpfr::real<1000, MPFR_RNDN>>
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
		PolynomialInequalityWithCoeffs(PolynomialInequalityWithCoeffs&&) = default;
		PolynomialInequalityWithCoeffs& operator=(PolynomialInequalityWithCoeffs&&) = default;
		~PolynomialInequalityWithCoeffs() = default;
		algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) const& override
		{
			return mat_[n].clone();
		}
		algebra::Matrix<algebra::Polynomial<Real>> matrix_polynomial(uint32_t n) && override
		{
			return std::move(mat_[n]);
		}
		algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() const& override { return target_.clone(); }
		algebra::Matrix<algebra::Polynomial<Real>> target_polynomial() && override { return std::move(target_); }
		algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) const& override
		{
			return mat_[n].eval(PolynomialInequality<Real>::sample_points()[k]);
		}
		algebra::Matrix<Real> matrix_eval(uint32_t n, uint32_t k) && override { return matrix_eval(n, k); }
		algebra::Matrix<Real> target_eval(uint32_t k) const& override
		{
			return target_.eval(PolynomialInequality<Real>::sample_points()[k]);
		}
		algebra::Matrix<Real> target_eval(uint32_t k) && override { return target_eval(k); }
	};

	// represents a polynomial matrix programming (maximize some linear quantity subject to semidefinite positivity)
	// maximize \sum_{n=0}^{N-1} b[n] y[n] + b[N]
	// over free (real) variables y[0], ..., y[N - 1]
	// j-th equalities:
	//   \sum_{n = 0}^{N - 1} y[n] M_{j}^{n} = C_{j}
	//   (each of M_{j}^{n} or C_{j} must be a real constant, which is independent from x)
	// j-th inequalities:
	//   \sum_{n = 0}^{N - 1} y[n] M_{j}^{n}(x) >= C_{j}(x) for all x >= 0
	//   (each of M_{j}^{n} or C_{j} can be a polynomial matrix of x)
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class PolynomialProgramming
	{
		// num of free variables
		uint32_t N_;
		// b[N]
		Real obj_const_{};
		// b[n]
		algebra::Vector<Real> obj_;
		// M_{j}^{n}
		std::vector<algebra::Vector<Real>> equality_{};
		// C_{j}
		std::vector<Real> equality_targets_{};
		// TODO: unique_ptr is the best choice?
		std::vector<std::unique_ptr<PolynomialInequality<Real>>> inequality_{};

	public:
		PolynomialProgramming(uint32_t num_of_vars) : N_(num_of_vars), obj_(num_of_vars) {}
		uint32_t num_of_variables() const { return N_; }
		Real& objective_constant() { return obj_const_; }
		const Real& objective_constant() const { return obj_const_; }
		const algebra::Vector<Real>& objectives() const { return obj_; }
		void objectives(algebra::Vector<Real>&& obj)
		{
			assert(obj.size() == N_);
			obj_ = std::move(obj);
		}
		// add a constraint \sum_{n = 0}^{N - 1} y[n] mat[n] = target
		void add_equality(algebra::Vector<Real>&& vec, const Real& target)
		{
			assert(vec.size() == N_);
			equality_.push_back(vec);
			equality_targets_.push_back(target);
		}
		void add_inequality(std::unique_ptr<PolynomialInequality<Real>>&& ineq)
		{
			assert(ineq->num_of_variables() == N_);
			inequality_.push_back(std::move(ineq));
		}
		// TODO: implement this function
		SDPBInput<Real> create_input() const { throw 0; }
	};
}  // namespace qboot

#endif  // QBOOT_POLYNOMIAL_PROGRAM_HPP_
