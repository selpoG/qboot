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
	class PolynomialInequality
	{
		uint32_t N_, sz_;
		std::unique_ptr<ScaleFactor> chi_;

	public:
		PolynomialInequality(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor>&& chi)
		    : N_(N), sz_(sz), chi_(std::move(chi))
		{
		}
		PolynomialInequality(const PolynomialInequality&) = delete;
		PolynomialInequality& operator=(const PolynomialInequality&) = delete;
		PolynomialInequality(PolynomialInequality&&) noexcept = default;
		PolynomialInequality& operator=(PolynomialInequality&&) noexcept = default;
		virtual ~PolynomialInequality();
		[[nodiscard]] const std::unique_ptr<ScaleFactor>& get_scale() const& { return chi_; }
		[[nodiscard]] std::unique_ptr<ScaleFactor> get_scale() && { return std::move(chi_); }
		// number of free variables N
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }

		// size of matrix
		[[nodiscard]] uint32_t size() const { return sz_; }

		// max degree D among all elements in M[n] / chi
		[[nodiscard]] uint32_t max_degree() const { return chi_->max_degree(); }

		// {q[0](x), ..., q[D / 2](x)}
		[[nodiscard]] algebra::Vector<algebra::Polynomial> bilinear_bases() { return chi_->bilinear_bases(); }

		// {x_0, ..., x_D}
		[[nodiscard]] algebra::Vector<mpfr::real> sample_points() { return chi_->sample_points(); }

		// {s_0, ..., s_D}
		[[nodiscard]] algebra::Vector<mpfr::real> sample_scalings() { return chi_->sample_scalings(); }

		// M[n] / chi (0 <= n < N)
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial> matrix_polynomial(uint32_t n) = 0;

		// M[N] / chi
		[[nodiscard]] virtual algebra::Matrix<algebra::Polynomial> target_polynomial() = 0;

		// evaluate M[n] at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] virtual algebra::Matrix<mpfr::real> matrix_eval_with_scale(uint32_t n, uint32_t k) = 0;

		// evaluate M[N] at x = x_k
		[[nodiscard]] virtual algebra::Matrix<mpfr::real> target_eval_with_scale(uint32_t k) = 0;

		// evaluate M[n] / chi at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] virtual algebra::Matrix<mpfr::real> matrix_eval_without_scale(uint32_t n, uint32_t k) = 0;

		// evaluate M[N] / chi at x = x_k
		[[nodiscard]] virtual algebra::Matrix<mpfr::real> target_eval_without_scale(uint32_t k) = 0;
	};

	class PolynomialInequalityEvaluated : public PolynomialInequality
	{
		// mat[n]: evaluated values of M[n]
		// mat[n][k]: M[n] evaluated at x = x[k]
		algebra::Vector<algebra::Vector<algebra::Matrix<mpfr::real>>> mat_;
		// M[N]
		algebra::Vector<algebra::Matrix<mpfr::real>> target_{};

		// get M[n] / chi or M[N] / chi as a polynomial matrix
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> as_polynomial(
		    const algebra::Vector<algebra::Matrix<mpfr::real>>& vals);

	public:
		PolynomialInequalityEvaluated(uint32_t N, std::unique_ptr<ScaleFactor>&& scale,
		                              algebra::Vector<algebra::Vector<mpfr::real>>&& mat,
		                              algebra::Vector<mpfr::real>&& target);
		PolynomialInequalityEvaluated(uint32_t N, algebra::Vector<mpfr::real>&& mat, mpfr::real&& target);
		PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, algebra::Vector<algebra::Matrix<mpfr::real>>&& mat,
		                              algebra::Matrix<mpfr::real>&& target);
		PolynomialInequalityEvaluated(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor>&& scale,
		                              algebra::Vector<algebra::Vector<algebra::Matrix<mpfr::real>>>&& mat,
		                              algebra::Vector<algebra::Matrix<mpfr::real>>&& target);
		PolynomialInequalityEvaluated(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated& operator=(const PolynomialInequalityEvaluated&) = delete;
		PolynomialInequalityEvaluated(PolynomialInequalityEvaluated&&) noexcept = default;
		PolynomialInequalityEvaluated& operator=(PolynomialInequalityEvaluated&&) noexcept = default;
		~PolynomialInequalityEvaluated() override;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> matrix_polynomial(uint32_t n) override
		{
			return as_polynomial(mat_[n]);
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> target_polynomial() override
		{
			return as_polynomial(target_);
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> matrix_eval_with_scale(uint32_t n, uint32_t k) override
		{
			return mat_[n][k].clone();
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> target_eval_with_scale(uint32_t k) override
		{
			return target_[k].clone();
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> matrix_eval_without_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
			return mat_[n][k] / chi->eval(chi->sample_point(k));
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> target_eval_without_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
			return target_[k] / chi->eval(chi->sample_point(k));
		}
	};

	class PolynomialInequalityWithCoeffs : public PolynomialInequality
	{
		// M[n] / chi
		algebra::Vector<algebra::Matrix<algebra::Polynomial>> mat_{};
		// M[N] / chi
		algebra::Matrix<algebra::Polynomial> target_{};

	public:
		PolynomialInequalityWithCoeffs(uint32_t N, std::unique_ptr<ScaleFactor>&& scale,
		                               algebra::Vector<algebra::Polynomial>&& mat, algebra::Polynomial&& target);
		PolynomialInequalityWithCoeffs(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor>&& scale,
		                               algebra::Vector<algebra::Matrix<algebra::Polynomial>>&& mat,
		                               algebra::Matrix<algebra::Polynomial>&& target);
		PolynomialInequalityWithCoeffs(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs& operator=(const PolynomialInequalityWithCoeffs&) = delete;
		PolynomialInequalityWithCoeffs(PolynomialInequalityWithCoeffs&&) noexcept = default;
		PolynomialInequalityWithCoeffs& operator=(PolynomialInequalityWithCoeffs&&) noexcept = default;
		~PolynomialInequalityWithCoeffs() override;
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> matrix_polynomial(uint32_t n) override
		{
			return mat_[n].clone();
		}
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> target_polynomial() override { return target_.clone(); }
		[[nodiscard]] algebra::Matrix<mpfr::real> matrix_eval_with_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
			auto x = chi->sample_point(k);
			return mul_scalar(chi->eval(x), mat_[n].eval(x));
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> target_eval_with_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
			auto x = chi->sample_point(k);
			return mul_scalar(chi->eval(x), target_.eval(x));
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> matrix_eval_without_scale(uint32_t n, uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
			return mat_[n].eval(chi->sample_point(k));
		}
		[[nodiscard]] algebra::Matrix<mpfr::real> target_eval_without_scale(uint32_t k) override
		{
			const auto& chi = PolynomialInequality::get_scale();
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
	class PolynomialProgram
	{
		// num of free variables
		uint32_t N_;
		// b[N]
		mpfr::real obj_const_{};
		// b[n]
		algebra::Vector<mpfr::real> obj_;
		// M_e[n]
		std::vector<algebra::Vector<mpfr::real>> equation_{};
		// M_e[N]
		std::vector<mpfr::real> equation_targets_{};
		// pivots for gaussian elimination
		// y[leading_indices[j]] are not free
		std::vector<uint32_t> leading_indices_{};
		// indices which are not eliminated
		// union of leading_indices_ and free_indices_ is always {0, ..., N - 1}
		std::vector<uint32_t> free_indices_{};
		// TODO(selpo): unique_ptr is the best choice?
		std::vector<std::unique_ptr<PolynomialInequality>> inequality_{};

	public:
		explicit PolynomialProgram(uint32_t num_of_vars) : N_(num_of_vars), obj_(num_of_vars)
		{
			for (uint32_t i = 0; i < N_; ++i) free_indices_.push_back(i);
		}
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }
		[[nodiscard]] mpfr::real& objective_constant() { return obj_const_; }
		[[nodiscard]] const mpfr::real& objective_constant() const { return obj_const_; }
		[[nodiscard]] const algebra::Vector<mpfr::real>& objectives() const { return obj_; }
		void objectives(algebra::Vector<mpfr::real>&& obj) &
		{
			assert(obj.size() == N_);
			obj_ = std::move(obj);
		}
		// add a constraint \sum_{n = 0}^{N - 1} y[n] mat[n] = target
		// we assume that all equations are linear independent
		// the order of call of this function may affects the resulting SDPB input
		// to guarantee the reproducibility, call this function in some fixed order
		void add_equation(algebra::Vector<mpfr::real>&& vec, mpfr::real&& target) &;
		void add_inequality(std::unique_ptr<PolynomialInequality>&& ineq) &
		{
			assert(ineq->num_of_variables() == N_);
			inequality_.push_back(std::move(ineq));
		}
		SDPBInput create_input() &&;
	};
}  // namespace qboot

#endif  // QBOOT_POLYNOMIAL_PROGRAM_HPP_
