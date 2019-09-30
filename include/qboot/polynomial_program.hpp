#ifndef QBOOT_POLYNOMIAL_PROGRAM_HPP_
#define QBOOT_POLYNOMIAL_PROGRAM_HPP_

#include <cstdint>  // for uint32_t
#include <memory>   // for unique_ptr
#include <vector>   // for vector

#include "qboot/algebra/matrix.hpp"      // for Matrix, Vector
#include "qboot/algebra/polynomial.hpp"  // for Polynomial
#include "qboot/mp/real.hpp"             // for real
#include "qboot/scale_factor.hpp"        // for ScaleFactor
#include "qboot/sdpb_input.hpp"          // for SDPBInput
#include "qboot/task_queue.hpp"          // for _event_base
#include "qboot/xml_input.hpp"           // for XMLInput

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
		// mat[n]: evaluated values of M[n]
		// mat[n][k]: M[n] evaluated at x = x[k]
		algebra::Vector<algebra::Vector<algebra::Matrix<mp::real>>> mat_;
		// M[N]
		algebra::Vector<algebra::Matrix<mp::real>> target_{};

		// get M[n] / chi or M[N] / chi as a polynomial matrix
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> as_polynomial(
		    const algebra::Vector<algebra::Matrix<mp::real>>& vals);

	public:
		PolynomialInequality(uint32_t N, std::unique_ptr<ScaleFactor>&& scale,
		                     algebra::Vector<algebra::Vector<mp::real>>&& mat, algebra::Vector<mp::real>&& target);
		PolynomialInequality(uint32_t N, algebra::Vector<mp::real>&& mat, mp::real&& target);
		PolynomialInequality(uint32_t N, uint32_t sz, algebra::Vector<algebra::Matrix<mp::real>>&& mat,
		                     algebra::Matrix<mp::real>&& target);
		PolynomialInequality(uint32_t N, uint32_t sz, std::unique_ptr<ScaleFactor>&& scale,
		                     algebra::Vector<algebra::Vector<algebra::Matrix<mp::real>>>&& mat,
		                     algebra::Vector<algebra::Matrix<mp::real>>&& target);
		PolynomialInequality(const PolynomialInequality&) = delete;
		PolynomialInequality& operator=(const PolynomialInequality&) = delete;
		PolynomialInequality(PolynomialInequality&&) noexcept = default;
		PolynomialInequality& operator=(PolynomialInequality&&) noexcept = default;
		~PolynomialInequality() = default;
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
		[[nodiscard]] algebra::Vector<mp::real> sample_points() { return chi_->sample_points(); }

		// {s_0, ..., s_D}
		[[nodiscard]] algebra::Vector<mp::real> sample_scalings() { return chi_->sample_scalings(); }

		// M[n] / chi (0 <= n < N)
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> matrix_polynomial(uint32_t n)
		{
			return as_polynomial(mat_[n]);
		}

		// M[N] / chi
		[[nodiscard]] algebra::Matrix<algebra::Polynomial> target_polynomial() { return as_polynomial(target_); }

		// evaluate M[n] at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] algebra::Matrix<mp::real> matrix_eval_with_scale(uint32_t n, uint32_t k)
		{
			return mat_[n][k].clone();
		}

		// evaluate M[N] at x = x_k
		[[nodiscard]] algebra::Matrix<mp::real> target_eval_with_scale(uint32_t k) { return target_[k].clone(); }

		// evaluate M[n] / chi at x = x_k (0 <= n < N, 0 <= k <= D)
		[[nodiscard]] algebra::Matrix<mp::real> matrix_eval_without_scale(uint32_t n, uint32_t k)
		{
			const auto& chi = PolynomialInequality::get_scale();
			return mat_[n][k] / chi->eval(chi->sample_point(k));
		}

		// evaluate M[N] / chi at x = x_k
		[[nodiscard]] algebra::Matrix<mp::real> target_eval_without_scale(uint32_t k)
		{
			const auto& chi = PolynomialInequality::get_scale();
			return target_[k] / chi->eval(chi->sample_point(k));
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
		mp::real obj_const_{};
		// b[n]
		algebra::Vector<mp::real> obj_;
		// M_e[n]
		std::vector<algebra::Vector<mp::real>> equation_{};
		// M_e[N]
		std::vector<mp::real> equation_targets_{};
		// pivots for gaussian elimination
		// y[leading_indices[j]] are not free
		std::vector<uint32_t> leading_indices_{};
		// indices which are not eliminated
		// union of leading_indices_ and free_indices_ is always {0, ..., N - 1}
		std::vector<uint32_t> free_indices_{};
		std::vector<std::unique_ptr<PolynomialInequality>> inequality_{};

	public:
		explicit PolynomialProgram(uint32_t num_of_vars) : N_(num_of_vars), obj_(num_of_vars)
		{
			for (uint32_t i = 0; i < N_; ++i) free_indices_.push_back(i);
		}
		[[nodiscard]] uint32_t num_of_variables() const { return N_; }
		[[nodiscard]] mp::real& objective_constant() { return obj_const_; }
		[[nodiscard]] const mp::real& objective_constant() const { return obj_const_; }
		[[nodiscard]] const algebra::Vector<mp::real>& objectives() const { return obj_; }
		void objectives(algebra::Vector<mp::real>&& obj) &
		{
			assert(obj.size() == N_);
			obj_ = std::move(obj);
		}
		// add a constraint \sum_{n = 0}^{N - 1} y[n] mat[n] = target
		// we assume that all equations are linear independent
		// the order of call of this function may affects the resulting SDPB input
		// to guarantee the reproducibility, call this function in some fixed order
		void add_equation(algebra::Vector<mp::real>&& vec, mp::real&& target) &;
		void add_inequality(std::unique_ptr<PolynomialInequality>&& ineq) &
		{
			assert(ineq->num_of_variables() == N_);
			inequality_.push_back(std::move(ineq));
		}
		SDPBInput create_input(uint32_t parallel = 1, _event_base* event = nullptr) &&;
		XMLInput create_xml(uint32_t parallel = 1, _event_base* event = nullptr) &&;
	};
}  // namespace qboot

#endif  // QBOOT_POLYNOMIAL_PROGRAM_HPP_
