#include "conformal_scale.hpp"

#include <utility>  // for move

using algebra::Vector, algebra::Polynomial, algebra::Matrix;
using mp::real, mp::log, mp::gamma_inc;
using std::move;

inline static Matrix<real> anti_band_to_inverse(const Vector<real>& ab)
{
	auto dim = (1 + ab.size()) / 2;
	Matrix<real> A(dim, dim);
	for (uint32_t i = 0; i < dim; ++i)
		for (uint32_t j = 0; j < dim; ++j) A.at(i, j) = ab[i + j];
	return algebra::lower_triangular_inverse(algebra::cholesky_decomposition(A));
}

// let pole_position = -p, base = exp(-k)
// calculate \int_{0}^{\infty} e ^ {-k x} x ^ n / (x + p) dx, for n = 0, ..., pole_order_max
// this integral equals to
// n! p ^ n e ^ {p k} \Gamma(-n, p k)
// = (-p) ^ n e ^ {p k} \Gamma(0, p k)
//   + (1 / k ^ n) \sum_{i = 0}^{n - 1} (n - i - 1)! (-p k) ^ i
// incomplete_gamma_factor = e ^ {p k} \Gamma(0, p k)
inline static Vector<real> simple_pole_integral(uint32_t pole_order_max, const real& base, const real& pole_position)
{
	real incomplete_gamma = pole_position == 0
	                            ? real(mp::global_prec)
	                            : mp::pow(base, pole_position) * gamma_inc(real(0), pole_position * log(base));
	Vector<real> result(pole_order_max + 1);
	result[0] = incomplete_gamma;
	real tmp{}, pow = pole_position * incomplete_gamma;
	real minus_pole_position = -pole_position;
	real factorial(1);
	real minus_log_base = -1 / log(base);
	real log_base_power = minus_log_base;
	for (uint32_t j = 1; j <= pole_order_max; ++j)
	{
		// factorial == (j - 1)!;
		// pow == incomplete_gamma_factor * pow(pole_position, j);
		// log_base_power == pow(minus_log_base, j);
		tmp = factorial * log_base_power + tmp * pole_position;
		// tmp == sum((j - k - 1)! * pow(pole_position, k) * pow(minus_log_base, j - k), 0 <= k < j);
		result[j] = tmp + pow;
		// result[j] == sum((j - k - 1)! * pow(pole_position, k) * pow(minus_log_base, j - k), 0 <= k < j)
		//              + incomplete_gamma_factor * pow(pole_position, j);

		if (j < pole_order_max)
		{
			pow *= pole_position;
			log_base_power *= minus_log_base;
			factorial *= j;
		}
	}
	return result;
}

inline static Vector<real> fast_partial_fraction(const Vector<real>& poles)
{
	Vector<real> result(poles.size());
	for (uint32_t i = 0; i < poles.size(); ++i)
	{
		result[i] = 1;
		for (uint32_t j = 0; j < poles.size(); ++j)
			if (i != j) result[i] *= poles[i] - poles[j];
		result[i] = 1 / result[i];
	}
	return result;
}

namespace qboot
{
	ConformalScale::~ConformalScale() = default;
	void ConformalScale::_set_bilinear_bases() &
	{
		if (bilinear_bases_.has_value()) return;
		if (end_.has_value())
		{
			bilinear_bases_ = Vector<Polynomial>{max_degree() / 2 + 1};
			for (uint32_t i = 0; i < bilinear_bases_.value().size(); ++i) bilinear_bases_.value().at(i) = Polynomial(i);
		}
		else
		{
			// orthogonal polynomial of weight function (4 rho) ^ {Delta} / \prod_i (Delta - poles[i])
			Vector<real> shifted_poles(poles_.size());
			for (uint32_t i = 0; i < poles_.size(); ++i) shifted_poles[i] = poles_[i] - gap_;
			auto weight = fast_partial_fraction(shifted_poles);
			auto deg = max_degree() / 2;
			// inner_prods[i] = \int_{0}^{\infty} dx (4 rho) ^ x x ^ i / \prod_i (x - poles[i])
			Vector<real> inner_prods(2 * deg + 1);
			for (uint32_t i = 0; i < poles_.size(); ++i)
				inner_prods += mul_scalar(weight[i], simple_pole_integral(2 * deg, 4 * *rho_, shifted_poles[i]));
			inner_prods *= mp::pow(4 * *rho_, gap_);
			auto mat = anti_band_to_inverse(inner_prods);
			bilinear_bases_ = Vector<Polynomial>{deg + 1};
			for (uint32_t i = 0; i <= deg; ++i)
			{
				Vector<real> v(i + 1);
				for (uint32_t j = 0; j <= i; ++j) v[j] = mat.at(i, j);
				bilinear_bases_.value().at(i) = Polynomial(move(v));
			}
		}
	}
	ConformalScale::ConformalScale(uint32_t cutoff, uint32_t spin, const Context& context, bool include_odd)
	    : odd_included_(include_odd),
	      spin_(spin),
	      lambda_(context.lambda()),
	      epsilon_(&context.epsilon()),
	      rho_(&context.rho()),
	      poles_(cutoff)
	{
		// type 1 or 3 PoleData vanishes when delta12 == 0 or delta34 == 0 and k is odd
		auto get_pols = [this, include_odd](uint32_t type) {
			return PoleSequence(type, this->spin_, *this->epsilon_, include_odd);
		};
		auto pole_seq = Merged(Merged(get_pols(1), get_pols(2)), get_pols(3));
		uint32_t pos = 0;
		while (pos < cutoff)
		{
			poles_[pos++] = real(pole_seq.get());
			pole_seq.next();
		}
	}
}  // namespace qboot
