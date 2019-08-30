#ifndef QBOOT_POLE_DATA_HPP_
#define QBOOT_POLE_DATA_HPP_

#include <cassert>  // for assert
#include <cstdint>  // for uint32_t, int32_t
#include <utility>  // for move

#include "context_variables.hpp"  // for Context
#include "matrix.hpp"             // for Vector
#include "polynomial.hpp"         // for Polynomial
#include "real.hpp"               // for real, pow, log

namespace qboot
{
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	algebra::Vector<Real> sample_points(uint32_t degree)
	{
		algebra::Vector<Real> v(degree + 1);
		for (uint32_t i = 0; i <= degree; ++i) v[i] = mpfr::pow(Real::pi() * (-1 + 4 * int32_t(i)), 2uL);
		v /= (-64 * int32_t(degree + 1)) * mpfr::log(3 - mpfr::sqrt(Real(8)));
		return v;
	}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class PoleSequence
	{
		bool include_odd;
		uint32_t type, k, spin;
		Real epsilon;

	public:
		using value_type = Real;
		PoleSequence(uint32_t type, uint32_t spin, const Real& epsilon, bool include_odd = true)
		    : include_odd(include_odd), type(type), k(1), spin(spin), epsilon(epsilon)
		{
			assert(1 <= type && type <= 3);
			if (type != 2 && !include_odd) k = 2;
		}
		[[nodiscard]] bool valid() const { return type != 3 || k <= spin; }
		[[nodiscard]] Real get() const
		{
			switch (type)
			{
			case 1: return Real(-int32_t(spin + k - 1));
			case 2: return epsilon - (k - 1);
			default: return 2 * epsilon + (1 + spin - k);  // note: k <= spin -> 1 + spin - k >= 1
			}
		}
		void next() & { k += type == 2 || include_odd ? 1 : 2; }
	};
	template <class L, class R, class Real = typename L::value_type>
	class Merged
	{
		L seql;
		R seqr;
		bool next_l{};
		void update() & { next_l = !seqr.valid() || (seql.valid() && seql.get() >= seqr.get()); }

	public:
		using value_type = Real;
		Merged(L&& l, R&& r) : seql(std::move(l)), seqr(std::move(r)) { update(); }
		void next() &
		{
			if (next_l)
				seql.next();
			else
				seqr.next();
			update();
		}
		[[nodiscard]] bool valid() const { return seql.valid() || seqr.valid(); }
		[[nodiscard]] Real get() const { return next_l ? seql.get() : seqr.get(); }
	};
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	bool include_odd(const Real& d1, const Real& d2, const Real& d3, const Real& d4)
	{
		return d1 != d2 && d3 != d4;
	}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	bool include_odd(const Real& d12, const Real& d34)
	{
		return d12 != 0 && d34 != 0;
	}
	template <class Real = mpfr::real<1000, MPFR_RNDN>>
	class RationalApproxData
	{
		uint32_t spin, lambda;
		const Real &epsilon, &rho;
		Real unitarity_bound{};
		algebra::Vector<Real> poles;

	public:
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& d1, const Real& d2,
		                   const Real& d3, const Real& d4)
		    : RationalApproxData(cutoff, spin, context, d1 - d2, d3 - d4)
		{
		}
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Context<Real>& context, const Real& delta12,
		                   const Real& delta34)
		    : RationalApproxData(cutoff, spin, context, include_odd(delta12, delta34))
		{
		}
		RationalApproxData(uint32_t cutoff, uint32_t spin, const Context<Real>& context, bool include_odd)
		    : spin(spin), lambda(context.lambda()), epsilon(context.epsilon()), rho(context.rho()), poles(cutoff)
		{
			// type 1 or 3 PoleData vanishes when delta12 == 0 or delta34 == 0 and k is odd
			auto get_pols = [this, include_odd](uint32_t type) {
				return PoleSequence<Real>(type, this->spin, this->epsilon, include_odd);
			};
			auto pole_seq = Merged(Merged(get_pols(1), get_pols(2)), get_pols(3));
			uint32_t pos = 0;
			while (pos < cutoff)
			{
				poles[pos++] = pole_seq.get();
				pole_seq.next();
			}
		}
		[[nodiscard]] uint32_t max_degree() const { return poles.size() + lambda; }
		[[nodiscard]] const algebra::Vector<Real>& get_poles() const { return poles; }
		[[nodiscard]] Real get_scale(const Real& delta) const
		{
			Real ans = mpfr::pow(4 * rho, delta);
			for (uint32_t i = 0; i < poles.size(); ++i) ans /= delta - poles[i];
			return ans;
		}
		[[nodiscard]] algebra::Vector<algebra::Polynomial<Real>> get_bilinear_basis(const Real& gap) const
		{
			// orthogonal polynomial of weight function (4 rho) ^ {Delta} / \prod_i (Delta - poles[i])
			algebra::Vector<Real> shifted_poles(poles.size());
			for (uint32_t i = 0; i < poles.size(); ++i) shifted_poles[i] = poles[i] - gap;
			auto weight = fast_partial_fraction(shifted_poles);
			auto deg = max_degree() / 2;
			// inner_prods[i] = \int_{0}^{\infty} dx (4 rho) ^ x x ^ i / \prod_i (x - poles[i])
			algebra::Vector<Real> inner_prods(2 * deg + 1);
			for (uint32_t i = 0; i < poles.size(); ++i)
				inner_prods += mul_scalar(weight[i], simple_pole_integral(2 * deg, 4 * rho, shifted_poles[i]));
			auto mat = anti_band_to_inverse(inner_prods);
			algebra::Vector<algebra::Polynomial<Real>> q(deg + 1);
			for (uint32_t i = 0; i <= deg; ++i)
			{
				algebra::Vector<Real> v(i + 1);
				for (uint32_t j = 0; j <= i; ++j) v[j] = mat.at(i, j);
				q[i] = algebra::Polynomial<Real>(std::move(v));
			}
			return q;
		}
		[[nodiscard]] algebra::Vector<Real> sample_points() const { return qboot::sample_points<Real>(max_degree()); }
	};
}  // namespace qboot

#endif  // QBOOT_CONTEXT_VARIABLES_HPP_
