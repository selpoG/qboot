#include "qboot/bootstrap_equation.hpp"

#include <algorithm>   // for any_of
#include <functional>  // for function
#include <iostream>    // for cout, endl
#include <memory>      // for unique_ptr
#include <string>      // for string
#include <utility>     // for move
#include <vector>      // for vector

#include "qboot/task_queue.hpp"  // for _parallel_evaluate, _event_base

using qboot::algebra::Vector, qboot::algebra::Matrix, qboot::algebra::Polynomial;
using qboot::mp::real, qboot::mp::rational, qboot::mp::integer;
using std::move, std::unique_ptr, std::make_unique, std::cout, std::endl, std::string, std::string_view, std::vector,
    std::optional, std::map, std::function, std::tuple, std::shared_ptr, std::array;

namespace qboot
{
	class ExactPolynomial
	{
		// coeff_[N_] must be positive
		vector<integer> coeff_{};
		// degree
		uint32_t N_;
		// binom(n, m) for 0<=n<=deg, 0<=m<=n/2
		shared_ptr<vector<vector<integer>>> binom_{};
		void set_binom()
		{
			binom_ = std::make_shared<vector<vector<integer>>>();
			binom_->emplace_back();
			binom_->at(0).emplace_back(1);
			binom_->at(0).emplace_back(0);
			if (N_ < 1) return;
			binom_->emplace_back();
			binom_->at(1).emplace_back(1);
			if (N_ < 2) return;
			binom_->emplace_back();
			binom_->at(2).emplace_back(1);
			binom_->at(2).emplace_back(2);
			for (uint32_t n = 3; n <= N_; ++n)
			{
				binom_->emplace_back();
				binom_->at(n).emplace_back(1);
				for (uint32_t m = 1; 2 * m < n; ++m)
					binom_->at(n).emplace_back(binom_->at(n - 1)[m] + binom_->at(n - 1)[m - 1]);
				if (n % 2 == 0) binom_->at(n).emplace_back(binom_->at(n - 1)[n / 2 - 1] << 1);
			}
		}
		[[nodiscard]] const integer& binom(uint32_t n, uint32_t m) const
		{
			if (m > n) return binom_->at(0)[1];
			if (m > n - m) return binom(n, n - m);
			return binom_->at(n)[m];
		}
		ExactPolynomial(vector<integer>&& coeff, const shared_ptr<vector<vector<integer>>>& binom)
		    : coeff_(move(coeff)), N_(0), binom_(binom)
		{
			assert(!coeff_.empty() && coeff_.back() != 0);
			N_ = uint32_t(coeff_.size()) - 1;
			if (coeff_.back() < 0)
				for (auto& c : coeff_) c.negate();
		}

	public:
		explicit ExactPolynomial(vector<integer>&& coeff) : coeff_(move(coeff)), N_(0)
		{
			assert(!coeff_.empty() && coeff_.back() != 0);
			N_ = uint32_t(coeff_.size()) - 1;
			if (coeff_.back() < 0)
				for (auto& c : coeff_) c.negate();
			set_binom();
		}
		explicit ExactPolynomial(const Polynomial& d) : N_(0)
		{
			assert(!d.iszero());
			vector<rational> coeff;  // d as rational
			N_ = uint32_t(d.degree());
			for (uint32_t i = 0; i <= N_; ++i) coeff.push_back(d.at(i).to_rational());
			integer l(1);  // lcm of denominators of coeff
			for (uint32_t i = 0; i <= N_; ++i) l = lcm(l, coeff[i].den());
			vector<integer> coeff2;  // coeff as integer
			for (uint32_t i = 0; i <= N_; ++i) coeff_.push_back(coeff[i].num() * (l / coeff[i].den()));
			integer g(0);  // gcd of coeff2
			for (uint32_t i = 0; i <= N_; ++i) g = gcd(g, coeff_[i]);
			for (uint32_t i = 0; i <= N_; ++i) coeff_[i] /= g;
			assert(coeff_.back() != 0);
			if (coeff_.back() < 0)
				for (auto& c : coeff_) c.negate();
			set_binom();
		}
		[[nodiscard]] uint32_t degree() const noexcept { return N_; }
		[[nodiscard]] rational eval(const rational& x) const
		{
			auto s = rational(coeff_[N_]);
			for (uint32_t i = N_ - 1; i <= N_; --i)
			{
				s *= x;
				s += coeff_[i];
			}
			return s;
		}
		[[nodiscard]] integer eval(const integer& x) const
		{
			auto s = coeff_[N_];
			for (uint32_t i = N_ - 1; i <= N_; --i)
			{
				s *= x;
				s += coeff_[i];
			}
			return s;
		}
		[[nodiscard]] integer roots_upper_bound() const
		{
			integer k(1);
			while (count_roots_pos(k) > 0) k *= 2;
			return k + 1;
		}
		[[nodiscard]] integer roots_lower_bound() const
		{
			integer k(-1);
			while (count_roots_neg(k) > 0) k *= 2;
			return k - 1;
		}
		// num of roots in (lb, ub)
		[[nodiscard]] uint32_t count_roots(const rational& lb, const rational& ub) const
		{
			return normalize(lb, ub).count_roots();
		}
		[[nodiscard]] vector<array<rational, 2>> isolate() const
		{
			auto lb = rational(roots_lower_bound());
			return isolate(rational(roots_lower_bound()));
		}
		[[nodiscard]] vector<array<rational, 2>> isolate(const rational& lb) const
		{
			auto ub = rational(roots_upper_bound());
			vector<array<rational, 2>> ans;
			ans.reserve(N_);
			normalize(lb, ub).isolate_helper(&ans);
			for (auto&& t : ans)
				for (uint32_t i = 0; i < 2; ++i) t.at(i) = lb + (ub - lb) * t.at(i);
			return ans;
		}

	private:
		[[nodiscard]] ExactPolynomial normalize(const rational& lb, const rational& ub) const
		{
			assert(ub > lb);
			// p(lb + (ub - lb) * x)
			// c^n p((a * x + b) / c) = sum(p[i] (ax+b)^i c^(n-i), 0<=i<=n)
			// = sum(p[i] C(i,j) a^j x^j b^(i-j) c^(n-i), 0<=j<=i<=n)
			// = sum(a^j sum(p[i] C(i,j) b^(i-j) c^(n-i), j<=i<=n) x^j, 0<=j<=n)
			// = sum(q[j] x^j, 0<=j<=n)
			// q[n] = a^n p[n]
			vector<integer> q(N_ + 1);
			// b/c = lb, a/c = ub-lb = d
			auto d = ub - lb;
			auto c = lcm(lb.den(), d.den());
			auto a = d.num() * (c / d.den());
			auto b = lb.num() * (c / lb.den());
			auto g = gcd(gcd(a, b), c);
			if (g > 1)
			{
				a /= g;
				b /= g;
				c /= g;
			}
			vector<integer> pow_b, pow_c;
			pow_b.emplace_back(1);
			pow_c.emplace_back(1);
			for (uint32_t i = 1; i <= N_; ++i) pow_b.emplace_back(pow_b.back() * b);
			for (uint32_t i = 1; i <= N_; ++i) pow_c.emplace_back(pow_c.back() * c);
			integer pow_a(1);
			for (uint32_t j = 0; j <= N_; ++j)
			{
				for (uint32_t i = j; i <= N_; ++i) q[j] += coeff_[i] * binom(i, j) * pow_b[i - j] * pow_c[N_ - i];
				q[j] *= pow_a;
				pow_a *= a;
			}
			return ExactPolynomial(move(q), binom_);
		}
		// num of roots in (a, infty)
		[[nodiscard]] uint32_t count_roots_pos(const integer& a) const
		{
			integer q, t;
			vector<integer> pow_a;
			pow_a.emplace_back(1);
			for (uint32_t i = 1; i <= N_; ++i) pow_a.emplace_back(pow_a.back() * a);
			bool positive = true;
			uint32_t cnt = 0;
			for (uint32_t i = N_ - 1; i <= N_; --i)
			{
				q = 0;
				for (uint32_t k = 0; k <= N_ - i; ++k) q += coeff_[i + k] * binom(i + k, k) * pow_a[k];
				if (q != 0)
				{
					auto f = q > 0;
					if (f != positive) ++cnt;
					positive = f;
				}
				if (cnt >= 2) return cnt;
			}
			return cnt;
		}
		// num of roots in (-infty, a)
		[[nodiscard]] uint32_t count_roots_neg(const integer& a) const
		{
			integer q, t;
			vector<integer> pow_a;
			pow_a.emplace_back(1);
			for (uint32_t i = 1; i <= N_; ++i) pow_a.emplace_back(pow_a.back() * a);
			bool positive = true;
			uint32_t cnt = 0;
			for (uint32_t i = N_ - 1; i <= N_; --i)
			{
				q = 0;
				for (uint32_t k = 0; k <= N_ - i; ++k) q += coeff_[i + k] * binom(i + k, k) * pow_a[k];
				if (q != 0)
				{
					bool f = ((N_ - i) % 2 == 1) ^ (q > 0);
					if (f != positive) ++cnt;
					positive = f;
				}
				if (cnt >= 2) return cnt;
			}
			return cnt;
		}
		// num of roots in (0, 1)
		[[nodiscard]] uint32_t count_roots() const
		{
			// num of roots of p(1/(1+x)) in (0,infty)
			// (1+x)^n p(1/(1+x)) = sum(p[i] (1+x)^(n-i), 0<=i<=n)
			// = sum(p[i] C(n-i,j) x^j, 0<=i<=n, 0<=j<=n-i)
			// = sum(sum(p[i] C(n-i,j), 0<=i<=n-j) x^j, 0<=j<=n)
			// r[n] = p[0]
			bool positive = true, zero = true;
			uint32_t cnt = 0;
			integer r;
			for (uint32_t i = N_; i <= N_; --i)
			{
				r = 0;
				for (uint32_t j = 0; j <= N_ - i; ++j) addmul(r, coeff_[j], binom(N_ - j, i));
				if (r == 0) continue;
				if (zero)
				{
					zero = false;
					positive = r > 0;
				}
				else
				{
					auto f = r > 0;
					if (f != positive) ++cnt;
					positive = f;
				}
				if (cnt >= 2) return cnt;
			}
			return cnt;
		}
		// 2^n p(x/2) = sum(p[i] x^i 2^(n-i), 0<=i<=n)
		void left_child()
		{
			for (uint32_t i = 0; i <= N_; ++i) coeff_[i] <<= N_ - i;
		}
		// p(x+1) = sum(p[i] (x+1)^i, 0<=i<=n)
		// = sum(p[i] C(i,j) x^j, 0<=j<=i<=n)
		// = sum(sum(p[i] C(i,j), j<=i<=n) x^j, 0<=j<=n)
		[[nodiscard]] ExactPolynomial shift_one() const
		{
			vector<integer> q(N_ + 1);
			for (uint32_t j = 0; j <= N_; ++j)
				for (uint32_t i = j; i <= N_; ++i) addmul(q[j], coeff_[i], binom(i, j));
			return ExactPolynomial(move(q), binom_);
		}
		void isolate_helper(vector<array<rational, 2>>* ans) &&
		{
			rational m(1, 2u);
			auto p_m = eval(m);
			auto cnt = count_roots();
			if (p_m == 0) ans->push_back({m, m});
			if (cnt == 1 && p_m != 0) ans->push_back({rational(0), rational(1)});
			if (cnt >= 2)
			{
				left_child();
				auto right = shift_one();
				auto it = ans->end();
				move(*this).isolate_helper(ans);
				for (; it != ans->end(); ++it)
					for (uint32_t i = 0; i < 2; ++i) it->at(i) /= 2;
				vector<integer>().swap(coeff_);
				move(right).isolate_helper(ans);
				for (; it != ans->end(); ++it)
					for (uint32_t i = 0; i < 2; ++i) it->at(i) = (1 + it->at(i)) / 2;
			}
		}

	public:
		friend std::ostream& operator<<(std::ostream& out, const ExactPolynomial& p)
		{
			auto f = false;
			{
				auto& val = p.coeff_[0];
				if (!val.iszero())
				{
					out << val;
					f = true;
				}
			}
			for (uint32_t i = 1; i <= p.N_; ++i)
			{
				auto& val = p.coeff_[i];
				auto sgn = mp::sgn(val);
				if (sgn == 0) continue;
				if (f)
				{
					if (val == 1)
						out << " + ";
					else if (val == -1)
						out << " - ";
					else if (sgn > 0)
						out << " + " << val << " * ";
					else
						out << " - " << -val << " * ";
				}
				else if (val == -1)
					out << "-";
				else if (val != 1)
					out << val << " * ";
				out << "x";
				if (i > 1) out << "^" << i;
				f = true;
			}
			return out;
		}
	};

	SDPMode::~SDPMode() = default;
	_find_contradiction::~_find_contradiction() = default;
	_extremal_OPE::~_extremal_OPE() = default;
	PolynomialProgram _find_contradiction::_prepare(uint32_t N, const BootstrapEquation& boot) const
	{
		PolynomialProgram prg(N);
		prg.objective_constant() = real(0);
		prg.objectives(Vector(prg.num_of_variables(), real(0)));
		prg.add_equation(boot.make_disc_mat_v(norm_), real(1));
		return prg;
	}
	PolynomialProgram _extremal_OPE::_prepare(uint32_t N, const BootstrapEquation& boot) const
	{
		PolynomialProgram prg(N);
		prg.objective_constant() = real(0);
		prg.objectives(boot.make_disc_mat_v(norm_));
		prg.add_equation(boot.make_disc_mat_v(target_), real(maximize_ ? 1 : -1));
		return prg;
	}

	Vector<real> read_raw_functional(const fs::path& y_txt)
	{
		std::ifstream is(y_txt);
		uint32_t sz, M;
		is >> sz >> M;
		assert(M == 1);
		Vector<real> v(sz);
		for (uint32_t i = 0; i < sz; ++i) is >> v.at(i);
		return v;
	}
	using FixedBlock = ConformalBlock<PrimaryOperator>;
	using GeneralBlock = GeneralConformalBlock;
	Matrix<real> BootstrapEquation::apply_functional(const Vector<real>& func, string_view disc_sector) const
	{
		auto id = get_id(disc_sector);
		const auto& sec = sector(id);
		assert(sec.type() == SectorType::Discrete);
		auto mats = make_disc_mat(id);
		auto sz = sec.size();
		Matrix<real> mat{sz, sz};
		for (uint32_t n = 0; n < N_; ++n) mat += mul_scalar(func[n], move(mats[n]));
		return mat;
	}
	[[nodiscard]] Vector<Matrix<real>> BootstrapEquation::make_disc_mat(uint32_t id) const
	{
		assert(sector(id).type() == SectorType::Discrete);
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
			include_odd |= std::any_of(eqs_[i][id].begin(), eqs_[i][id].end(), [](const auto& term) {
				return std::get<GeneralBlock>(term.block()).include_odd();
			});
		return make_unique<ConformalScale>(op, cont_, include_odd);
	}
	[[nodiscard]] Vector<Vector<Matrix<real>>> BootstrapEquation::make_cont_mat(
	    uint32_t id, const GeneralPrimaryOperator& op, const unique_ptr<ConformalScale>& ag) const
	{
		assert(sector(id).type() == SectorType::Continuous);
		auto sz = sector(id).size();
		auto sp = ag->sample_points();
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
					auto delta = ag->get_delta(sp[k]);
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

	map<string, vector<PrimaryOperator>, std::less<>> BootstrapEquation::get_spectrum(
	    const Vector<real>& recovered_func, uint32_t parallel, const unique_ptr<_event_base>& event) const
	{
		assert(N_ > 0);
		map<string, vector<PrimaryOperator>, std::less<>> spectrum_dict;
		_memoized<Matrix<real>(uint32_t)> inv(
		    [](uint32_t deg) { return algebra::interpolation_matrix(ConformalScale::sample_points(deg)); }, parallel);
		vector<function<tuple<string, vector<PrimaryOperator>>()>> spectrums;
		for (const auto& [sec, id] : sector_id_)
		{
			if (sector(id).type_ != SectorType::Continuous) continue;
			auto sz = sector(id).size();
			for (const auto& op : sector(id).ops_)
				spectrums.emplace_back([this, id = id, sec = sec, &op, sz, &recovered_func, &inv, &event]() {
					vector<PrimaryOperator> spectrum;
					auto tag = op.str();
					tag += " in ";
					tag += sec;
					cout << tag << endl;
					_scoped_event scope(tag, event);
					auto ag = common_scale(id, op);
					auto ps = ag->sample_points();
					auto mat = make_cont_mat(id, op, ag);
					auto num_pts = ps.size();
					assert(N_ == mat.size());
					Vector<Matrix<real>> mats(num_pts);
					for (uint32_t k = 0; k < num_pts; ++k)
					{
						mats[k] = {sz, sz};
						for (uint32_t n = 0; n < N_; ++n) mats[k] += mul_scalar(recovered_func[n], move(mat[n][k]));
						mats[k] /= ag->eval(ps[k]);
					}
					auto mat_pol = polynomial_interpolate(mats, inv(ag->max_degree()));
					auto det = determinant(mat_pol);
					auto d = det.derivative();
					auto d2 = d.derivative();
					// TODO(selpo): roots of pol may be multiple (square-free factorization by Yun's algorithm)
					ExactPolynomial pol(d), pol2(d2);
					for (auto&& t : pol.isolate(rational(0)))
					{
						if (t[0] == t[1])
						{
							spectrum.push_back(op.fix_delta(ag->get_delta(real(t[0]))));
							continue;
						}
						rational lb(t[0]), ub(t[1]);
						auto sgn = pol.eval(lb) > 0;
						auto flag = true;
						while (flag)
						{
							auto x = (lb + ub) / 2;
							auto ev = pol.eval(x);
							if (sgn == (ev > 0))
								lb = x;
							else if (ev != 0)
								ub = x;
							else if (pol2.eval(x) > 0)
							{
								spectrum.push_back(op.fix_delta(ag->get_delta(real(x))));
								flag = false;
								break;
							}
							if (pol2.count_roots(lb, ub) == 0) break;
						}
						if (!flag) continue;
						if (pol2.eval((lb + ub) / 2) < 0) continue;
						real x((lb + ub) / 2);
						for (uint32_t i = 0; i < 15; ++i) x = x - d.eval(x) / d2.eval(x);
						assert(t[0] <= x && x <= t[1]);
						spectrum.push_back(op.fix_delta(ag->get_delta(x)));
						cout << x << ": " << det.eval(x) / d2.eval(x) << endl;
					}
					cout << "at zero: " << det.eval(0) / d.eval(0) << endl;
					return tuple{sec, spectrum};
				});
		}
		for (auto&& [sec, ops] : _parallel_evaluate(spectrums, parallel))
			for (auto& op : ops) spectrum_dict[sec].push_back(op);
		return spectrum_dict;
	}

	void BootstrapEquation::finish() &
	{
		for (const auto& eq : eqs_) N_ += eq.dimension();
	}

	[[nodiscard]] Vector<real> BootstrapEquation::recover(const unique_ptr<SDPMode>& mode,
	                                                      const Vector<real>& raw_func) const
	{
		assert(N_ > 0);
		return mode->_prepare(N_, *this).recover(raw_func);
	}
	[[nodiscard]] PolynomialProgram BootstrapEquation::convert(const unique_ptr<SDPMode>& mode, uint32_t parallel,
	                                                           const unique_ptr<_event_base>& event) const
	{
		assert(N_ > 0);
		auto prg = mode->_prepare(N_, *this);
		auto filter = mode->_get_filter();
		vector<function<optional<PolynomialInequality>()>> ineqs;
		for (const auto& [sec, id] : sector_id_)
		{
			if (!filter(sec)) continue;
			auto sz = sector(id).size();
			if (sector(id).type() == SectorType::Discrete)
				if (sector(id).is_matrix())
					ineqs.emplace_back([this, id = id, sz, sec = sec, &event] {
						_scoped_event scope(sec, event);
						return optional{PolynomialInequality(N_, sz, make_disc_mat(id), Matrix<real>(sz, sz))};
					});
				else
					ineqs.emplace_back([this, id = id, sec = sec, &event] {
						_scoped_event scope(sec, event);
						return optional{PolynomialInequality(N_, make_disc_mat_v(id), real(0))};
					});
			else
				for (const auto& op : sector(id).ops_)
					ineqs.emplace_back([this, id = id, sz, &op, sec = sec, &event] {
						auto tag = op.str();
						tag += " in ";
						tag += sec;
						_scoped_event scope(tag, event);
						auto ag = common_scale(id, op);
						auto mat = make_cont_mat(id, op, ag);
						auto deg = ag->max_degree();
						return optional{
						    PolynomialInequality(N_, sz, move(ag), move(mat), Vector<Matrix<real>>(deg + 1, {sz, sz}))};
					});
		}
		for (auto&& x : _parallel_evaluate(ineqs, parallel)) prg.add_inequality(move(x));
		return prg;
	}
}  // namespace qboot
