#include <array>
#include <cassert>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <variant>
#include <vector>

#include "qboot/mp/real.hpp"
#include "qboot/my_filesystem.hpp"
#include "qboot/task_queue.hpp"

namespace mp = qboot::mp;
namespace fs = qboot::fs;
using fs::path;
using mp::real, mp::rational, mp::integer;
using std::ostringstream, std::ifstream, std::endl, std::cout;
using std::string, std::vector, std::unique_ptr, std::move, std::variant, std::array, std::optional;
template <class T>
using my_result = variant<T, string>;
using my_error = my_result<std::monostate>;

template <class T>
bool is_ok(const my_result<T>& v)
{
	return v.index() == 0;
}

template <class T>
string& as_err(my_result<T>& v)
{
	return std::get<1>(v);
}
template <class T>
const string& as_err(const my_result<T>& v)
{
	return std::get<1>(v);
}
template <class T>
string&& as_err(my_result<T>&& v)
{
	return std::get<1>(move(v));
}

template <class T>
T& as_val(my_result<T>& v)
{
	return std::get<0>(v);
}
template <class T>
const T& as_val(const my_result<T>& v)
{
	return std::get<0>(v);
}
template <class T>
T&& as_val(my_result<T>&& v)
{
	return std::get<0>(move(v));
}

namespace qboot::mp
{
	inline real log2(const real& x) { return log(x) / const_log2(); }
}  // namespace qboot::mp

class ElementBase
{
	uint32_t r_, c_;

public:
	constexpr ElementBase(uint32_t r, uint32_t c) noexcept : r_(r), c_(c) {}
	[[nodiscard]] constexpr uint32_t row() const noexcept { return r_; }
	[[nodiscard]] constexpr uint32_t column() const noexcept { return c_; }
	[[nodiscard]] string str() const noexcept { return "r=" + std::to_string(r_) + ", c=" + std::to_string(c_); }
};

class EOFElement : public ElementBase
{
public:
	constexpr EOFElement(uint32_t r, uint32_t c) noexcept : ElementBase(r, c) {}
	[[nodiscard]] string str() const noexcept { return "EOF(" + ElementBase::str() + ")"; }
};

class SpaceElement : public ElementBase
{
public:
	static bool contains(char c) noexcept { return c == ' '; }
	constexpr SpaceElement(uint32_t r, uint32_t c) noexcept : ElementBase(r, c) {}
	[[nodiscard]] string str() const noexcept { return "WhiteSpace(" + ElementBase::str() + ")"; }
};

class NewLineElement : public ElementBase
{
public:
	static bool contains(char c) noexcept { return c == '\n'; }
	constexpr NewLineElement(uint32_t r, uint32_t c) noexcept : ElementBase(r, c) {}
	[[nodiscard]] string str() const noexcept { return "NewLine(" + ElementBase::str() + ")"; }
};

class NumberElement : public ElementBase
{
	real num_;

public:
	static bool contains(char c) noexcept
	{
		return ('0' <= c && c <= '9') || c == '+' || c == 'e' || c == '-' || c == '.';
	}
	NumberElement(uint32_t r, uint32_t c, const real& num) noexcept : ElementBase(r, c), num_(num) {}
	NumberElement(uint32_t r, uint32_t c, real&& num) noexcept : ElementBase(r, c), num_(move(num)) {}
	[[nodiscard]] const real& number() const noexcept { return num_; }
	[[nodiscard]] string str() const noexcept
	{
		return "Number(" + real(num_).str(10) + ", " + ElementBase::str() + ")";
	}
};

using Element = variant<EOFElement, SpaceElement, NewLineElement, NumberElement>;

static string to_str(fs::file_type t);
static string to_str(const Element& x) noexcept;
static real error(const real& p, const real& q) noexcept;
static string error_msg(const Element& x, const Element& y) noexcept;
static bool equal(const Element& x, const Element& y, const real& error_bound) noexcept;
static my_error check_file(const path& d1, const path& d2, const real& error_bound) noexcept;
static my_result<fs::file_type> check_type(const path& x, const path& y);
static my_error check_dir(const path& d1, const path& d2, const real& error_bound, qboot::TaskQueue* q,
                          vector<std::future<my_error>>* futs) noexcept;
static my_error check_dir(const path& d1, const path& d2, const real& error_bound, uint32_t parallel) noexcept;

static string to_str(fs::file_type t)
{
	using fs::file_type;
	switch (t)
	{
	case file_type::none: return "none";
	case file_type::not_found: return "not_found";
	case file_type::regular: return "regular";
	case file_type::directory: return "directory";
	case file_type::symlink: return "symlink";
	case file_type::block: return "block";
	case file_type::character: return "character";
	case file_type::fifo: return "fifo";
	case file_type::socket: return "socket";
	case file_type::unknown: return "unknown";
	default: return "implementation-defined (" + std::to_string(int(t)) + ")";
	}
}

string to_str(const Element& x) noexcept
{
	return std::visit([](auto t) { return t.str(); }, x);
}

real error(const real& p, const real& q) noexcept
{
	auto d = mp::abs(p - q);
	return mp::cmpabs(real(1), p) >= 0 ? d : d / mp::abs(p);
}

string error_msg(const Element& x, const Element& y) noexcept
{
	if (x.index() != y.index()) return "different type of Elements";
	ostringstream os;
	auto p = std::get<NumberElement>(x).number();
	auto q = std::get<NumberElement>(y).number();
	auto err = error(p, q);
	os << "error = " << err << ", log2(err) = " << mp::log2(err);
	return os.str();
}

bool equal(const Element& x, const Element& y, const real& error_bound) noexcept
{
	if (x.index() != y.index()) return false;
	if (std::holds_alternative<NumberElement>(x))
	{
		auto p = std::get<NumberElement>(x).number();
		auto q = std::get<NumberElement>(y).number();
		return error(p, q) <= error_bound;
	}
	return true;
}

class buf_ifstream
{
	static constexpr uint32_t buf_size = 4096;
	static_assert(buf_size >= 1);
	ifstream fin_;
	array<char, buf_size> buf_{};
	uint32_t pos_ = 0, eof_ = 0, r_ = 0, c_ = 0;

public:
	explicit buf_ifstream(const path& path) noexcept : fin_(path, std::ios::binary)
	{
		eof_ = uint32_t(fin_.read(buf_.data(), buf_size).gcount());
	}
	uint32_t row() const noexcept { return r_; }
	uint32_t column() const noexcept { return c_; }
	char peek() const noexcept { return buf_.at(pos_); }
	char get() noexcept
	{
		auto ret = peek();
		next();
		return ret;
	}
	void next() noexcept
	{
		assert(!eof());
		if (peek() == '\n')
		{
			++r_;
			c_ = 0;
		}
		else
			++c_;
		if (++pos_ >= buf_size)
		{
			eof_ = uint32_t(fin_.read(buf_.data(), buf_size).gcount());
			pos_ = 0;
		}
	}
	bool eof() const noexcept { return pos_ == eof_; }
};

class Content
{
	vector<Element> elems_;

	explicit Content(vector<Element>&& elems) noexcept : elems_(move(elems)) {}

public:
	static my_result<Content> read(const path& path) noexcept
	{
		buf_ifstream fin(path);
		vector<Element> elems;
		while (!fin.eof())
		{
			auto ch = fin.get();
			uint32_t r = fin.row(), c = fin.column();
			if (SpaceElement::contains(ch))
				elems.emplace_back(SpaceElement(r, c));
			else if (NewLineElement::contains(ch))
				elems.emplace_back(NewLineElement(r, c));
			else if (NumberElement::contains(ch))
			{
				string str;
				do
					str.push_back(ch);
				while (!fin.eof() && NumberElement::contains(ch = fin.get()));
				elems.emplace_back(NumberElement(r, c, real(str)));
			}
			else
			{
				ostringstream os;
				os << "unrecognized character '" << ch << "' (" << int(ch) << ") at row = " << r << ", column = " << c
				   << " in " << path;
				return os.str();
			}
		}
		elems.emplace_back(EOFElement(fin.row(), 0));
		return Content(move(elems));
	}
	static optional<array<Element, 2>> compare(const Content& x, const Content& y, const real& error_bound) noexcept
	{
		for (uint32_t i = 0; i < x.elems_.size() && i < y.elems_.size(); ++i)
			if (!equal(x.elems_[i], y.elems_[i], error_bound)) return array{x.elems_[i], y.elems_[i]};
		return {};
	}
};

my_error check_file(const path& d1, const path& d2, const real& error_bound) noexcept
{
	for (const auto& d : {d1, d2})
		if (!fs::is_regular_file(d)) return d.string() + " is not a file";
	auto c1 = Content::read(d1);
	if (!is_ok(c1)) return as_err(move(c1));
	auto c2 = Content::read(d2);
	if (!is_ok(c2)) return as_err(move(c2));
	auto diff = Content::compare(as_val(c1), as_val(c2), error_bound);
	if (!diff.has_value()) return std::monostate{};
	auto [x, y] = diff.value();
	ostringstream os;
	os << "x (in " << d1 << ") != y (in " << d2 << ")" << endl;
	os << "x = " << to_str(x) << endl;
	os << "y = " << to_str(y) << endl;
	os << error_msg(x, y);
	return os.str();
}

my_result<fs::file_type> check_type(const path& x, const path& y)
{
	auto tx = fs::status(x).type();
	auto ty = fs::status(y).type();
	if (tx == ty) return tx;
	ostringstream os;
	os << x << " and " << y << " has different file type" << endl;
	os << x << " has a type:" << to_str(tx) << endl;
	os << y << " has a type:" << to_str(ty);
	return os.str();
}

my_error check_dir(const path& d1, const path& d2, const real& error_bound, qboot::TaskQueue* q,
                   vector<std::future<my_error>>* futs) noexcept
{
	for (const auto& d : {d1, d2})
		if (!fs::is_directory(d)) return d.string() + " is not a directory";
	for (const auto& x : fs::directory_iterator(d1))
	{
		auto name = x.path().filename();
		auto y = d2 / name;
		auto _type = check_type(x.path(), y);
		if (!is_ok(_type)) return as_err(_type);
		auto type = as_val(_type);
		if (type == fs::file_type::regular)
			futs->push_back(q->push([x = x, y = y, &error_bound] { return check_file(x.path(), y, error_bound); }));
		else if (type == fs::file_type::directory)
			if (auto err = check_dir(x.path(), y, error_bound, q, futs); !is_ok(err)) return err;
	}
	for (const auto& y : fs::directory_iterator(d2))
		if (auto type = check_type(y.path(), d1 / y.path().filename()); !is_ok(type)) return as_err(type);
	return std::monostate{};
}

my_error check_dir(const path& d1, const path& d2, const real& error_bound,
                   uint32_t parallel = std::thread::hardware_concurrency()) noexcept
{
	qboot::TaskQueue q(parallel);
	vector<std::future<my_error>> futs;
	if (auto err = check_dir(d1, d2, error_bound, &q, &futs); !is_ok(err)) return err;
	for (auto&& x : futs)
		if (auto err = x.get(); !is_ok(err)) return err;
	return std::monostate{};
}

int main(int argc, char* argv[])
{
	mp::global_prec = 1200;
	mp::global_rnd = MPFR_RNDN;
	auto error_bound = real(mp::parse("100e-240").value());
	if (argc != 3 && argc != 4)
	{
		cout << "usage: ./comp dir1 dir2 [error]" << endl;
		return 1;
	}
	path d1, d2;
	{
		unique_ptr<char*[]> args(argv);
		d1 = args[1];
		d2 = args[2];
		if (argc == 4) error_bound = mp::parse(args[3]).value();
		args.release();
	}
	auto err = check_dir(d1, d2, error_bound);
	if (!is_ok(err))
	{
		cout << as_err(err) << endl;
		return 1;
	}
	return 0;
}
