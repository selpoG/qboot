#include "io.hpp"
#include "multi_array.hpp"
#include "my_hash.hpp"

[[maybe_unused]] static void test()
{
	util::MyHash()(2);
	util::MyHash()(-1.0);
	util::MyHash()(
	    std::tuple(4, std::pair(2, std::optional(-1)), std::variant<size_t, double>(-0.0), std::array{3, 1, 4}));
}
