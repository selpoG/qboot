#ifndef QBOOT_MY_FILESYSTEM_HPP_
#define QBOOT_MY_FILESYSTEM_HPP_

#if __has_include(<filesystem>)
#include <filesystem>  // for path, create_directory
namespace qboot
{
	namespace fs = std::filesystem;
}  // namespace qboot
#else
#include <experimental/filesystem>
namespace qboot
{
	namespace fs = std::experimental::filesystem;
}  // namespace qboot
#endif

#endif  // QBOOT_MY_FILESYSTEM_HPP_
