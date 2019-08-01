# qboot

Converter from Conformal Bootstrap Equations to Semidefinite Programming,
which admits Mobius transformation of Delta.

Some codes are taken from [cboot](https://github.com/tohtsky/cboot).

## Requirements

- [cmake](https://cmake.org/) (`3.13.0+`)

### Unix

- [gmp](https://gmplib.org/) (You must configure `gmp` with `--enable-cxx`)

- [mpfr](http://mpfr.org/)

- [gcc](http://gcc.gnu.org/) (`7.4.0+`) or [clang](http://clang.llvm.org/) (`8.0.0+`)

### Windows

- [mpir](https://github.com/BrianGladman/mpir) (You must build `lib_mpir_gc` (or `lib_mpir_{some_architecture}`) **and** `lib_mpir_cxx`)

- [mpfr](https://github.com/BrianGladman/mpfr)

- [Visual Studio](https://visualstudio.microsoft.com/) (`2017+`)

After building these libraries with Visual Studio,
you will have `gmp.h`, `gmpxx.h`, `mpir.lib` and `mpirxx.lib` in, for example, `C:\somewhere\mpir\lib\x64\Debug\`
and `mpfr.lib`, `mpfr.h` in, for example, `C:\somewhere\mpfr\lib\x64\Debug\`.

If you are working with `WSL`, please follow the `Unix` section above.

## Build

### Unix

1. `mkdir qboot/build && cd qboot/build`

2. `cmake ..`

3. `make`

If you installed requirements in custom location, You may need to tell `cmake` some paths.
In Unix system, add options `-DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp`.

(For developers: If you want to use `clang-tidy`, add `-DCMAKE_CXX_CLANG_TIDY="clang-tidy"` to `cmake` options.)

### Windows

In Windows, use cmake-gui to build.
Options to cmake can be passed by `Add Entry` button.
You have to add 2 `PATH`s,
`GMP_ROOT` (must contain `gmp.h`, `gmpxx.h` and `mpir.lib`. Typical value is `C:\somewhere\mpir\lib\x64\Debug`) and
`MPFR_ROOT` (must contain `mpfr.h` and `mpfr.lib`. Typical value is `C:\somewhere\mpfr\lib\x64\Debug`).
