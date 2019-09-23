# qboot

Converter from Conformal Bootstrap Equations to Semidefinite Programming,
which admits Mobius transformation of Delta.

Some codes are taken from [cboot](https://github.com/tohtsky/cboot.git).

## Requirements

- [cmake](https://cmake.org/) (`3.13.0+`)

### Unix (or WSL)

- [gmp](https://gmplib.org/) (You must configure `gmp` with `--enable-cxx`)

- [mpfr](http://mpfr.org/)

- [gcc](http://gcc.gnu.org/) (`7.4.0+`) or [clang](http://clang.llvm.org/) (`8.0.0+`)

### Windows (MSVC)

- [mpir](https://github.com/BrianGladman/mpir.git) (You must build `lib_mpir_gc` (or `lib_mpir_{some_architecture}`) **and** `lib_mpir_cxx`)

- [mpfr](https://github.com/BrianGladman/mpfr.git)

- [Visual Studio](https://visualstudio.microsoft.com/) (`2017+`)

Please build `mpir` and `mpfr` in both `Debug` and `Release` mode.

After building these libraries with Visual Studio,
you will have `gmp.h`, `gmpxx.h`, `mpir.lib` and `mpirxx.lib`
in, for example, `C:\somewhere\mpir\lib\x64\Debug\`, `C:\somewhere\mpir\lib\x64\Release\`
and `mpfr.lib`, `mpfr.h` in, for example, `C:\somewhere\mpfr\lib\x64\Debug\`, `C:\somewhere\mpfr\lib\x64\Release\`.

## Build

### Unix

1. `mkdir qboot/build && cd qboot/build`

2. `cmake ..`

3. `make`

To build in `Debug` or `Release` mode, call `cmake` with `-DCMAKE_BUILD_TYPE=Debug` or `-DCMAKE_BUILD_TYPE=Release`.

If you installed requirements in custom location, You may need to tell `cmake` some paths.
In Unix system, add options `-DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp`.

(For developers: If you want to use `clang-tidy`, add `-DCMAKE_CXX_CLANG_TIDY="clang-tidy"` to `cmake` options.)

### Windows

In Windows, use cmake-gui to build.
Options to cmake can be passed by `Add Entry` button.
You have to add 2 `PATH`s,
`GMP_ROOT` (must contain `Debug` and `Release` folders which contain `gmp.h`, `gmpxx.h` and `mpir.lib`. Typical value is `C:\somewhere\mpir\lib\x64`) and
`MPFR_ROOT` (must contain `Debug` and `Release` folders which contain `mpfr.h` and `mpfr.lib`. Typical value is `C:\somewhere\mpfr\lib\x64`).

## Data Structures

A Polynomial Matrix Program (PMP) is a problem
to maximize `b[N] + \sum_{n = 0}^{N - 1} b[n] y[n]`
over free real variables `y[0], ..., y[N - 1]`
such that for all `0 <= j < J` and `x >= 0`, `\sum_{n = 0}^{N - 1} y[n] M_j[n](x) \succeq M_j[N](x)`.

`b[0], ..., b[N]` are real constants
and `M_j[0], ..., M_j[N]` are symmetric matrices whose elements `P_j[n][r, c]` are real polynomials of `x`.

We generalize a PMP to a Function Matrix Program (FMP),
which allows each elements in `M_j[0], ..., M_j[N]` to be in the form `\chi_j(x) Q_j[n][r, c](x)`
in which `Q_j[n][r, c](x)` is a real polynomial of `x`
and `\chi_j(x)` is an arbitrary real function of `x` which is positive in `x >= 0`.
For example, in the conformal bootstrap, the typical from of `\chi_j(x)` is `exp(-A x) / ((x + a) ... (x + z))`
(`A, a, ..., z >= 0`).

This generalization seems to be redundant, because we can divide inequality by `\chi_j(x)` to get a PMP.
But, as discussed in section 3.3 in [arXiv:1502.02033](https://arxiv.org/abs/1502.02033),
preserving these factors contributes to numerical stability of [SDPB](https://github.com/davidsd/sdpb.git).
