# qboot

Converter from Conformal Bootstrap Equations to Semidefinite Programming.
Sum of operators may run over not only infinite range `[lb, infty)` but also finite one `[lb, ub)`.

For more information, see
[An Automated Generation of Bootstrap Equations for Numerical Study of Critical Phenomena](https://arxiv.org/abs/2006.04173).

Some codes are taken from [cboot](https://github.com/tohtsky/cboot.git).

- [qboot](#qboot)
	- [Requirements](#requirements)
		- [Unix (or WSL)](#unix-or-wsl)
		- [Windows (MSVC)](#windows-msvc)
	- [Install](#install)
		- [Unix](#unix)
		- [Windows](#windows)
	- [Use installed `qboot`](#use-installed-qboot)
	- [Docker](#docker)

## Requirements

- [cmake](https://cmake.org/) (`3.12.4+`)

### Unix (or WSL)

- [gmp](https://gmplib.org/) (configured with `--enable-cxx`)

- [mpfr](http://mpfr.org/)

- [gcc](http://gcc.gnu.org/) (`7.4.0+`) or [clang](http://clang.llvm.org/) (`8.0.0+`)

### Windows (MSVC)

- [mpir](https://github.com/BrianGladman/mpir.git) (You must build `lib_mpir_gc` **and** `lib_mpir_cxx`)

- [mpfr](https://github.com/BrianGladman/mpfr.git)

- [Visual Studio](https://visualstudio.microsoft.com/) (`2017+`)

Please build `mpir` and `mpfr` in both `Debug` and `Release` mode.

If you cloned `mpir` and `mpfr` in `C:\somewhere\mpir` and `C:\somewhere\mpfr`,
you will have

- `C:\somewhere\mpir\lib\x64\Debug\gmpxx.h`
- `C:\somewhere\mpir\lib\x64\Release\gmpxx.h`
- `C:\somewhere\mpfr\lib\x64\Debug\mpfr.h`
- `C:\somewhere\mpir\lib\x64\Release\gmpxx.h`

## Install

### Unix

```sh
mkdir qboot/build && cd qboot/build
cmake .. -DMPFR_ROOT=/path/to/mpfr -DGMP_ROOT=/path/to/gmp -DCMAKE_INSTALL_PREFIX=/somewhere/you/like
make
make install
```

To build in `Debug` or `Release` mode, call `cmake` with `-DCMAKE_BUILD_TYPE=Debug` or `-DCMAKE_BUILD_TYPE=Release`.

You have to ensure the exsistence of:

- `MPFR_ROOT/include/mpfr.h`
- `MPFR_ROOT/lib/libmpfr.a`
- `GMP_ROOT/include/gmpxx.h`
- `GMP_ROOT/lib/libgmpxx.a`

### Windows

In Windows, use cmake-gui to build.

1. Create `qboot\build` folder.
2. Set "Where is the source code:" to the path of `qboot` folder,
"Where to build the binaries:" to the path of `qboot\build` folder.
3. Push `Add Entry` to set `GMP_ROOT`, `MPFR_ROOT` and `CMAKE_INSTALL_PREFIX`.
4. Push `Configure` button and `Generate` button.
5. Open `qboot.sln` in `qboot\build` folder and build `INSTALL` project.

You have to ensure the exsistence of:

- `MPFR_ROOT\Debug\mpfr.h`, `MPFR_ROOT\Release\mpfr.h`
- `MPFR_ROOT\Debug\mpfr.lib`, `MPFR_ROOT\Release\mpfr.lib`
- `GMP_ROOT\Debug\gmpxx.h`, `GMP_ROOT\Release\gmpxx.h`
- `GMP_ROOT\Debug\mpirxx.lib`, `GMP_ROOT\Release\mpirxx.lib`

The typical value for `GMP_ROOT` is `C:\somewhere\mpir\lib\x64`
and for `MPFR_ROOT` `C:\somewhere\mpfr\lib\x64`.
You can set `CMAKE_INSTALL_PREFIX` to anywhere you like.

## Use installed `qboot`

[sample](/sample) folder is a sample project which use installed `qboot`.
You can copy this folder to anywhere you like.

Once you have installed `qboot` with `-DCMAKE_INSTALL_PREFIX=/some/where`, you can build this sample by:

```sh
cd sample
mkdir build
cd build
cmake .. -DQBoot_ROOT=/some/where -DCMAKE_BUILD_TYPE=Debug
make
```

And you can execute `sample/build/bin/sample`.

Of course, you can build also in Windows.

1. Create `sample\build` folder.
2. Set "Where is the source code:" to the path of `sample` folder,
"Where to build the binaries:" to the path of `sample\build` folder.
3. Push `Add Entry` to set `QBoot_DIR` to `/some/where`.
4. Push `Configure` button and `Generate` button.

You can build `sample.sln` in `sample\build` folder.

## Docker

Docker for qboot can be obtained from [selpo/qboot](https://hub.docker.com/r/selpo/qboot),
which was built using [Dockerfile](Dockerfile).
We also have a simple sample using this container in [sample/Dockerfile](sample/Dockerfile).

```sh
cd qboot/sample
docker build -t qboot.sample .
docker run -it --rm qboot.sample
./build/bin/sample-debug
```
