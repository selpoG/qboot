# Implementation Details

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
