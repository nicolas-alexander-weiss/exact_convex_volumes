# Exact Volumes of Semi-Algebraic Convex Bodies

[![DOI](https://zenodo.org/badge/1033727412.svg)](https://doi.org/10.5281/zenodo.18472470) [![arXiv](https://img.shields.io/badge/arXiv-2602.04707-red.svg)](https://arxiv.org/abs/2602.04707)

## Description

We provide a SageMath package to compute the exact volume of semi-algebraic convex bodies. The input consists of $k$ polynomials
$`f_1,f_2,\ldots, f_k  \in  \mathbb{Q}[x_1,\ldots,x_n],`$
that each define a _concave_ function on $`\mathbb{R}^n`$.

The output is a (real/complex) ball of radius $`2^{-\text{prec}}`$ which contains the volume of the convex semi-algebraic set
$`
C = \left\{x \in \mathbb{R}^n \mid f_1(x) > 0\right\} \cap \ldots \cap \left\{x \in \mathbb{R}^n \mid f_k(x) > 0\right\}.
`$

The computation of exact volumes of semi-algebraic convex bodies defined by finitely many concave polynomials is  motivated by geometric statistics, where intersections of convex bodies arise as maximum likelihood estimator (MLE) sets.

This package implements the algorithms described in

  [1] L. Ramesh, N. Weiss (2026): "Exact Volumes of Semi-Algebraic Convex Bodies", ISSAC'26: [doi.org/10.1145/3815436.3815478](https://doi.org/10.1145/3815436.3815478) [arXiv:2602.04707](https://arxiv.org/abs/2602.04707).

which rests on work by [Lairez, Mezzarobba, and Safey El Din](http://dx.doi.org/10.1145/3326229.3326262).

The experimental results of [1] can be reproduced using the Jupyter notebook:

- "examples/Examples from Ramesh, Weiss (2026).ipynb"

This is a project of [Lakshmi Ramesh](https://lax202.github.io) and [Nicolas Weiss](https://nicolas-alexander-weiss.github.io).


## Requirements

This package was developed using SageMath (10.6) and using Python (3.12.5). It requires installation of:

SageMath Packages
- [ore_algebra](https://github.com/mkauers/ore_algebra/)

External Software:
- [msolve](https://msolve.lip6.fr/downloads/msolve-tutorial.pdf)


(Optional) Julia Package:
- [HypersurfaceRegions.jl](https://github.com/JuliaAlgebra/HypersurfaceRegions.jl)
- [MultivariateCreativeTelescoping.jl](https://github.com/HBrochet/MultivariateCreativeTelescoping.jl)


## Installation

Begin by installing the requirements above. 

Instructions for installing the ore_algebra package can be found on their [Github page](https://github.com/mkauers/ore_algebra/). Note that Macaulay2 and msolve have to be accessible from the command line via "M2" and "msolve" respectively.

To install our package, download this repository and then install it via

``sage -pip install <path-to-"exact_convex_volumes"-repo>``

If you would like to avoid installing the package you can simply add the folder "src" to the sys.path.

## Using MultivariateCreativeTelescoping.jl

By default, our package use the implementation of creative telescoping (Chyzak's algorithm). However, we also added support for the more recent [MultivariateCreativeTelescoping.jl](https://github.com/HBrochet/MultivariateCreativeTelescoping.jl). You can install it using:

```
julia
using Pkg
Pkg.add("MultivariateCreativeTelescoping")
```

To use it as the default creative telescoping algorithm in the computation, simply set the "use_julia_for_CT" parameter to True:

```
R =  PolynomialRing(QQ, "x", 2); x = R.gens()
volume2([1- x[0]^2 - x[1]^2, 1 - (x[0]-1)^2 - x[1]^2], prec=200, use_julia_for_CT=True)
```

## Using HypersurfaceRegions.jl

You can install it using:

```
julia
using Pkg
Pkg.add("HypersurfaceRegions.jl")
```


## Generate Docs

If you would like to browse the documentation to this package in an html page you can generate it with:

``TODO: Add numpydoc code.``

## Features

- [x] First complete implementation.
- [ ] Generate the docs automatically and host online.
- [x] Volumes as objects: 
      - [x] Store computed PF operators
      - [x] Store computed slice / deformed volumes
      - [ ] Resume computations
      - [ ] Recompute with increased precision without recomputing PF operators
- [x] CT using [MultivariateCreativeTelescoping.jl](https://github.com/HBrochet/MultivariateCreativeTelescoping.jl)
- [ ] Parallelization
