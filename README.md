# Exact Volumes of Semi-Algebraic Convex Bodies

## Description

We provide a SageMath package to compute the exact volume of semi-algebraic convex bodies. The input consists of $k$ concave polynomials
$`f_1,f_2,\ldots, f_k  \in  \mathbb{Q}[x_1,\ldots,x_n],`$
where concave is defined by
`` TODO definition.
``
The output is a (real/complex) ball of radius $`2^{-\text{prec}}`$ which contains the volume of the convex semi-algebraic set
$`
C = \left\{x \in \mathbb{R}^n \mid f_1(x) > 0\right\} \cap \ldots \cap \left\{x \in \mathbb{R}^n \mid f_k(x) > 0\right\}.
`$

The computation of exact volumes of semi-algebraic convex bodies defined by finitely many concave polynomials is  motivated by geometric statistics, where intersections of convex bodies arise as maximum likelihood estimator (MLE) sets.

This package implements the algorithms described in

  [1] L. Ramesh, N. Weiss (2026): "Exact Volumes of Semi-Algebraic Convex Bodies". [arXiv:to_appear](https://arxiv.org)

which rests on work by [Lairez, Mezzarobba, and Safey El Din](http://dx.doi.org/10.1145/3326229.3326262).

The experimental results of [1] can be reproduced using the Jupyter notebooks:

- Notebook 1
- Notebook 2

This is a project of [Lakshmi Ramesh](https://lax202.github.io) and [Nicolas Weiss](https://nicolas-alexander-weiss.github.io).


## Requirements

This package was developed using SageMath (10.6) and using Python (3.12.5). It requires installation of:

SageMath Packages
- [ore_algebra](https://github.com/mkauers/ore_algebra/)

External Software:
- [Macaulay2](https://github.com/Macaulay2/M2/wiki)
- [msolve](https://msolve.lip6.fr/downloads/msolve-tutorial.pdf)


(Optional) Julia Package:
- [HypersurfaceRegions.jl](https://github.com/JuliaAlgebra/HypersurfaceRegions.jl)


## Installation

Begin by installing the requirements above. 

Instructions for installing the ore_algebra package can be found on their [Github page](https://github.com/mkauers/ore_algebra/). Note that Macaulay2 and msolve have to be accessible from the command line via "M2" and "msolve" respectively.

To install our package, download this repository and then install it via

``sage -pip install <path-to-"exact_convex_volumes"-repo>``

If you would like to avoid installing the package you can simply add the folder "src" to the sys.path.

## Generate Docs

If you would like to browse the documentation to this package in an html page you can generate it with:

``TODO: Add numpydoc code.``

**Features**: 
- [x] First complete implementation.
- [ ] Generate the docs automatically and host online.
- [ ] Volumes as objects: 
      - [ ] Store computed PF operators
      - [ ] Resume computations
      - [ ] Recompute with increased precision without recomputing PF operators
- [ ] Parallelization
