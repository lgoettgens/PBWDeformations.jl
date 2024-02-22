# PBWDeformations

## Introduction
The package [PBWDeformations](https://github.com/PBWDeformations/PBWDeformations.jl) will provide both a general framework and specialized functions in order to
- classify PBW deformations of certain smash products and
- study their representations.

To solve classification problems efficiently, we use representation theoretic ideas.

## Features
- Construct Lie algebras and their modules.
- Construct smash products of the form ``TV \rtimes U(L)`` for a Lie algbra ``L`` and a module ``V``.
- Construct deformations of such smash products.
- Compute a normal form for elements of smash products and their deformations.
- Check, if a given deformation is a PBW-deformation (using [WW14](@cite)).
- For some smash product, compute a basis of all PBW-deformations up to a given degree (using [WW14](@cite)). It is possible to give a basis of the relevant part of the deformation space, which is then used in the computation.
- For some modules of ``\mathfrak{so}_n``, give an explicit basis using arc diagrams or pseudographs (cf. [FM22](@cite)).

## Installation
o install this package in Julia, clone it from github and then run the following command in the Julia REPL from the package directory:
```
using Pkg
Pkg.activate(".")
include(joinpath(pwd(), "etc", "add_oscar.jl"))
using PBWDeformations, Oscar
```

This package depends on a development version of the [Oscar](https://oscar.computeralgebra.de/) package. The `add_oscar.jl` script will add the Oscar package to the current environment. If you want to use the package in a different environment, you can run the `add_oscar.jl` script int the other environment to obtain the specific version of Oscar.

## Outline
```@contents
Pages = [
    "smash_product_lie.md",
    "smash_product_deform_lie.md",
    "smash_product_pbwdeform_lie.md",
    "arc_diagrams.md",
    "util.md",
]
```
- [References](@ref)

### [Index](@id main-index)
```@index
```
