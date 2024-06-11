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
As this package heavily relies on [Oscar](https://oscar.computeralgebra.de/), it is recommended to install Oscar first ([installation instructions](https://oscar.computeralgebra.de/install/)). Then, install this package via the Julia package manager:
```
] add PBWDeformations
```

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
