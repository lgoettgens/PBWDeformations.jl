# PBWDeformations

## Introduction
Documentation for [PBWDeformations](https://gitlab.com/johannesflake/PBWDeformations.jl).

## Features
- Construct smash products of Lie algebras and their modules, either for highest weight modules or (only for ``\mathfrak{so}_n`` and ``\mathfrak{sp}_{2n}``) for symmetric and exterior powers of the standard module.
- Construct deformations of such smash products.
- Compute a normal form for elements of smash products and their deformations.
- Check, if a given deformation is a PBW-deformation (using [WW14](@cite)).
- For some smash product, compute a basis of all PBW-deformations up to a given degree (using [WW14](@cite)). It is possible to give a basis of the relevant part of the deformation space, which is then used in the computation.
- For some modules of ``\mathfrak{so}_n``, give an explicit basis using arc diagrams (cf. [FM22](@cite)).

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
    "structure_constants.md",
    "util.md",
]
```
- [References](@ref)

### [Index](@id main-index)
```@index
```
