# PBWDeformations.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lgoettgens.github.io/PBWDeformations.jl/stable/) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://lgoettgens.github.io/PBWDeformations.jl/dev/) | [![](https://github.com/lgoettgens/PBWDeformations.jl/actions/workflows/OscarTests.yml/badge.svg?branch=master)](https://github.com/lgoettgens/PBWDeformations.jl/actions/workflows/OscarTests.yml) [![](https://codecov.io/gh/lgoettgens/PBWDeformations.jl/branch/master/graph/badge.svg?token=J9XN35I1WU)](https://app.codecov.io/gh/lgoettgens/PBWDeformations.jl) |

# PBWDeformations.jl Julia package

## Install

To install this package in Julia, clone it from github and then run the following command in the Julia REPL from the package directory:
```
using Pkg
Pkg.activate(".")
include(joinpath(pwd(), "etc", "add_oscar.jl"))
using PBWDeformations, Oscar
```

This package depends on a development version of the [Oscar](https://oscar.computeralgebra.de/) package. The `add_oscar.jl` script will add the Oscar package to the current environment. If you want to use the package in a different environment, you can run the `add_oscar.jl` script int the other environment to obtain the specific version of Oscar.

## Functionality

The package will provide both a general framework and specialized functions in order to
- classify PBW deformations of certain smash products and
- study their representations.

To solve classification problems efficiently, we use representation theoretic ideas.


## Basic usage

Please consult the [example jupyter notebook](https://nbviewer.org/urls/raw.githubusercontent.com/lgoettgens/PBWDeformations.jl/master/examples/PBWDeformationsNotebook.ipynb) for v0.1.
We expect documentation to be found at some point in future at [https://lgoettgens.github.io/PBWDeformations.jl/](https://lgoettgens.github.io/PBWDeformations.jl/).

## General Disclaimer

All code in this repository is preliminary work.

It comes with absolutely no warranty and will most likely have errors. If you use it for computations, please check the correctness of the result very carefully.

Also, everything in this repository might change in the future, so currently any update can break the code you wrote upon functionality from packages in this repository.

This software is licensed under the GPL, version 3, or any later version.

## Funding

The development of this Julia package is supported by the Deutsche Forschungsgemeinschaft DFG within the [Collaborative Research Center TRR 195](https://www.computeralgebra.de/sfb/).
