# PBWDeformations.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pbwdeformations.github.io/pbwdeformations.jl/stable/) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pbwdeformations.github.io/pbwdeformations.jl/dev/) | [![](https://github.com/PBWDeformations/pbwdeformations.jl/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/PBWDeformations/pbwdeformations.jl/actions/workflows/Tests.yml) [![](https://codecov.io/gh/PBWDeformations/pbwdeformations.jl/branch/master/graph/badge.svg?token=J9XN35I1WU)](https://app.codecov.io/gh/PBWDeformations/pbwdeformations.jl) |

# PBWDeformations.jl Julia package

## Install

To install this package in Julia:
```
using Pkg; Pkg.add("PBWDeformations")
```

## Functionality

The package will provide both a general framework and specialized functions in order to
- classify PBW deformations of certain smash products and
- study their representations.

To solve classification problems efficiently, we use representation theoretic ideas.


## Basic usage

Please consult the [example jupyter notebook](https://nbviewer.org/urls/gitlab.com/johannesflake/pbwdeformations.jl/-/raw/master/examples/PBWDeformationsNotebook.ipynb) for v0.1.
We expect documentation to be found at some point in future at [https://pbwdeformations.github.io/pbwdeformations.jl/](https://pbwdeformations.github.io/pbwdeformations.jl/).

## General Disclaimer

All code in this repository is preliminary work.

It comes with absolutely no warranty and will most likely have errors. If you use it for computations, please check the correctness of the result very carefully.

Also, everything in this repository might change in the future, so currently any update can break the code you wrote upon functionality from packages in this repository.

This software is licensed under the GPL, version 3, or any later version.

## Funding

The development of this Julia package is supported by the Deutsche Forschungsgemeinschaft DFG within the [Collaborative Research Center TRR 195](https://www.computeralgebra.de/sfb/).
