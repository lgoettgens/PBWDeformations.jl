using Test
using TestSetExtensions

using Combinatorics
using Random
using SymPy

using PBWDeformations

PD = PBWDeformations

lie = PD.lie
mod = PD.mod
grp = PD.grp


numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]
