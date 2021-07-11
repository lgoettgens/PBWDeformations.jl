using Test
using TestSetExtensions

using Combinatorics
using Oscar
using Random
using SymPy

using PBWDeformations

PD = PBWDeformations

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

lie = PD.lie
mod = PD.mod
grp = PD.grp


numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]
