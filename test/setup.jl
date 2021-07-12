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

Coefficient = PD.Coefficient
BasisElement = PD.BasisElement
Monomial{T} = PD.Monomial{T}
LinearCombination{T} = PD.LinearCombination{T}
AlgebraElement = PD.AlgebraElement
algebraElement = PD.algebraElement

sameSum = PD.sameSum
normalForm = PD.normalForm

numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]
