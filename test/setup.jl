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
test = PD.test

BasisElement = PD.BasisElement
Monomial = PD.Monomial
AlgebraElement = PD.AlgebraElement

sameSum = PD.sameSum
normalForm = PD.normalForm
comm = PD.comm
≐ = PD.:(≐)

numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]

randLength(start=0) = rand(start:10)
randNums(quantity) = rand(-20:20, quantity)
randNum() = randNums(1)[1]
randMonomial(basis) = Monomial([basis[rand(1:length(basis))] for _ in 1:randLength()])
randAlgebraElement(basis) = AlgebraElement{Rational{Int64}}([(Rational{Int64}(randNum()), randMonomial(basis)) for _ in 1:randLength()])