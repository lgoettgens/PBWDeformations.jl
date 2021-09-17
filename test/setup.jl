using Test
using TestSetExtensions

using Combinatorics
using Oscar
using Random
using SparseArrays

using PBWDeformations

PD = PBWDeformations

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

lie = PD.lie
mod = PD.mod
grp = PD.grp
test = PD.test

ScalarTypes = PD.ScalarTypes
DefaultScalarType = PD.DefaultScalarType

BasisElement = PD.BasisElement
Monomial = PD.Monomial
AlgebraElement = PD.AlgebraElement

unpack = PD.unpack
issamesum = PD.issamesum
normal_form = PD.normal_form
comm = PD.comm
≐ = PD.:(≐)

numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]
indexRange = -20:20

randLength(start=0) = rand(start:10)
randNums(quantity) = rand(indexRange, quantity)
randNum() = randNums(1)[1]

randBasisElement() = test(randNum())

randMonomial(basis) = Monomial{DefaultScalarType}([basis[rand(1:length(basis))] for _ in 1:randLength()])
randMonomial() = randMonomial(test(indexRange))

randAlgebraElement(basis) = AlgebraElement{DefaultScalarType}([(DefaultScalarType(randNum()), randMonomial(basis)) for _ in 1:randLength()])
randAlgebraElement() = randAlgebraElement(test(indexRange))