struct Group
    groupName :: String
    order :: Int64
    permRep :: Vector{String}
end

GroupAlgebra = Group

function Base.:(==)(ga1::GroupAlgebra, ga2::GroupAlgebra) :: Bool
    (ga1.groupName, ga1.order, ga1.permRep) ==
    (ga2.groupName, ga2.order, ga2.permRep)
end

function Base.show(io::IO, ga::GroupAlgebra) :: Nothing
    println(io, "Group algebra of the finite group ", ga.groupName, " with ", ga.order, " elements")
end


function _groupAlgebra(group #= :: Gap.Group =#, groupName::String; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert GAP.IsGroup(group)
    @assert GAP.Order(group) != GAP.infinity
    
    order = GAP.Order(group)

    relTable = Dict{Tuple{BasisElement{C}, BasisElement{C}}, AlgebraElement{C}}()
        # (grp(i), grp(j)) => [(coeff, [grp(k)])]

    multTable = fromGAP(GAP.MultiplicationTable(group))

    for i in 1:order, j in 1:order
        relTable[(grp(i; C), grp(j; C))] = [(C(1), [grp(multTable[i][j]; C)])]
    end

    permRep = fromGAP(GAP.List(GAP.Elements(group), GAP.String))
    extraData = Group(groupName, order, permRep)
    basis = [grp(i; C) for i in 1:order] :: Vector{BasisElement{C}}

    return QuadraticAlgebra{C, GroupAlgebra}(basis, relTable, extraData)
end

function groupAlgebraCyclicGroup(n::Int64; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    return _groupAlgebra(GAP.CyclicGroup(GAP.IsPermGroup, n), "C"*string(n); C)
end

function groupAlgebraDihedralGroup(n::Int64; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert n % 2 == 0
    return _groupAlgebra(GAP.DihedralGroup(GAP.IsPermGroup, n), "D"*string(n); C)
end

function groupAlgebraDicyclicGroup(n::Int64; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert n % 4 == 0
    return _groupAlgebra(GAP.DicyclicGroup(GAP.IsPermGroup, n), "Dic"*string(n); C)
end

function groupAlgebraAlternatingGroup(n::Int64; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    return _groupAlgebra(GAP.AlternatingGroup(GAP.IsPermGroup, n), "A"*string(n); C)
end

function groupAlgebraSymmetricGroup(n::Int64; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, GroupAlgebra}
    return _groupAlgebra(GAP.SymmetricGroup(GAP.IsPermGroup, n), "S"*string(n); C)
end
