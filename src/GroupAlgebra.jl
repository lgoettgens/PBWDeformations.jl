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


function group_algebra(group #= :: Gap.Group =#, groupName::String; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert GAP.IsGroup(group)
    @assert GAP.Order(group) != GAP.infinity
    
    order = GAP.Order(group)

    relTable = Dict{Tuple{BasisElement{C}, BasisElement{C}}, AlgebraElement{C}}()
        # (grp(i), grp(j)) => [(coeff, [grp(k)])]

    multTable = fromGAP(GAP.MultiplicationTable(group))

    for i in 1:order, j in 1:order
        relTable[(grp(i; C), grp(j; C))] = AlgebraElement{C}([(one(C), Monomial{C}(grp(multTable[i][j]; C)))])
    end

    permRep = fromGAP(GAP.List(GAP.Elements(group), GAP.String))
    extraData = Group(groupName, order, permRep)
    basis = [grp(i; C) for i in 1:order] :: Vector{BasisElement{C}}

    return QuadraticAlgebra{C, GroupAlgebra}(basis, relTable, extraData)
end

function group_algebra_cyclic_group(n::Int64; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    return group_algebra(GAP.CyclicGroup(GAP.IsPermGroup, n), "C"*string(n); C)
end

function group_algebra_dihedral_group(n::Int64; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert n % 2 == 0
    return group_algebra(GAP.DihedralGroup(GAP.IsPermGroup, n), "D"*string(n); C)
end

function group_algebra_dicyclic_group(n::Int64; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    @assert n % 4 == 0
    return group_algebra(GAP.DicyclicGroup(GAP.IsPermGroup, n), "Dic"*string(n); C)
end

function group_algebra_alternating_group(n::Int64; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    return group_algebra(GAP.AlternatingGroup(GAP.IsPermGroup, n), "A"*string(n); C)
end

function group_algebra_symmetric_group(n::Int64; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, GroupAlgebra}
    return group_algebra(GAP.SymmetricGroup(GAP.IsPermGroup, n), "S"*string(n); C)
end
