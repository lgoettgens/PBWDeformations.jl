struct Group
    groupName :: String
    order :: Int64
    permRep :: Vector{String}
end

GroupAlgebra = Group

function Base.:(==)(ga1::GroupAlgebra, ga2::GroupAlgebra)
    (ga1.groupName, ga1.order, ga1.permRep) ==
    (ga2.groupName, ga2.order, ga2.permRep)
end

function Base.show(io::IO, ga::GroupAlgebra)
    println(io, "Group algebra of the finite group ", ga.groupName, " with ", ga.order, " elements")
end


function _groupAlgebra(group, groupName::String) :: QuadraticAlgebra{GroupAlgebra}
    order = GAP.Order(group)

    relTable = Dict{Tuple{BasisElement, BasisElement}, AlgebraElement}()
        # (grp(i), grp(j)) => [(c, [grp(k)])]

    multTable = fromGAP(GAP.MultiplicationTable(group))

    for i in 1:order, j in 1:order
        relTable[(grp(i), grp(j))] = [(1, [grp(multTable[i][j])])]
    end

    permRep = fromGAP(GAP.List(GAP.Elements(group), GAP.String))
    extraData = Group(groupName, order, permRep)
    basis = [grp(i) for i in 1:order] :: Vector{BasisElement}
    return QuadraticAlgebra{GroupAlgebra}(basis, relTable, extraData)
end

function groupAlgebraCyclicGroup(n::Int64) :: QuadraticAlgebra{GroupAlgebra}
    return _groupAlgebra(GAP.CyclicGroup(GAP.IsPermGroup, n), "C"*string(n))
end

function groupAlgebraDihedralGroup(n::Int64) :: QuadraticAlgebra{GroupAlgebra}
    @assert n % 2 == 0
    return _groupAlgebra(GAP.DihedralGroup(GAP.IsPermGroup, n), "D"*string(n))
end

function groupAlgebraDicyclicGroup(n::Int64) :: QuadraticAlgebra{GroupAlgebra}
    @assert n % 4 == 0
    return _groupAlgebra(GAP.DicyclicGroup(GAP.IsPermGroup, n), "Dic"*string(n))
end

function groupAlgebraAlternatingGroup(n::Int64) :: QuadraticAlgebra{GroupAlgebra}
    return _groupAlgebra(GAP.AlternatingGroup(GAP.IsPermGroup, n), "A"*string(n))
end

function groupAlgebraSymmetricGroup(n::Int64) :: QuadraticAlgebra{GroupAlgebra}
return _groupAlgebra(GAP.SymmetricGroup(GAP.IsPermGroup, n), "S"*string(n))
end
