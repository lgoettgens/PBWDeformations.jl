const ActingGroupCache = WeakKeyIdDict{LieAlgebraModule, Tuple{PermGroup, GAPGroupHomomorphism{PermGroup, PermGroup}}}()

function acting_group_with_sgn(V::LieAlgebraModule)
    return get!(ActingGroupCache, V) do
        G, gens_and_sgn = _acting_group_with_gens_and_sgn(V)
        @assert sub(G, first.(gens_and_sgn))[1] == G # TODO: remove this check
        codom = symmetric_group(2)
        h = hom(G, codom, [g for (g, sgn) in gens_and_sgn], [signbit(sgn) ? codom[1] : codom[0] for (g, sgn) in gens_and_sgn]; check=true) # TODO: set check=false
        return G, h
    end
end

function _acting_group_with_gens_and_sgn(V::LieAlgebraModule)
    if is_tensor_generator(V)
        G = symmetric_group(1)
        gens_and_sgn = Tuple{PermGroupElem, Int}[]
        return G, gens_and_sgn
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        inner_stuff = [_acting_group_with_gens_and_sgn(inner_mod) for inner_mod in inner_mods]
        inner_groups = first.(inner_stuff)
        inner_gens_and_sgns = last.(inner_stuff)
        G, injs, _ = inner_direct_product(inner_groups; morphisms=true)::Tuple{eltype(inner_groups), Vector{Oscar.GAPGroupHomomorphism{eltype(inner_groups), eltype(inner_groups)}}, Any}
        @assert degree(G) == sum(degree, inner_groups)
        @assert order(G) == prod(order, inner_groups)
        gens_and_sgn = reduce(vcat, [(inj(g), sgn_g) for (inner_gens_and_sgn, inj) in zip(inner_gens_and_sgns, injs) for (g, sgn_g) in inner_gens_and_sgn]; init=Tuple{PermGroupElem, Int}[])::Vector{Tuple{PermGroupElem, Int}}
        return G, gens_and_sgn
    elseif ((fl_ext, inner_mod, power) = _is_exterior_power(V); fl_ext) || ((fl_sym, inner_mod, power) = _is_symmetric_power(V); fl_sym)
        inner_group, inner_gens_and_sgn = _acting_group_with_gens_and_sgn(inner_mod)
        if power == 1
            return inner_group, inner_gens_and_sgn
        end
        top_group = symmetric_group(power)
        W = wreath_product(symmetric_group(degree(inner_group)), top_group) # maybe replace by GAP.Globals.WreathProductImprimitiveAction
        first_inj = canonical_injection(W, 1)::Oscar.GAPGroupHomomorphism{typeof(inner_group), typeof(W)}
        top_inj = canonical_injection(W, power + 1)::Oscar.GAPGroupHomomorphism{typeof(top_group), typeof(W)}
        W_, _ = sub(W, [(first_inj(g) for g in gens(inner_group))..., (top_inj(h) for h in gens(top_group))...])
        iso = isomorphism(PermGroup, W_)
        G = codomain(iso)
        @assert degree(G) == degree(inner_group) * power
        @assert order(G) == order(inner_group)^power * factorial(power)
        gens_and_sgn = let fl_ext = fl_ext, inner_gens_and_sgn = inner_gens_and_sgn
            Tuple{PermGroupElem, Int}[
                ((iso(first_inj(g)), sgn_g) for (g, sgn_g) in inner_gens_and_sgn)...;
                ((iso(top_inj(h)), fl_ext ? sign(h) : 1) for h in gens(top_group))...
            ]
        end::Vector{Tuple{PermGroupElem, Int}}
        return G, gens_and_sgn
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        inner_group, inner_gens_and_sgn = _acting_group_with_gens_and_sgn(inner_mod)
        if power == 1
            return inner_group, inner_gens_and_sgn
        end
        G, injs, _ = inner_direct_product(fill(inner_group, power); morphisms=true)::Tuple{typeof(inner_group), Vector{Oscar.GAPGroupHomomorphism{typeof(inner_group), typeof(inner_group)}}, Any}
        @assert degree(G) == degree(inner_group)*power
        @assert order(G) == order(inner_group)^power
        gens_and_sgn = let inner_gens_and_sgn = inner_gens_and_sgn
            reduce(vcat, [(inj(g), sgn_g) for inj in injs for (g, sgn_g) in inner_gens_and_sgn]; init=Tuple{PermGroupElem, Int}[])#
        end
        return G, gens_and_sgn
    else
        error("Not implemented.")
    end
end
