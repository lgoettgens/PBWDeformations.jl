using PBWDeformations
using Oscar


function main(args)
    @req length(args) == 3 "Usage: julia DBComputeEntry.jl <V_id> <n> <maxdeg>"

    V_id = parse(Int, args[1])
    n = parse(Int, args[2])
    maxdeg = parse(Int, args[3])

    Oscar.set_verbosity_level(:PBWDeformationsDatabase, 1)

    @vprint :PBWDeformationsDatabase "Constructing Lie algebra..."
    L = general_linear_lie_algebra(QQ, n)
    @vprintln :PBWDeformationsDatabase " Done"

    @vprint :PBWDeformationsDatabase "Constructing module..."
    V = module_by_id(L, V_id)
    @vprintln :PBWDeformationsDatabase " Done"

    @vprint :PBWDeformationsDatabase "Constructing smash product..."
    sp = smash_product(L, V)
    @vprintln :PBWDeformationsDatabase " Done"

    PBWDeformations.Database.compute_and_save_instance("db", sp, maxdeg)
end

function module_by_id(L::LinearLieAlgebra, id::Int)
    if id == 100
        return direct_sum(standard_module(L), dual(standard_module(L)))

    elseif id == 112
        return exterior_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 113
        return exterior_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 114
        return exterior_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 115
        return exterior_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 122
        return symmetric_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 123
        return symmetric_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 124
        return symmetric_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 125
        return symmetric_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 132
        return tensor_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 133
        return tensor_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 134
        return tensor_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 135
        return tensor_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 142
        return direct_sum(exterior_power_obj(standard_module(L), 2), exterior_power_obj(dual(standard_module(L)), 2))
    elseif id == 143
        return direct_sum(exterior_power_obj(standard_module(L), 3), exterior_power_obj(dual(standard_module(L)), 3))
    elseif id == 144
        return direct_sum(exterior_power_obj(standard_module(L), 4), exterior_power_obj(dual(standard_module(L)), 4))
    elseif id == 145
        return direct_sum(exterior_power_obj(standard_module(L), 5), exterior_power_obj(dual(standard_module(L)), 5))

    elseif id == 152
        return direct_sum(symmetric_power_obj(standard_module(L), 2), symmetric_power_obj(dual(standard_module(L)), 2))
    elseif id == 153
        return direct_sum(symmetric_power_obj(standard_module(L), 3), symmetric_power_obj(dual(standard_module(L)), 3))
    elseif id == 154
        return direct_sum(symmetric_power_obj(standard_module(L), 4), symmetric_power_obj(dual(standard_module(L)), 4))
    elseif id == 155
        return direct_sum(symmetric_power_obj(standard_module(L), 5), symmetric_power_obj(dual(standard_module(L)), 5))

    elseif id == 162
        return direct_sum(tensor_power_obj(standard_module(L), 2), tensor_power_obj(dual(standard_module(L)), 2))
    elseif id == 163
        return direct_sum(tensor_power_obj(standard_module(L), 3), tensor_power_obj(dual(standard_module(L)), 3))
    elseif id == 164
        return direct_sum(tensor_power_obj(standard_module(L), 4), tensor_power_obj(dual(standard_module(L)), 4))
    elseif id == 165
        return direct_sum(tensor_power_obj(standard_module(L), 5), tensor_power_obj(dual(standard_module(L)), 5))


    elseif id == 200
        return tensor_product(standard_module(L), dual(standard_module(L)))

    elseif id == 212
        return exterior_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 213
        return exterior_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 214
        return exterior_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 215
        return exterior_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 222
        return symmetric_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 223
        return symmetric_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 224
        return symmetric_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 225
        return symmetric_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 232
        return tensor_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 2)
    elseif id == 233
        return tensor_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 3)
    elseif id == 234
        return tensor_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 4)
    elseif id == 235
        return tensor_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 5)

    elseif id == 242
        return tensor_product(exterior_power_obj(standard_module(L), 2), exterior_power_obj(dual(standard_module(L)), 2))
    elseif id == 243
        return tensor_product(exterior_power_obj(standard_module(L), 3), exterior_power_obj(dual(standard_module(L)), 3))
    elseif id == 244
        return tensor_product(exterior_power_obj(standard_module(L), 4), exterior_power_obj(dual(standard_module(L)), 4))
    elseif id == 245
        return tensor_product(exterior_power_obj(standard_module(L), 5), exterior_power_obj(dual(standard_module(L)), 5))

    elseif id == 252
        return tensor_product(symmetric_power_obj(standard_module(L), 2), symmetric_power_obj(dual(standard_module(L)), 2))
    elseif id == 253
        return tensor_product(symmetric_power_obj(standard_module(L), 3), symmetric_power_obj(dual(standard_module(L)), 3))
    elseif id == 254
        return tensor_product(symmetric_power_obj(standard_module(L), 4), symmetric_power_obj(dual(standard_module(L)), 4))
    elseif id == 255
        return tensor_product(symmetric_power_obj(standard_module(L), 5), symmetric_power_obj(dual(standard_module(L)), 5))

    elseif id == 262
        return tensor_product(tensor_power_obj(standard_module(L), 2), tensor_power_obj(dual(standard_module(L)), 2))
    elseif id == 263
        return tensor_product(tensor_power_obj(standard_module(L), 3), tensor_power_obj(dual(standard_module(L)), 3))
    elseif id == 264
        return tensor_product(tensor_power_obj(standard_module(L), 4), tensor_power_obj(dual(standard_module(L)), 4))
    elseif id == 265
        return tensor_product(tensor_power_obj(standard_module(L), 5), tensor_power_obj(dual(standard_module(L)), 5))

    else
        error("Unknown module id: $id")
    end
end

main(ARGS)
