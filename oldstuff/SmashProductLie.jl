"""
    smash_product_lie_sp_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{sp}_{2n}`` and the
e-th symmetric power of the standard module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_sp_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # sp_2n, e-th symm power of standard module
    symbL = liealgebra_sp_symbols(n)
    scL = liealgebra_sp_struct_const(n)
    symbV = liealgebra_sp_symmpowers_standard_module_symbols(n, e)
    scV = liealgebra_sp_symmpowers_standard_module_struct_const(n, e)

    info = SmashProductLieInfo(dynkin='C', n=n, constructive_basis=true, power_of_std_mod=e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{sp}_{2n}`` and the
e-th exterior power of the fundamental module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # sp_2n, e-th exterior power of standard module
    symbL = liealgebra_sp_symbols(n)
    scL = liealgebra_sp_struct_const(n)
    symbV = liealgebra_sp_extpowers_standard_module_symbols(n, e)
    scV = liealgebra_sp_extpowers_standard_module_struct_const(n, e)

    info = SmashProductLieInfo(dynkin='C', n=n, constructive_basis=true, power_of_std_mod=-e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end
