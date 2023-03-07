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

    sp, basis = smash_product_lie(coeff_ring, symbL, symbV, scL, scV)

    set_attribute!(sp, :dynkin, 'C')
    set_attribute!(sp, :n, n)
    set_attribute!(sp, :constructive_basis, true)
    set_attribute!(sp, :power_of_std_mod, e)

    return sp, basis
end

"""
    smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{sp}_{2n}`` and the
e-th exterior power of the standard module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # sp_2n, e-th exterior power of standard module
    symbL = liealgebra_sp_symbols(n)
    scL = liealgebra_sp_struct_const(n)
    symbV = liealgebra_sp_extpowers_standard_module_symbols(n, e)
    scV = liealgebra_sp_extpowers_standard_module_struct_const(n, e)

    sp, basis = smash_product_lie(coeff_ring, symbL, symbV, scL, scV)

    set_attribute!(sp, :dynkin, 'C')
    set_attribute!(sp, :n, n)
    set_attribute!(sp, :constructive_basis, true)
    set_attribute!(sp, :power_of_std_mod, -e)

    return sp, basis
end
