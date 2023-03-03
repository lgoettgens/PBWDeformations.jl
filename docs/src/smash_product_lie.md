```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Smash products

## Constructors

### General case
```@docs
smash_product_lie
```

### Highest weight / GAP case
```@docs
smash_product_lie_highest_weight
```

### Constructive cases

#### ``\mathfrak{so}_n`` -- orthogonal Lie algebra
For the orthogonal Lie algebras ``\mathfrak{so}_n`` there is a different constructor, that results in the well-known basis of ``\mathfrak{so}_n`` given by `x_i_j` ``= E_{i,j} - E_{j,i}`` for ``1 \leq i < j \leq n``. See [`liealgebra_so_basis`](@ref).

```@docs
smash_product_lie_so_symmpowers_standard_module
smash_product_lie_so_extpowers_standard_module
```

#### ``\mathfrak{sp}_{2n}`` -- symplectic Lie algebra
For more details about the basis used for ``\mathfrak{sp}_{2n}``, refer to [`liealgebra_sp_basis`](@ref).

```@docs
smash_product_lie_sp_symmpowers_standard_module
smash_product_lie_sp_extpowers_standard_module
```

## SmashProductLie struct
```@docs
SmashProductLie
SmashProductLieInfo
```

## Functions
The [`SmashProductLie`](@ref) struct can be used as an argument for the following functions:
- `gens`
- `ngens`
- `change_base_ring`
