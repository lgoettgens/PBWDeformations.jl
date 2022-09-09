```@meta
CurrentModule = PBWDeformations
```

# Structure constants functions

To be able to use more sophisticated techniques, we need to have a description of Lie algebras in terms of matrices s.t. the standard representation operates by simple matrix-vector-multiplication. This is done by the functions in this section.

Since most of this package only considers structural constants of Lie algebras, modules etc. one needs to additionally save the knowledge about properties of the occurring objects and their bases, e.g. using [`SmashProductLieInfo`](@ref) for [`SmashProductLie`](@ref).

```@autodocs
Modules = [PBWDeformations]
Pages   = ["LieAlgebraStructConsts.jl"]
```
