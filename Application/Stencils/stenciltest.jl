import SpiderWebModel as SW
using SpiderWebModel.Stencils
##
plaquette = Moore(1)
##
mat = zeros(Int,7,7)
A = StencilArray(
    mat,
    plaquette;
    boundary = Wrap(),
    padding = Halo{:in}()
    )
A .= collect(LinearIndices(A))
Stencils.update_boundary!(A)
##
@time stencil(A, (5, 5))