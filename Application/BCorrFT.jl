using Cubature, CairoMakie, StaticArrays
function ω(kx,ky)
    sx,cx = sincos(kx)
    sy,cy = sincos(ky)
    w2 = (cx - cy)^2 + 4*(sx*sy)^2
    return sqrt(w2)
end
ω(k) = ω(k[1],k[2])
function getBBCorrFieldTheory(L,kx,ky)
    f(k) = ω(k)*ω(SA[kx,ky]-k)
    return hcubature(f, SA[-pi/2,-pi], SA[pi/2,pi], reltol = 1e-6, abstol = 1e-6)[1]
end

let 
    k = LinRange(-pi,pi,100)

    BB = fetch.([Threads.@spawn getBBCorrFieldTheory(10,ki,kj) for ki in k, kj in k])

    heatmap(k,k,BB)
    
end