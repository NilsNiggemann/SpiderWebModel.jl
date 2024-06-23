using CairoMakie
function SqExact(x,y)
    num = cos(x) - cos(y) +2sin(x)sin(y) 
    denom = (cos(x) - cos(y))^2 + (2sin(x)sin(y))^2
    return num^2/(sqrt(denom)+1e-30)
end

k = LinRange(0,2pi,100)
heatmap(k,k,[SqExact(x,y) for x in k, y in k],axis = (;aspect=1))
##
# perturbation
function V(x)
    s1,s2 = x

    S⁺S⁻(s) = 2 + s*(1-s)
    S⁻S⁺(s) = 2 - s*(1+s)

    V = S⁺S⁻(s1)*S⁻S⁺(s2) + S⁻S⁺(s1)*S⁺S⁻(s2)
end

[V(x) for x in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1)]]