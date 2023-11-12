function dot(x::AbstractVector, y::AbstractVector)
    sum = 0.
    for i in eachindex(x,y)
        sum += x[i] * y[i]
    end
    return sum
end
struct Point2D
    x::Float64
    y::Float64
end
function dot(x::Point2D, y::Point2D)
    return x.x * y.x + x.y * y.y
end

func(Object,arguments...)
@code_warntype dot([1,2,3],[4,5,6])
##
using StaticArrays
@code_warntype dot(SVector(1,2,3),SVector(0.4,0.5,0.6))
##
x = [] # can be of any type!
push!(x,1)
push!(x,2.0)

@code_warntype dot(x,[1,3])