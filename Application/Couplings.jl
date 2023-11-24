using SpinFRGLattices
using SpinFRGLattices.Octochlore:spin
using SpinFRGLattices.StructArrays
using CairoMakie
##

function getCouplingsToS1(RefSite = Rvec(0,0,1))
    L=4
    S1Terms = spin{Int,Rvec_2D}[]
    for n1 in -L+RefSite.n1:L+RefSite.n1,n2 in -L+RefSite.n2:L+RefSite.n2
        isodd(n1+n2) && continue
        sites = [Rvec(n1+1,n2,1),Rvec(n1+1,n2+1,1),Rvec(n1,n2+1,1),Rvec(n1-1,n2+1,1),Rvec(n1-1,n2,1),Rvec(n1-1,n2-1,1),Rvec(n1,n2-1,1),Rvec(n1+1,n2-1,1)]
        sgns = [1,1,-1,-1,1,1,-1,-1]
        plaq = [spin(r,s) for (r,s) in zip(sgns,sites)]
        for s1 in plaq
            for s2 in plaq
                if s1.site == RefSite
                    push!(S1Terms,spin(s1.fac*s2.fac,s2.site))
                elseif s2.site == RefSite
                    push!(S1Terms,spin(s1.fac*s2.fac,s1.site))
                end
            end
        end
    end
    return StructArray(S1Terms)
end

##
using DataFrames

##
p0 = Point2f(1,0)
a = SpinFRGLattices.Octochlore.reduceCouplings(getCouplingsToS1(Rvec(p0...,1)))
filter!(x->x.fac != 0,a)
filter!(x->x.site != Rvec(p0...,1),a)
aDF = DataFrame(;n1 = StructArray(a.site).n1, n2 = StructArray(a.site).n2,fac = a.fac)
##
let
    p0 = Point2f(0,0)
    a = SpinFRGLattices.Octochlore.reduceCouplings(getCouplingsToS1(Rvec(p0...,1)))
    filter!(x->x.fac != 0,a)
    filter!(x->x.site != Rvec(p0...,1),a)

    aDF = DataFrame(;n1 = StructArray(a.site).n1, n2 = StructArray(a.site).n2,fac = a.fac)

    p = [Point2f(n1,n2) for (n1,n2) in zip(aDF.n1,aDF.n2)]
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,xgridcolor = :black,ygridcolor = :black,backgroundcolor = :white)
    col = Dict(0=>:black,2 => :red, -2 => :blue,4 => :darkred, -4 => :darkblue)
    ls(i) = i >0 ? :solid : :dash
    for (i,p) in enumerate(p)
        fac = aDF.fac[i]
        lines!([p0,p],linewidth=5*abs(fac),color=col[fac],linestyle = ls(fac) ,label = string(fac))
    end
    scatter!([Point2f(n1,n2)+p0 for n1 in -2:2 for n2 in -2:2],color = :black,markersize = 10)
    scatter!(p,color = :black,markersize = 50,marker = '×')

    axislegend(ax,unique=true)

    save("couplings_$(Tuple(Int.(p0))).pdf",fig)
    fig
end

##

let
    L = 1
    AllPoints = [Point2f(n1,n2) for n1 in -L:L for n2 in -L:L]
    inbox(x) = abs(x[1]) <= 1 && abs(x[2]) <= 1
    inbox(x::Rvec_2D) = abs(x.n1) <= 1 && abs(x.n2) <= 1
    filter!(inbox,AllPoints)

    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1,xgridcolor = :black,ygridcolor = :black,backgroundcolor = :white,xticks = [-1,0,1],yticks = [-1,0,1])

    PlotPoints = Set(copy(AllPoints))
    for p0 in AllPoints
        a = SpinFRGLattices.Octochlore.reduceCouplings(getCouplingsToS1(Rvec(p0...,1)))
        filter!(x->x.fac != 0,a)
        filter!(x->inbox(x.site),a)
        aDF = DataFrame(;n1 = StructArray(a.site).n1, n2 = StructArray(a.site).n2,fac = a.fac)

        p = [Point2f(n1,n2) for (n1,n2) in zip(aDF.n1,aDF.n2)]

        col = Dict(0=>:black,2 => :red, -2 => :blue,4 => :darkred, -4 => :darkblue)

        for (i,p) in enumerate(p)
            lines!([p0,p],linewidth=1*abs(aDF.fac[i]),color=col[aDF.fac[i]],label = string(aDF.fac[i]))
            push!(PlotPoints,p)
        end
    end
    scatter!(collect(PlotPoints),color = :black,markersize = 25)

    axislegend(ax,unique=true)

    save("couplings_All_plaq.png",fig)
    fig
end
##
using LargeN
using LargeN.StaticArrays
function SpiderWebBasis()
    a1 = SA[1,1]
    a2 = SA[-1,1]
    b = [SA[0,0],SA[1,0]]
    return Basis_Struct_2D(;a1,a2,b,NNdist = 1.,NUnique =2)
end

let 
    B = SpiderWebBasis()
    L = 2
    
    points = [Point2f(getCartesian(Rvec(n1,n2,b),B)...) for n1 in -L:L for n2 in -L:L for b in 1:2]
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    scatter!(points,color = :black,markersize = 10)
    fig
end
##


function getJMatrix(q)
    qx,qy = q
    
    J11 = 2cos(2(qx+qy))+2cos(2(qx-qy))-4cos(2*qx)-4cos(2*qy)
    J12 = +2cos(2qx+qy) + 2cos(qx -2qy) - 2cos(2qx-qy) - 2cos(qx+2qy)
    J21 = J12
    J22 = -4cos(qx+qy)-4cos(qx-qy) +2cos(2qx) + 2cos(2qy)
    J= LargeN.Hermitian(SA[
        J11 J12;
        J21 J22
    ])
    return J
end

LargeN.getNCell(::typeof(getJMatrix),) = 2

##
let 
    res = 30
    q = LinRange(-pi,pi,res)
    eig = LargeN.ComputeEig2D(getJMatrix,res,ext = pi,2)
    fig = Figure()
    ax = Axis3(fig[1,1],aspect = (1,1,1/4))
    surface!(ax,q,q,eig.values[1,:,:],colorrange = extrema(eig.values),transparency = true)
    surface!(ax,q,q,eig.values[2,:,:],colorrange = extrema(eig.values),transparency = true,colormap = (:viridis,0.6))
    save("eig.png",fig)
    fig
end
##

let 
    # chi = getChiFunction(0.0001,getJMatrix,2,2;BZextent = float(2pi),nk=80)
    chi = getZeroTChi(getJMatrix)
    q = LinRange(-pi,pi,700)
    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    chiq = [chi(q1,q2) for q1 in q, q2 in q]
    heatmap!(ax,q,q,chiq)
    save("chi.png",fig)
    fig
end