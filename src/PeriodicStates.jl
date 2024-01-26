function getPeriodicState(UC,Lx,Ly,offset=0)
    Lx_UC,Ly_UC = size(UC)
    Mat = zeros(Float64,Lx,Ly)

    _remapIndices(i,j) = remapIndices(i,j,Lx_UC,Ly_UC,offset)
    
    for i in axes(Mat,1)
        for j in axes(Mat,2)
            x_i,y_j = _remapIndices(i,j)

            Mat[i,j] = UC[x_i,y_j]
        end
    end

    S = SpinConfig(Mat,1/2)
end

function remapIndices(x,y,Lx,Ly,offset)
    xRegion = (x-1)÷Lx
    
    y = y - offset*xRegion
    
    x = (x-1)%Lx+1
    y = (y-1)%Ly+1 
    if y <=0 
        y += Ly
    end
    return x,y
end

function horizontal_flip(UC)
    UCnew = stack([UC[:,i] for i in reverse(axes(UC,2))])
end
function vertical_flip(UC)
    UCnew = stack([UC[i,:] for i in reverse(axes(UC,1))])
end
rotate90(UC) = transpose(horizontal_flip(UC))

function getStairCase(L) 
    UC = SA[
        1 1 1 0;
        0 0 1 0;
        1 0 1 1;
        1 0 0 0;
    ] .- 1/2
    getPeriodicState(UC,L,L)
end

function periodicState5x5(L)
    UC = SA[
        1 0 1 0 1;
        0 1 1 0 1;
        0 1 0 1 1;
        1 1 0 1 0;
        1 0 1 1 0;
    ] .- 1/2 
    getPeriodicState(UC,L,L)
end


function periodicState6x6(L)

    UC = SA[
        0 0 1 1 1 1;
        0 0 0 0 0 1;
        1 1 1 0 0 0;
        0 1 1 1 1 1;
        1 0 0 0 1 1;
        1 1 1 1 0 0;
    ] .- 1/2
    UC = vertical_flip(UC)
    getPeriodicState(UC,L,L,-2)
end

