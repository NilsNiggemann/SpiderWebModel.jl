function horizontal_flip(UC)
    UCnew = stack([UC[:, i] for i in reverse(axes(UC, 2))])
end
function vertical_flip(UC)
    UCnew = stack([UC[i, :] for i in reverse(axes(UC, 1))])
end
rotate90(UC) = transpose(horizontal_flip(UC))

function getStairCase(L)
    UC = SA[
        1 1 1 0
        0 0 1 0
        1 0 1 1
        1 0 0 0
    ] .- 1 / 2
    getPeriodicState(UC, L, L)
end

function periodicState5x5(L)
    UC = SA[
        1 0 1 0 1
        0 1 1 0 1
        0 1 0 1 1
        1 1 0 1 0
        1 0 1 1 0
    ] .- 1 / 2
    getPeriodicState(UC, L, L)
end

function periodicState6x6(L)
    UC = SA[
        0 0 1 1 1 1
        0 0 0 0 0 1
        1 1 1 0 0 0
        0 1 1 1 1 1
        1 0 0 0 1 1
        1 1 1 1 0 0
    ] .- 1 / 2
    UC = vertical_flip(UC)
    getPeriodicState(UC, L, L, -2)
end

function periodicState6x6_3(L)
    UC = SA[
        0 0 0 1 0 0
        0 0 0 1 1 1
        0 0 0 0 0 1
        0 1 1 1 0 1
        0 1 0 0 1 0
        1 1 0 1 0 1
    ] .- 1 / 2
    getPeriodicState(UC, L, L, 0)
end

function periodicStateLoops(L)
    UC = SA[
        1 0 1 0;
        0 -1 0 0;
        0 0 0 0;
        0 -1 0 0;
    ]' .*2
    getPeriodicState(UC, L, L, 0)
end

function periodicStateDiag(L)
    UC = SA[
        0 1 0 -1;
        -1 0 1 0;
        0 -1 0 1;
        1 0 -1 0;
    ] .*2
    getPeriodicState(UC, L, L, 0)
end

function periodicStateWeb(L)
    UC = SA[
        -1 1 1 1;
        -1 1 -1 0;
        1 1 0 1;
        -1 0 -1 0;
    ] .*2
    getPeriodicState(UC, L, L, 0)
end

function periodicPlainWeave(L)
    UC = SA[
        1 0 1 0;
        1 -1 -1 -1;
        1 0 1 0;
        -1 -1 1 -1;
    ] .*2
    getPeriodicState(UC, L, L, 0)
end
