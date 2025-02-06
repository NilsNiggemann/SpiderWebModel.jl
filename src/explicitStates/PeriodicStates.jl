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
function periodicStateDenseLoops(L)
    UC = SA[
        0 0 0 0;
        -1 0 1 0;
        0 0 0 0;
        1 0 -1 0;
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

function periodicState6x6Condensate(L)
    UC = SA[
        0  0  -1  0  0  0;
        0  1   0  1  0  1;
        0  0  -1  0  0  0;
        0  0   0  0  0  0;
        0  0  -1  0  0  0;
        0  0   0  0  0  0;
    ]
    getPeriodicState(UC, L, L, 0)
end


SELECTED_S12CONFS = [
    # PeriodicMatrix(Int8[1 1 1 -1; 1 -1 -1 -1],2, 4, 2),
    PeriodicMatrix(Int8[1 1 1 -1; -1 -1 1 -1; 1 -1 1 1; 1 -1 -1 -1],4,4,0),
    PeriodicMatrix(Int8[-1 -1 -1 1 -1 -1; -1 -1 -1 1 1 1; -1 -1 -1 -1 -1 1; -1 1 1 1 -1 1; -1 1 -1 -1 1 -1; 1 1 -1 1 -1 1],6, 6, 6),
    PeriodicMatrix(Int8[-1 -1 -1 1; 1 1 -1 1; -1 1 -1 -1; -1 1 1 1; -1 -1 -1 1],5, 4, 5),
    PeriodicMatrix(Int8[1 1 -1 1 -1 -1; 1 -1 1 1 -1 -1; 1 -1 1 1 -1 1; 1 -1 1 -1 -1 1; -1 1 1 -1 -1 1; -1 1 1 -1 1 1],6, 6, -2),
    PeriodicMatrix(Int8[1 -1 1 -1 1; -1 1 1 -1 1; -1 1 -1 1 1; 1 1 -1 1 -1; 1 -1 1 1 -1],5,5,0),
    PeriodicMatrix(Int8[-1 -1 -1 1; 1 1 -1 1; -1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1; 1 1 -1 1],6, 4, 6),
]
"""Select one of the 28 precomputed inequivalent 4x4 sectors."""
function getSelectedS12PeriodicState(L,i)
    UC = SELECTED_S12CONFS[i]
    getPeriodicState(UC.UC, L, L, UC.offset)
end

const SELECTED_4x4_CONFS = Int8[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;;; 0 -2 0 2; 0 0 0 0; 0 2 0 -2; 0 0 0 0;;; 0 0 0 0; 0 -2 0 0; 0 0 0 0; 0 -2 0 0;;; -2 0 0 0; 2 0 -2 0; -2 0 0 0; -2 0 2 0;;; 0 0 -2 0; 0 0 0 0; 2 0 0 0; 0 0 0 0;;; 0 2 0 -2; 2 0 -2 0; 0 -2 0 2; -2 0 2 0;;; -2 0 -2 0; 0 0 0 0; -2 0 -2 0; 0 0 0 0;;; -2 0 -2 0; 0 2 0 2; -2 0 -2 0; 0 2 0 2;;; 0 -2 0 0; -2 0 0 0; 0 0 0 -2; 0 0 -2 0;;; 0 0 -2 0; 0 2 0 0; 0 0 -2 0; 0 2 0 0;;; -2 2 -2 -2; 0 0 0 0; -2 -2 -2 2; 0 0 0 0;;; 0 0 0 -2; 0 -2 -2 -2; 0 -2 0 0; -2 0 0 0;;; 0 0 -2 2; 2 0 0 0; 0 2 -2 0; 0 0 2 0;;; 0 0 -2 -2; 0 0 -2 0; 0 -2 -2 0; -2 2 0 2;;; -2 0 -2 2; 2 0 0 0; -2 2 -2 0; 0 0 2 0;;; -2 0 -2 0; 0 -2 0 -2; 0 0 0 0; 0 0 0 0;;; -2 -2 -2 0; 0 0 -2 0; -2 0 -2 -2; -2 0 0 0;;; -2 -2 -2 2; 0 -2 0 0; 0 2 0 -2; 0 -2 0 0;;; -2 0 -2 0; 0 -2 0 0; -2 0 -2 0; 0 0 0 2;;; -2 -2 0 0; 0 2 -2 0; -2 0 0 -2; -2 2 0 0;;; 0 -2 0 0; -2 2 0 0; 0 0 0 -2; 0 0 -2 -2;;; 0 0 0 -2; 0 2 -2 0; 0 -2 0 0; -2 0 0 -2;;; 0 0 -2 0; 0 -2 0 -2; 0 0 -2 0; 0 -2 0 -2;;; -2 0 -2 0; -2 -2 2 -2; -2 0 -2 0; 2 0 -2 0;;; 0 -2 -2 0; 0 -2 -2 0; 0 0 -2 -2; -2 -2 0 0;;; -2 2 -2 0; 2 -2 0 0; -2 0 -2 2; 0 0 2 2;;; 0 0 -2 2; 0 -2 2 -2; 0 2 -2 0; 2 -2 0 -2;;; 0 -2 -2 0; -2 -2 0 -2; 0 0 -2 -2; 0 -2 -2 -2]
"""Select one of the 28 precomputed inequivalent 4x4 sectors."""
function get4x4PeriodicState(L,i)
    UC = SELECTED_4x4_CONFS[:,:,i]
    getPeriodicState(UC, L, L, 0)
end
function get4x4PeriodicSpinConf(L,i;kwargs...)
    S0 = get4x4PeriodicState(L,i)
    S = stencilConfig(zeros(L,L),1,boundaryCondition= :periodic,kwargs...)
    S .= S0
    return S
end