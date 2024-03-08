using JuMP

function plaquetteIsInBounds(Conf::AbstractMatrix, iCenter::Integer, jCenter::Integer)
    i = (iCenter - 1):(iCenter + 1)
    j = (jCenter - 1):(jCenter + 1)
    return checkbounds(Bool, Conf, i, j)
end

function setUpSpiderWeb(optimizer,L;)
    model = Model(optimizer)
    setUpSpiderWeb!(model,L,optimizer())
end

function setUpSpiderWeb!(model,L,opt=nothing)
    INDEX = 1:L
    @variable(model, Sz[INDEX,INDEX], Bin)
    setConstraints!(model,L,opt)
end

function setConstraints!(model,L,opt=nothing)
    Sz = model[:Sz]
    for i in 1:L
        for j in 1:L
            iseven(i+j) || continue
            plaquetteIsInBounds(y, i, j) || continue
            setConstraint!(model,y,i,j,opt)
        end
    end
    return model
end

function setConstraint!(model,Sz,i,j,opt)

    @constraint(model,
    +Sz[i, j + 1]
    +Sz[i - 1, j + 1]
    -Sz[i - 1, j]
    -Sz[i - 1, j - 1]
    +Sz[i, j - 1]
    +Sz[i + 1, j - 1]
    -Sz[i + 1, j]
    -Sz[i + 1, j + 1] == 0)
end


function fixSpins!(model,fixInds,vals)
    
    Sz = model[:Sz]

    for i in eachindex(Sz)
        if is_fixed(Sz[i])
            unfix(Sz[i])
        end
    end
    for (I,v) in zip(fixInds,vals)
        fix(Sz[I], v; force = true)
    end
    return model
end

function getRandomSpins(L,initializeDenominator)
    numFixed = L^2÷initializeDenominator
    fixInds = rand(CartesianIndices((L,L)),numFixed)
    vals = rand(0:1,numFixed)
    return fixInds,vals
end