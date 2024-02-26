function HDF5.h5write(filename::AbstractString,key::AbstractString, P::PlaqMapping)
    ks = collect(keys(P.d))
    vs = collect(values(P.d))

    h5open(filename, "cw") do file
        file[key*"/PlaqMapping/i_site"] = getindex.(ks, 1)
        file[key*"/PlaqMapping/j_site"] = getindex.(ks, 2)
        file[key*"/PlaqMapping/plaqNumber"] = vs
    end
end

HDF5.h5write(filename::AbstractString, P::PlaqMapping) = h5write(filename, "", P)

function h5readPlaqMapping(filename,key="")
    h5open(filename, "r") do file
        i = file[key*"/PlaqMapping/i_site"] |> read
        j = file[key*"/PlaqMapping/j_site"] |> read
        plaqNumber = file[key*"/PlaqMapping/plaqNumber"] |> read
        
        keys = collect(zip(i, j))
        return PlaqMapping(OrderedDict([k => v for (k, v) in zip(keys, plaqNumber)]))
    end
end

function saveAllStates(filename::AbstractString,key::AbstractString, AllStates::AbstractVector{<:SBitVector})
    UIntStates = getproperty.(AllStates, :x)
    h5open(filename, "cw") do file
        file[key*"/AllStates",blosc=9] = UIntStates
    end
end

function readAllStates(filename::AbstractString,key::AbstractString="")
    h5open(filename) do file
        UIntStates = file[key*"/AllStates"] |> read
        return SBitVector.(UIntStates)
    end
end

function saveHamiltonian(filename::AbstractString,key::AbstractString, H::SparseMatrixCSC)
    h5open(filename, "cw") do file
        file[key*"/Hamiltonian/colptr",blosc=9] = H.colptr
        file[key*"/Hamiltonian/m"] = H.m
        file[key*"/Hamiltonian/n"] = H.n
        file[key*"/Hamiltonian/nzval",blosc=9] = H.nzval
        file[key*"/Hamiltonian/rowval",blosc=9] = H.rowval
    end
    return H
end
saveHamiltonian(filename::AbstractString,key::AbstractString, H::Symmetric) = saveHamiltonian(filename, key, H.data)

function readHamiltonian(filename::AbstractString,key::AbstractString="")
    h5open(filename) do file
        colptr = file[key*"/Hamiltonian/colptr"] |> read
        m = file[key*"/Hamiltonian/m"] |> read
        n = file[key*"/Hamiltonian/n"] |> read
        nzval = file[key*"/Hamiltonian/nzval"] |> read
        rowval = file[key*"/Hamiltonian/rowval"] |> read
        return Symmetric(SparseMatrixCSC(m, n, colptr, rowval, nzval))
    end
end


function h5saveHilbertSpace(filename::AbstractString,HSpace::HilbertSpace)
    h5write(filename, "HilbertSpace", HSpace.plaqMapping)
    saveHamiltonian(filename, "HilbertSpace", HSpace.H)
    saveAllStates(filename, "HilbertSpace", HSpace.AllStates)
    return nothing
end

function h5readHilbertSpace(filename::AbstractString,key="")
    plaqMapping = h5readPlaqMapping(filename, "HilbertSpace")
    H = readHamiltonian(filename, "HilbertSpace")
    AllStates = readAllStates(filename, "HilbertSpace")
    return HilbertSpace(AllStates, H, plaqMapping)
end