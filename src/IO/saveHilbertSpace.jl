function HDF5.h5write(filename, P::PlaqMapping)
    ks = collect(keys(P.d))
    vs = collect(values(P.d))

    h5open(filename, "w") do file
        file["PlaqMapping/i_site"] = getindex.(ks, 1)
        file["PlaqMapping/j_site"] = getindex.(ks, 2)
        file["PlaqMapping/plaqNumber"] = vs
    end
end

function h5readPlaqMapping(filename)
    h5open(filename, "r") do file
        i = file["PlaqMapping/i_site"] |> read
        j = file["PlaqMapping/j_site"] |> read
        keys = collect(zip(i, j))
        plaqNumber = file["PlaqMapping/plaqNumber"] |> read
        return PlaqMapping(OrderedDict([k => v for (k, v) in zip(keys, plaqNumber)]))
    end
end

# function h5saveHilbertSpace(filename::AbstractString,HSpace::HilbertSpace)
#     h5open(filename, "w") do file
#         h5write(file, "HilbertSpace/Dimension", HSpace.Dimension)
#         h5write(file, "HilbertSpace/NumberOfSites", HSpace.NumberOfSites)
#         h5write(file, "HilbertSpace/NumberOfParticles", HSpace.NumberOfParticles)
#         h5write(file, "HilbertSpace/NumberOfStates", HSpace.NumberOfStates)
#         h5write(file, "HilbertSpace/States", HSpace.States)
#     end
# end
