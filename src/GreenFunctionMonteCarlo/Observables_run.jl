
function saveObservables!(Observables::GFMCObservables_StructureFac_1,i,Walkers)
    (;SqBuffer,StructureFactors,TotalWeights) = Observables
    chunks = ChunkSplitters.chunks(eachindex(Walkers), n = Threads.nthreads())

    Threads.@threads for (i_chunk,αinds) in enumerate(chunks)
        SqFFT = SqFFT(size(get_config(Walkers[begin])))
        for α in αinds
            Config = Walkers[α]
            SqChunk = @view SqBuffer[:,α]
            SqFFT(SqChunk,Config)
        end
    end
end