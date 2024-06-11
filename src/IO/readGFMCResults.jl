function readResults(filename,binsize)
    energies_raw = h5read(filename,"energies")
    TotalWeights_raw = h5read(filename,"TotalWeights")
    reconfTable_raw = readMMapArray(filename,"reconfigurationTable")
    SaveConfigs_raw = readMMapArray(filename,"SaveConfigs")
    nBra = h5read(filename,"nBranch")
    getrange(i) = i*binsize+1:(i+1)*binsize
    @views function getRes(range)
        energies = energies_raw[range]
        TotalWeights = TotalWeights_raw[range]
        reconfTable = reconfTable_raw[:,range]
        SaveConfigs = SaveConfigs_raw[:,:,:,range]
        return (;energies,TotalWeights,SaveConfigs,reconfTable,nBra)
    end
    return [getRes(getrange(i)) for i in 0:length(energies_raw)÷binsize-1]
end