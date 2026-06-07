"""
Memory-efficient implementation of periodic tiling generation.
Uses lazy evaluation and callbacks to reduce memory footprint at the cost of recomputation.
"""
function constructAllConfigs_periodic(Lx, Ly, PlaquetteList, path = xdirecPath(cld(Lx, 2), cld(Ly, 2)))
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    El = PlaquetteList[begin]
    Mat = fill(NaN, Lx, Ly)
    P = SpinConfig(PeriodicMatrix(Mat,Lx,Ly), El.S)
    # P = SpinConfig(Mat, El.S)

    path = collect(plaquetteIterator(P))

    
    # correctPath!(path, P)
    # P = SpinConfig(PeriodicMatrix(Mat,Lx,Ly), El.S)

    boundary_plaquettes = [
        (i,j) for i in 1:Lx for j in 1:Ly if iseven(i + j) && (i == 1 || i == Lx || j == 1 || j == Ly)
    ]
    # return boundary_plaquettes
    # boundary_plaquettes = [
        # (3,1),
        # (4,2),
    # ]

    function checkBoundary(Conf)
        for I in boundary_plaquettes
            i,j = Tuple(I)
            Pij = getPlaquette(Conf, i, j)
            c = constraint(Pij)
            isnan(c) && continue
            # display(plotApplPlaquettes(Conf))
            if c ≠ 0
                return false
            end
        end

        return true
    end
    AllConfigs_current = [zeros(Int16, length(path))]

    AllConfigs_next = empty(AllConfigs_current)
    iter = 0
    TileListBuffer = collect(eachindex(PlaquetteList))

    while iter < lastindex(path)
        iter += 1
        i, j = path[iter]
        for history_buff in AllConfigs_current
            history = @view history_buff[1:(iter-1)]
            P = reconstructTiling!(P, history, PlaquetteList, path)
            # fulFillsConstraint_nonStrict(P) || continue
            # checkBoundary(P) || continue
            Pij = getPlaquette(P, i, j)
            Tiles = getFittingTiles!(TileListBuffer, Pij, PlaquetteList)
            for Tile in Tiles
                newhistory = copy(history_buff)
                newhistory[iter] = Tile
                push!(AllConfigs_next, newhistory)
            end
        end
        display(stack(AllConfigs_next))
        AllConfigs_current = AllConfigs_next
        AllConfigs_next = empty(AllConfigs_current)
        println(
            "progress: ",
            iter,
            "/",
            length(path),
            " = ",
            round(iter * 100 / length(path), digits = 1),
            "% \tnumber of Configs: ",
            length(AllConfigs_current),
        )
    end
    # AllConfs = [reconstructTiling!(copy(P), history, PlaquetteList, path) for history in AllConfigs_current]
    # AllConfs = fillEmptyStates_periodic(AllConfigs_current, Lx, Ly, PlaquetteList)
    # return AllConfs
    return AllConfigs_current
end

make_periodic(P::SpinConfig) = SpinConfig(PeriodicMatrix(parent(P), size(P)...), P.S)

function fillEmptyStates_periodic(States, Lx, Ly, PlaquetteList)
    path = collect(plaquetteIterator(SpinConfig(PeriodicMatrix(fill(NaN, Lx, Ly), Lx, Ly), PlaquetteList[begin].S)))
    return fillEmptyStates_periodic_path(States, PlaquetteList, path)
end
function fillEmptyStates_periodic_path(States, PlaquetteList, path)
    println(length(path))
    reconst(x) = make_periodic(reconstructTiling_periodic(x, PlaquetteList, path))

    newStates = empty([reconst(first(States))])
    for State in States
        spinState = reconst(State)
        FilledStates = fillEmptyState(spinState)
        # filter!(fulFillsConstraint_nonStrict, FilledStates)
        append!(newStates, FilledStates)
    end
    return newStates
end

function reconstructTiling_periodic(tilingHistory, PlaquetteList, path)
    Lx = maximum(p[1] for p in path) + 1
    Ly = maximum(p[2] for p in path) + 1
    Spin = _getSpinFromPlaquetteList(PlaquetteList)
    P = SpinConfig(PeriodicMatrix(fill(NaN, Lx, Ly), Lx, Ly), Spin)
    path = correctPath(path, P)
    reconstructTiling!(P, tilingHistory, PlaquetteList, path)
end
# ============================================================================
# Iterator-based lazy generation
# ============================================================================

"""
    ConfigIterator(Lx, Ly, PlaquetteList)

Lazily generate all valid tiling configurations one at a time using depth-first search.
Memory usage: O(path_length) instead of O(total_configs).

# Example
```julia
for config in ConfigIterator(Lx, Ly, PlaquetteList)
    state = reconstructTiling_xDirec(Lx, Ly, config, PlaquetteList)
    # Process state...
end
```
"""
struct ConfigIterator
    Lx::Int
    Ly::Int
    PlaquetteList::Vector
    path::Vector{Tuple{Int,Int}}
    P::SpinConfig
end

function ConfigIterator(Lx, Ly, PlaquetteList)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    El = PlaquetteList[begin]
    P = SpinConfig(fill(NaN, Lx, Ly), El.S)
    path = correctPath!(xdirecPath(LPx, LPy), P)
    ConfigIterator(Lx, Ly, PlaquetteList, path, P)
end

function Base.iterate(iter::ConfigIterator, state=nothing)
    if state === nothing
        # Initialize DFS with stack: [(depth, history)]
        history = zeros(Int, length(iter.path))
        stack = [(1, history)]
        state = stack
    else
        stack = state
    end
    
    isempty(stack) && return nothing
    
    TileListBuffer = collect(eachindex(iter.PlaquetteList))
    
    while !isempty(stack)
        depth, history = pop!(stack)
        
        if depth > length(iter.path)
            # Found complete configuration
            return (copy(history), stack)
        end
        
        # Reconstruct tiling up to current depth
        reconstructTiling!(iter.P, @view(history[1:depth-1]), iter.PlaquetteList, iter.path)
        
        i, j = iter.path[depth]
        Pij = getPlaquette(iter.P, i, j)
        Tiles = getFittingTiles!(TileListBuffer, Pij, iter.PlaquetteList)
        
        # Push children onto stack (reverse to maintain left-to-right order)
        for tile_idx in reverse(Tiles)
            new_history = copy(history)
            new_history[depth] = tile_idx
            push!(stack, (depth + 1, new_history))
        end
    end
    
    return nothing
end

Base.IteratorSize(::Type{ConfigIterator}) = Base.SizeUnknown()

# ============================================================================
# Callback-based generation
# ============================================================================

"""
    forEachConfig(f::Function, Lx, Ly, PlaquetteList; verbose=true)

Generate all valid configurations and apply callback function `f` to each.
No memory accumulation - processes configs immediately.

# Arguments
- `f`: Callback function that takes a configuration history as input
- `Lx, Ly`: System dimensions
- `PlaquetteList`: List of allowed plaquette tiles
- `verbose`: Print progress information

# Example
```julia
periodic_configs = SpinConfig[]
forEachConfig(Lx, Ly, PlaquetteList) do config
    state = reconstructTiling_xDirec(Lx, Ly, config, PlaquetteList)
    if isPeriodicTiling(state, Lx, Ly)
        push!(periodic_configs, state)
    end
end
```
"""
function forEachConfig(f::Function, Lx, Ly, PlaquetteList; verbose=true)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    path = xdirecPath(LPx, LPy)
    
    El = PlaquetteList[begin]
    P = SpinConfig(fill(NaN, Lx, Ly), El.S)
    correctPath!(path, P)
    
    history = zeros(Int, length(path))
    TileListBuffer = collect(eachindex(PlaquetteList))
    
    n_configs = Ref(0)
    
    function dfs(depth)
        if depth > length(path)
            # Found complete config - apply callback
            n_configs[] += 1
            f(copy(history))
            return
        end
        
        reconstructTiling!(P, @view(history[1:depth-1]), PlaquetteList, path)
        i, j = path[depth]
        Pij = getPlaquette(P, i, j)
        Tiles = getFittingTiles!(TileListBuffer, Pij, PlaquetteList)
        
        for tile_idx in Tiles
            history[depth] = tile_idx
            dfs(depth + 1)
        end
        
        if verbose && depth == 1
            println("Explored branch $tile_idx at depth 1")
        end
    end
    
    dfs(1)
    
    verbose && println("Total configurations generated: $(n_configs[])")
    return n_configs[]
end

# ============================================================================
# Chunked processing
# ============================================================================

"""
    constructAllConfigs_chunked(Lx, Ly, PlaquetteList; chunk_size=1000, process_fn=identity, verbose=true)

Process configurations in chunks to balance memory usage vs recomputation.

# Arguments
- `chunk_size`: Maximum number of configs to keep in memory at once
- `process_fn`: Function to apply to each complete configuration
- `verbose`: Print progress information

# Returns
Vector of processed results (can be configs or filtered/transformed data)
"""
function constructAllConfigs_chunked(Lx, Ly, PlaquetteList; 
                                     chunk_size=1000, 
                                     process_fn=identity,
                                     verbose=true)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    path = xdirecPath(LPx, LPy)
    
    El = PlaquetteList[begin]
    P = SpinConfig(fill(NaN, Lx, Ly), El.S)
    correctPath!(path, P)
    
    results = []
    current_chunk = [zeros(Int, length(path))]
    TileListBuffer = collect(eachindex(PlaquetteList))
    
    iter = 0
    while iter < length(path)
        iter += 1
        i, j = path[iter]
        
        next_chunk = empty(current_chunk)
        
        for history_buff in current_chunk
            history = @view history_buff[1:(iter-1)]
            reconstructTiling!(P, history, PlaquetteList, path)
            Pij = getPlaquette(P, i, j)
            Tiles = getFittingTiles!(TileListBuffer, Pij, PlaquetteList)
            
            for Tile in Tiles
                newhistory = copy(history_buff)
                newhistory[iter] = Tile
                push!(next_chunk, newhistory)
            end
        end
        
        current_chunk = next_chunk
        
        # Process and clear completed configs
        if iter == length(path) && length(current_chunk) > 0
            # Process current batch
            batch_results = process_fn.(current_chunk)
            append!(results, batch_results)
            empty!(current_chunk)
        end
        
        if verbose
            println(
                "Progress: $iter/$(length(path)) = ",
                round(iter * 100 / length(path), digits=1),
                "% \tConfigs in memory: $(length(current_chunk))"
            )
        end
    end
    
    # Process any remaining configs
    if !isempty(current_chunk)
        append!(results, process_fn.(current_chunk))
    end
    
    return results
end

# ============================================================================
# Optimized periodic state generation
# ============================================================================

"""
    getAllPeriodicStates_lazy(Lx, Ly, PlaquetteList; offset=0, verbose=true)

Memory-efficient version of getAllPeriodicStates.
Only stores configurations that satisfy periodicity constraint.

# Arguments
- `Lx, Ly`: Unit cell dimensions
- `PlaquetteList`: List of allowed plaquette tiles
- `offset`: Periodic boundary offset
- `verbose`: Print progress information

# Returns
Vector of SpinConfig objects that are periodic
"""
function getAllPeriodicStates_lazy(Lx, Ly, PlaquetteList; offset=0, verbose=true)
    periodic_states = SpinConfig[]
    n_checked = Ref(0)
    n_periodic = Ref(0)
    
    forEachConfig(Lx, Ly, PlaquetteList, verbose=false) do history
        n_checked[] += 1
        
        # Only reconstruct and store if periodic
        if isPeriodicTiling(Lx, Ly, history, PlaquetteList, offset)
            state = reconstructTiling_xDirec(Lx, Ly, history, PlaquetteList)
            push!(periodic_states, state)
            n_periodic[] += 1
            
            if verbose && n_periodic[] % 10 == 0
                println("Found $(n_periodic[]) periodic states (checked $(n_checked[]))")
            end
        end
    end
    
    if verbose
        println("\nTotal: $(n_periodic[]) periodic out of $(n_checked[]) configurations")
        println("Filtering ratio: $(round(100 * n_periodic[] / n_checked[], digits=2))%")
    end
    
    return periodic_states
end

"""
    getAllPeriodicStates_lazy_filled(Lx, Ly, PlaquetteList; offset=0, verbose=true)

Like getAllPeriodicStates_lazy but also fills empty (NaN) sites.
"""
function getAllPeriodicStates_lazy_filled(Lx, Ly, PlaquetteList; offset=0, verbose=true)
    periodic_histories = Vector{Int}[]  # Changed from Int[] to Vector{Int}[]
    
    # First collect periodic histories
    forEachConfig(Lx, Ly, PlaquetteList, verbose=false) do history
        if isPeriodicTiling(Lx, Ly, history, PlaquetteList, offset)
            push!(periodic_histories, history)
        end
    end

    if verbose
        println("Found $(length(periodic_histories)) periodic tilings")
        println("Filling empty sites...")
    end
    # Then fill empty states (this is the memory-intensive part)
    filled_states = fillEmptyStates(periodic_histories, Lx, Ly, PlaquetteList)
    # filter!(fulFillsConstraint_nonStrict, filled_states)
    if verbose
        println("Generated $(length(filled_states)) total states after filling")
    end
    
    return filled_states
end

# ============================================================================
# Iterator with filtering
# ============================================================================

"""
    PeriodicConfigIterator(Lx, Ly, PlaquetteList; offset=0)

Iterator that only yields periodic configurations.
Combines ConfigIterator with on-the-fly periodicity checking.

# Example
```julia
for config in PeriodicConfigIterator(Lx, Ly, PlaquetteList)
    state = reconstructTiling_xDirec(Lx, Ly, config, PlaquetteList)
    # state is guaranteed to be periodic
end
```
"""
struct PeriodicConfigIterator
    base_iter::ConfigIterator
    offset::Int
end

function PeriodicConfigIterator(Lx, Ly, PlaquetteList; offset=0)
    PeriodicConfigIterator(ConfigIterator(Lx, Ly, PlaquetteList), offset)
end

function Base.iterate(iter::PeriodicConfigIterator, state=nothing)
    base_state = state
    
    while true
        result = iterate(iter.base_iter, base_state)
        result === nothing && return nothing
        
        config, base_state = result
        
        # Check if periodic
        if isPeriodicTiling(iter.base_iter.Lx, iter.base_iter.Ly, config, 
                           iter.base_iter.PlaquetteList, iter.offset)
            return (config, base_state)
        end
    end
end

Base.IteratorSize(::Type{PeriodicConfigIterator}) = Base.SizeUnknown()

# ============================================================================
# Utility functions
# ============================================================================

"""
    countConfigs(Lx, Ly, PlaquetteList; periodic_only=false, offset=0)

Count total number of configurations without storing them.
"""
function countConfigs(Lx, Ly, PlaquetteList; periodic_only=false, offset=0)
    if periodic_only
        count = 0
        for _ in PeriodicConfigIterator(Lx, Ly, PlaquetteList, offset=offset)
            count += 1
        end
        return count
    else
        return forEachConfig(Lx, Ly, PlaquetteList, verbose=false) do _
            nothing
        end
    end
end

"""
    estimateMemoryUsage(Lx, Ly, n_configs)

Estimate memory usage for storing n_configs configurations.
"""
function estimateMemoryUsage(Lx, Ly, n_configs)
    LPx = cld(Lx, 2)
    LPy = cld(Ly, 2)
    config_length = LPx * LPy
    bytes_per_config = sizeof(Int) * config_length
    total_bytes = bytes_per_config * n_configs
    
    # Convert to human-readable
    if total_bytes < 1024
        return "$total_bytes B"
    elseif total_bytes < 1024^2
        return "$(round(total_bytes / 1024, digits=2)) KB"
    elseif total_bytes < 1024^3
        return "$(round(total_bytes / 1024^2, digits=2)) MB"
    else
        return "$(round(total_bytes / 1024^3, digits=2)) GB"
    end
end
