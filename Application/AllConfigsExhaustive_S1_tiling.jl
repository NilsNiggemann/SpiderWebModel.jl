using Pkg
Pkg.activate(".")
import SpiderWebModel as SW
using CairoMakie
using SpiderWebModel.HDF5

"""
    generate_all_layer_configs(Lx, Ly, PlaquetteList)

Generate all possible configurations for a layer with dimensions Lx × Ly.
Returns vector of configuration arrays.
"""
function generate_all_layer_configs(Lx, Ly, PlaquetteList)
    println("Generating all $(Lx)×$(Ly) configurations...")
    
    @time raw_configs = SW.constructAllConfigs(Lx, Ly, PlaquetteList)
    println("Generated $(length(raw_configs)) raw configurations")
    
    @time filled_configs = SW.stencilConfig.(SW.fillEmptyStates(raw_configs, Lx, Ly, PlaquetteList))
    config_arrays = parent.(parent.(filled_configs))
    println("After filling: $(length(config_arrays)) configurations")
    
    return config_arrays
end

"""
    check_x_boundary_constraints(config, S=1.0)

Check if constraints are satisfied at x-direction boundaries (periodic wrapping).
"""
function check_x_boundary_constraints(config::AbstractMatrix, S=1.0)
    Lx, Ly = size(config)
    
    # Create a periodic version by wrapping
    periodic_config = SW.stencilConfig(config, S, boundaryCondition=:periodic)
    
    # Check plaquettes that span the x-boundary (wrap from Lx to 1)
    for j in 1:Ly
        for i in [Lx, 1]  # Check plaquettes at boundary
            if iseven(i + j)  # Only even plaquettes matter
                P = SW.getPlaquette(periodic_config, i, j)
                if SW.constraint(P) ≠ 0
                    return false
                end
            end
        end
    end
    
    return true
end

"""
    filter_x_boundary_valid_configs(configs, S=1.0)

Filter configurations to keep only those that satisfy x-boundary constraints.
"""
function filter_x_boundary_valid_configs(configs::Vector, S=1.0)
    println("\nFiltering for x-boundary constraint satisfaction...")
    
    x_boundary_valid = filter(configs) do config
        check_x_boundary_constraints(config, S)
    end
    
    println("Found $(length(x_boundary_valid)) configs with valid x-boundaries")
    return x_boundary_valid
end

"""
    check_vertical_compatibility(bottom_layer, top_layer, buffer_config, boundary_plaquettes)

Check if two layers can be stacked vertically while satisfying constraints at the boundary.
"""
function check_vertical_compatibility(bottom_layer::AbstractMatrix, top_layer::AbstractMatrix, 
                                     buffer_config, boundary_plaquettes::Vector)
    Lx, Ly = size(bottom_layer)
    buffer_parent = parent(parent(buffer_config))
    
    # Place layers in buffer
    buffer_parent[:, 1:Ly] .= bottom_layer
    buffer_parent[:, Ly+1:2Ly] .= top_layer
    
    # Check boundary plaquettes
    for (i, j) in boundary_plaquettes
        P = SW.getPlaquette(buffer_config, i, j)
        if SW.constraint(P) ≠ 0
            return false
        end
    end
    
    return true
end

"""
    find_valid_6x4_stackings(x_boundary_valid_configs, Lx, Ly, S=1.0)

Find all valid ways to stack two 6×2 layers to form 6×4 configurations.
Returns:
- valid_stackings: Vector of tuples (bottom_idx, top_idx)
- bottom_to_top_dict: Dict mapping each bottom index to vector of valid top indices
"""
function find_valid_6x4_stackings(x_boundary_valid_configs::Vector, Lx, Ly, S=1.0)
    println("\nFinding valid 6×4 stackings (two 6×2 layers)...")
    
    # Create buffer for testing
    buffer_config = SW.stencilConfig(zeros(Lx, 2*Ly), S, boundaryCondition=:periodic)
    
    # Get boundary plaquettes between the two layers (at y=Ly boundary)
    boundary_plaquettes = Tuple{Int,Int}[]
    for i in 1:Lx
        for j in [Ly, Ly+1]
            if iseven(i + j)  # Only check even plaquettes
                push!(boundary_plaquettes, (i, j))
            end
        end
    end
    boundary_plaquettes = unique(boundary_plaquettes)
    
    println("Checking $(length(boundary_plaquettes)) boundary plaquettes")
    
    # Find all valid stackings
    valid_stackings = Tuple{Int,Int}[]
    bottom_to_top_dict = Dict{Int, Vector{Int}}()
    
    n_configs = length(x_boundary_valid_configs)
    println("Testing $(n_configs^2) combinations...")
    
    for (bottom_idx, bottom_layer) in enumerate(x_boundary_valid_configs)
        valid_tops = Int[]
        
        for (top_idx, top_layer) in enumerate(x_boundary_valid_configs)
            if check_vertical_compatibility(bottom_layer, top_layer, buffer_config, boundary_plaquettes)
                push!(valid_stackings, (bottom_idx, top_idx))
                push!(valid_tops, top_idx)
            end
        end
        
        bottom_to_top_dict[bottom_idx] = valid_tops
        
        if bottom_idx % 100 == 0
            println("Progress: $(bottom_idx)/$(n_configs) - Found $(length(valid_tops)) valid tops for bottom $bottom_idx")
            flush(stdout)
        end
    end
    
    println("\nFound $(length(valid_stackings)) valid 6×4 stackings")
    return valid_stackings, bottom_to_top_dict
end

"""
    check_three_layer_compatibility(layer1, layer2, layer3, buffer_config, boundary_12, boundary_23)

Check if three layers can be stacked while satisfying all internal constraints.
"""
function check_three_layer_compatibility(layer1::AbstractMatrix, layer2::AbstractMatrix, layer3::AbstractMatrix,
                                        buffer_config, boundary_12::Vector, boundary_23::Vector)
    Lx, Ly = size(layer1)
    buffer_parent = parent(parent(buffer_config))
    
    # Place all three layers
    buffer_parent[:, 1:Ly] .= layer1
    buffer_parent[:, Ly+1:2Ly] .= layer2
    buffer_parent[:, 2Ly+1:3Ly] .= layer3
    
    # Check boundary 1-2 (already checked, but include for completeness)
    for (i, j) in boundary_12
        P = SW.getPlaquette(buffer_config, i, j)
        if SW.constraint(P) ≠ 0
            return false
        end
    end
    
    # Check boundary 2-3
    for (i, j) in boundary_23
        P = SW.getPlaquette(buffer_config, i, j)
        if SW.constraint(P) ≠ 0
            return false
        end
    end
    
    return true
end

"""
    check_y_boundary_constraints(layer_bottom, layer_top, buffer_config, boundary_plaquettes)

Check if constraints are satisfied at y-direction boundary (layer_top wraps to layer_bottom).
"""
function check_y_boundary_constraints(layer_bottom::AbstractMatrix, layer_top::AbstractMatrix,
                                     buffer_config, boundary_plaquettes::Vector)
    Lx, Ly = size(layer_bottom)
    buffer_parent = parent(parent(buffer_config))
    
    # Place top layer followed by bottom layer to simulate y-periodic wrapping
    buffer_parent[:, 1:Ly] .= layer_top    # Top layer (at y = 3Ly in full config, wraps to y=0)
    buffer_parent[:, Ly+1:2Ly] .= layer_bottom  # Bottom layer (at y=0 in full config)
    
    # Check boundary plaquettes where y wraps around
    for (i, j) in boundary_plaquettes
        P = SW.getPlaquette(buffer_config, i, j)
        if SW.constraint(P) ≠ 0
            return false
        end
    end
    
    return true
end

"""
    find_all_6x6_configs_layered(; S=1.0, save_output=true, plot_examples=true)

Find all valid 6×6 periodic configurations using optimized layered approach.

Steps:
1. Generate all 6×2 configurations with valid x-boundary constraints
2. Find all valid 6×4 stackings and create lookup dictionary
3. Iterate through 6×4 stackings and find valid third layers with y-boundary constraints
"""
function find_all_6x6_configs_layered(; S=1.0, save_output=true, plot_examples=true)
    println("="^70)
    println("Finding all 6×6 configurations using optimized layered approach")
    println("="^70)
    
    Lx = 6
    Ly = 2
    ALLGS_1 = SW.getAllGS_noMissing(1)
    
    # Step 1: Generate all x-boundary-valid 6×2 configurations
    println("\nStep 1: Generate all 6×2 configurations with valid x-boundaries")
    all_layer_configs = generate_all_layer_configs(Lx, Ly, ALLGS_1)
    x_boundary_valid_configs = filter_x_boundary_valid_configs(all_layer_configs, S)
    
    # Step 2: Find all valid 6×4 stackings
    println("\nStep 2: Find all valid 6×4 stackings")
    valid_6x4_stackings, bottom_to_top_dict = find_valid_6x4_stackings(x_boundary_valid_configs, Lx, Ly, S)
    
    # Step 3: Build 6×6 configurations by adding third layer
    println("\nStep 3: Building 6×6 configurations with y-boundary constraints")
    
    # Create buffers
    buffer_3layer = SW.stencilConfig(zeros(Lx, 3*Ly), S, boundaryCondition=:periodic)
    buffer_2layer = SW.stencilConfig(zeros(Lx, 2*Ly), S, boundaryCondition=:periodic)
    
    # Get boundary plaquettes
    boundary_12 = Tuple{Int,Int}[]
    boundary_23 = Tuple{Int,Int}[]
    boundary_y = Tuple{Int,Int}[]
    
    for i in 1:Lx
        # Boundary between layer 1 and 2 (at y=Ly)
        for j in [Ly, Ly+1]
            if iseven(i + j)
                push!(boundary_12, (i, j))
            end
        end
        # Boundary between layer 2 and 3 (at y=2*Ly)
        for j in [2*Ly, 2*Ly+1]
            if iseven(i + j)
                push!(boundary_23, (i, j))
            end
        end
    end
    boundary_12 = unique(boundary_12)
    boundary_23 = unique(boundary_23)
    
    # For y-boundary check: boundary between layer 3 and layer 1 (wrapping)
    # This has the same structure as boundary_12 but checks wrapping
    boundary_y = copy(boundary_12)
    
    println("Boundary 1-2: $(length(boundary_12)) plaquettes")
    println("Boundary 2-3: $(length(boundary_23)) plaquettes")
    println("Boundary y-wrap: $(length(boundary_y)) plaquettes")
    
    # Iterate through all valid 6×4 stackings
    valid_6x6_configs = []
    n_checked = 0
    n_total = length(valid_6x4_stackings)
    
    println("\nIterating through $(n_total) valid 6×4 stackings...")
    
    for (idx, (bottom_idx, middle_idx)) in enumerate(valid_6x4_stackings)
        layer1 = x_boundary_valid_configs[bottom_idx]
        layer2 = x_boundary_valid_configs[middle_idx]
        
        # Get valid top layers for middle layer (from 6×4 stacking dictionary)
        valid_top_indices = bottom_to_top_dict[middle_idx]
        
        # Check each potential top layer
        for top_idx in valid_top_indices
            layer3 = x_boundary_valid_configs[top_idx]
            
            # Check if all three layers are compatible (internal boundaries)
            if check_three_layer_compatibility(layer1, layer2, layer3, 
                                              buffer_3layer, boundary_12, boundary_23)
                # Check y-boundary constraint (layer3 wraps to layer1)
                if check_y_boundary_constraints(layer1, layer3, buffer_2layer, boundary_y)
                    # Valid configuration found!
                    full_config = zeros(Lx, 3*Ly)
                    full_config[:, 1:Ly] .= layer1
                    full_config[:, Ly+1:2Ly] .= layer2
                    full_config[:, 2Ly+1:3Ly] .= layer3
                    push!(valid_6x6_configs, full_config)
                end
            end
            
            n_checked += 1
        end
        
        # Progress reporting
        if idx % 100 == 0
            progress = 100 * idx / n_total
            println("Progress: $(round(progress, digits=1))% - Found $(length(valid_6x6_configs)) valid 6×6 configs")
            flush(stdout)
        end
    end
    
    # Results
    println("\n" * "="^70)
    println("RESULTS")
    println("="^70)
    println("Total 6×6 configurations found: $(length(valid_6x6_configs))")
    println("="^70)
    return valid_6x6_configs
    # Save results
    if save_output && !isempty(valid_6x6_configs)
        save_configs_to_hdf5(valid_6x6_configs, "all_6x6_configs.h5", Lx, 3*Ly, S)
    end
    
    # Verify and plot
    if plot_examples && !isempty(valid_6x6_configs)
        verify_and_plot_configs(valid_6x6_configs, Lx, 3*Ly, S)
    end
    
    return valid_6x6_configs
end

"""
    save_configs_to_hdf5(configs, filename, Lx, Ly, S)

Save configurations to HDF5 file.
"""
function save_configs_to_hdf5(configs::Vector, filename::String, Lx, Ly, S)
    println("\nSaving results to $filename...")
    
    h5open(filename, "w") do file
        file["n_configs"] = length(configs)
        file["Lx"] = Lx
        file["Ly"] = Ly
        file["S"] = S
        
        for (idx, config) in enumerate(configs)
            file["config_$idx"] = config
        end
    end
    
    println("Saved $(length(configs)) configurations")
end

"""
    verify_and_plot_configs(configs, Lx, Ly, S; n_plot=3)

Verify configurations and plot examples.
"""
function verify_and_plot_configs(configs::Vector, Lx, Ly, S; n_plot=3)
    println("\nVerifying configurations...")
    
    # Convert to SpinConfig objects
    spinconfigs = [SW.stencilConfig(config, S, boundaryCondition=:periodic) for config in configs]
    
    # Verify all satisfy constraints
    n_verified = count(SW.fulFillsConstraint, spinconfigs)
    println("Verification: $n_verified / $(length(spinconfigs)) satisfy constraints")
    
    # Plot examples
    if !isempty(spinconfigs)
        println("\nPlotting first few configurations...")
        fig = Figure(resolution=(400*n_plot, 400))
        n_to_plot = min(n_plot, length(spinconfigs))
        
        for i in 1:n_to_plot
            ax = Axis(fig[1, i], aspect=DataAspect(), title="Config $i")
            SW.plotSpinConfig!(ax, spinconfigs[i])
        end
        
        plot_file = "example_$(Lx)x$(Ly)_configs.png"
        save(plot_file, fig)
        println("Saved example plots to: $plot_file")
    end
    
    return spinconfigs
end

# Run the main function
find_all_6x6_configs_layered()
##
import SpiderWebModel as SW
ALLGS_1 = SW.getAllGS_noMissing(1)
# ALLGS_1 = SW.getAllGS(1)
Lx = 4
Ly = 4
@time a = SW.constructAllConfigs_periodic(Lx, Ly, ALLGS_1)
a_rec = SW.fillEmptyStates_periodic(a, Lx, Ly, ALLGS_1)
##
function SC1(x)
    SW.stencilConfig(x, 1.0, boundaryCondition=:periodic)
end
filter(SW.fulFillsConstraint ∘ SC1, a)