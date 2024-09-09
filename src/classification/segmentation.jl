using ADRIA
using WGLMakie, GeoMakie, GraphMakie

using CSV, DataFrames

using Statistics

features_fn = "raster/shelf_interpolated_results_MPA.csv"
features = CSV.read(features_fn, DataFrame)

function plot_var(df, var_name)
    shelf_names = ["Inner", "Mid", "Outer"]
    region_names = ["North", "Central", "South"]
    f = Figure(; size=(900, 900))
    for shelf_pos in 1:3
        axs = []
        for region in 1:3
            loc_mask = df.shelf_position .== shelf_pos .&& df.ltmp_region .== region
            ax = Axis(f[region, shelf_pos]; xlabel=String(var_name), ylabel="location count", title="$(region_names[region]), $(shelf_names[shelf_pos]), $(String(var_name))")
            hist!(df[loc_mask, var_name]; strokewidth=1, strokecolor=:black)
            push!(axs, ax)
        end
        linkxaxes!(axs...)
    end
    return f
end


"""
The classification of locations is done using base 3.
- 3^0 : bathy
- 3^1 : turbid
- 3^2 : ubed
- 3^3 : shelf position
- 3^4 : ltmp region

add one after converting to base 3.
"""

bathy_bounds = quantile(features.bathy_mean, [0.66, 0.33])
ub_bounds = quantile(features.ub_mean, [0.33, 0.66])
turbid_bounds = quantile(features.turbid_mean, [0.33, 0.66])


function interval_idx(elmt::Float64, intervals::Vector{Float64}; increasing::Bool=true)::Int64
    if increasing
        for (idx, val) in enumerate(intervals)
            elmt <= val && return idx - 1
        end
        return length(intervals)
    else
        for (idx, val) in enumerate(intervals)
            elmt >= val && return idx - 1
        end
        return length(intervals)
    end
end

function classify_location(df_row)::Int64
    classification::Int64 = 0
    classification += interval_idx(df_row.bathy_mean, bathy_bounds; increasing=false)
    classification += interval_idx(df_row.turbid_mean, turbid_bounds; increasing=false) * 3
    classification += interval_idx(df_row.ub_mean, ub_bounds; increasing=false) * 3^2
    classification += (df_row.shelf_position - 1) * 3^3
    classification += (df_row.ltmp_region - 1) * 3^4

    return classification + 1
end

features[!, :classification] .= classify_location.(eachrow(features))
features[!, :consecutive_classification] .= [findfirst(x -> x == class, unique(features.classification)) for class in features.classification]
