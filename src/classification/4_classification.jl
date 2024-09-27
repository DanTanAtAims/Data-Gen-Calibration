"""
Segment locations into classes and add classification to the output csv.
"""

include("common.jl")

using ADRIA

using CSV, DataFrames, Statistics

features = CSV.read(OUTPUT_CSV, DataFrame)

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

@info "Classifying Locations"
features[!, :classification] .= classify_location.(eachrow(features))
# Some classes contain no locations. New column gives consecutive class ids consecutive.
features[!, :consecutive_classification] .= [findfirst(x -> x == class, unique(features.classification)) for class in features.classification]

features[!, :RME_UNIQUE_ID] .= dom.site_data.UNIQUE_ID

@info "Writing location classification to $(OUTPUT_CSV)"
CSV.write(OUTPUT_CSV, features)
