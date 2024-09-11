include("common.jl")

using ADRIA
using CSV, DataFrames

import ArchGDAL as AG

if !isdefined(Main, :dom)
    dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_DIR, "45")
end

res = CSV.read(OUTPUT_CSV, DataFrame)

function norm(x::Float64, y::Float64)::Float64
    # Prioritise shelf position
    return sqrt(x^2 + y^2)
end

"""
    interpolate!(df::DataFrame, missing_idx::Int64, non_missing::DataFrame)::Nothing

Fill in missing values as a weighted average of closest locations in the same shelf position.
Weightings are calculated as the normalisation of 10 closest distances.
"""
function interpolate!(df::DataFrame, missing_idx::Int64, non_missing::DataFrame)::Nothing

    shelf_position = df.shelf_position[missing_idx]
    shelf_pos_mask = non_missing.shelf_position .== shelf_position
    masked_df = non_missing[shelf_pos_mask, :]

    distances::Vector{Float64} = norm.(
        masked_df.lons .- df.lons[missing_idx], masked_df.lats .- df.lats[missing_idx]
    )
    n_neighbours = 10
    closest = sortperm(distances)[1:n_neighbours]
    weightings = distances[closest] ./ sum(distances[closest])
    df[missing_idx, :] .= sum.(eachcol(masked_df[closest, :] .* weightings))

    return nothing
end

# Locations with/without remote sensing data
missing_data_mask = res.slp_prop .== 0.0 .&& res.flt_prop .== 0.0
missing_idxs = (1:3806)[missing_data_mask]

non_missing_data_mask = (!).(missing_data_mask)

centroids = AG.centroid.(dom.site_data.geom)

lons = AG.getx.(centroids, 0)
lats = AG.gety.(centroids, 0)

res[!, :lons] .= lons
res[!, :lats] .= lats

non_missing_df = res[non_missing_data_mask, :]

@info "Filling Missing Values"
for idx in missing_idxs
    interpolate!(res, idx, non_missing_df)
end

# Reassign lat and lons
res[!, :lons] .= lons
res[!, :lats] .= lats

@info "Writing updated location data to $(OUTPUT_CSV)"
CSV.write(OUTPUT_CSV, res)
