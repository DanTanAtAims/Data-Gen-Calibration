using Infiltrator
using ProgressMeter

using ADRIA
using WGLMakie, GeoMakie, GraphMakie

import ArchGDAL as AG

using CSV, DataFrames

if !isdefined(Main, :dom)
    dom = ADRIA.load_domain(RMEDomain, "C:\\Users\\dtan\\repos\\rme_ml_2024_01_08", "45")
end

df = CSV.read("raster/shelf_result_habitable_MPA.csv", DataFrame)

function weighted_euclidean_norm(x::Float64, y::Float64)::Float64
    # Prioritise shelf position
    return sqrt(x^2 + y^2)
end

function norm(x::Float64, y::Float64)::Float64
    return sqrt(x^2 + y^2)
end

# Locations with/without remote sensing data
missing_data_mask = df.slp_prop .== 0.0 .&& df.flt_prop .== 0.0
missing_idxs = (1:3806)[missing_data_mask]

non_missing_data_mask = (!).(missing_data_mask)

centroids = AG.centroid.(dom.site_data.geom)

lons = AG.getx.(centroids, 0)
lats = AG.gety.(centroids, 0)

df[!, :lons] .= lons
df[!, :lats] .= lats

non_missing_df = df[non_missing_data_mask, :]

function interpolate!(df::DataFrame, missing_idx::Int64, non_missing::DataFrame)::Nothing

    shelf_position = df.shelf_position[missing_idx]
    shelf_pos_mask = non_missing.shelf_position .== shelf_position
    masked_df = non_missing[shelf_pos_mask, :]

    distances::Vector{Float64} = weighted_euclidean_norm.(
        masked_df.lons .- df.lons[missing_idx], masked_df.lats .- df.lats[missing_idx]
    )
    n_neighbours = 10
    closest = sortperm(distances)[1:n_neighbours]
    weightings = distances[closest] ./ sum(distances[closest])
    df[missing_idx, :] .= sum.(eachcol(masked_df[closest, :] .* weightings))

    return nothing
end

@showprogress for idx in missing_idxs
    interpolate!(df, idx, non_missing_df)
end

df[!, :lons] .= lons
df[!, :lats] .= lats

using GeoDataFrames
import GeoDataFrames as GDF
if !isdefined(Main, :regions_GDA2020)
    regions_GDA2020 = GDF.read("C:\\Users\\dtan\\data\\GDA2020-Data-for-PDP\\Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg")
    rename!(regions_GDA2020, Dict(:SHAPE => :geometry))
    GDA2020_crs = crs(regions_GDA2020[1,:geometry])
    ports = GeoDataFrames.read("C:\\Users\\dtan\\data\\QLD_ports_mercator_via_MP\\ports_QLD_merc.shp")
    ports.geometry = AG.reproject(ports.geometry, crs(ports[1, :geometry]), GDA2020_crs; order=:trad)
    portxs = AG.getx.(ports.geometry, 0)
    portys = AG.gety.(ports.geometry, 0)
end

function closest_port(x::Float64, y::Float64)::Float64

    distances::Vector{Float64} = norm.(
        portxs .- x, portys .- y
    )

    return minimum(distances)
end

df[!, :closest_port] .= closest_port.(df.lons, df.lats)

border_shp_fn = "C:\\Users\\dtan\\data\\STE_2021_AUST_SHP_GDA2020\\STE_2021_AUST_GDA2020.shp"
border_shp    = GDF.read(border_shp_fn)
border = border_shp.geometry[3]

dists = Vector{Float64}([AG.distance(border, cent) for cent in centroids])

function shelf_position(centroids, dists, m, n, q1, q2, memory)
    ys = AG.gety.(centroids, 0)
    min_y = minimum(ys)
    max_y = maximum(ys)

    width::Float64 = (max_y - min_y) / m
    lbs = min_y:(max_y - min_y) / n:max_y - width
    ubs = min_y + width:(max_y - min_y) / n:max_y

    counts::Vector{Int64} = zeros(Int64, 3806)
    inner_counts::Vector{Int64} = zeros(Int64, 3806)
    outer_counts::Vector{Int64} = zeros(Int64, 3806)
    mask::BitVector = BitVector([true for i in 1:3806])
    lower_mask::BitVector = BitVector([true for i in 1:3806])
    upper_mask::BitVector = BitVector([true for i in 1:3806])
    upper_bound::Float64 = 0.0
    for (lb, ub) in reverse(collect(zip(lbs, ubs)))
        mask    .= ys .>= lb .&& ys .<= ub
        if !any(mask)
            continue
        end
        quants = quantile(dists[mask], [q1, q2])
        upper_bound = upper_bound == 0.0 ? q2 : upper_bound * memory + quants[2] * (1 - memory)
        counts[mask] .+= 1
        lower_mask .= dists .> quants[1] .&& mask
        upper_mask .= dists .> memory .&& mask
        inner_counts[lower_mask] .+= 1
        outer_counts[upper_mask] .+= 1
    end

    counts[counts .== 0] .= 1

    return inner_counts ./ counts, outer_counts ./ counts
end
