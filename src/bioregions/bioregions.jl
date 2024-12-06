"""
Visualise the spatial overlap and data availability between AIMS Long Term Monitoring
Program Manta Tow Sites and GBRMPA Bioregions.
"""


using ADRIA
using DataFrames
using ProgressMeter
using Statistics

using CairoMakie, GeoMakie, GraphMakie

import ArchGDAL as AG
import GeoDataFrames as GDF
import GeoInterface as GI

include("common.jl")

dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_DIR, "45")

bioregions_df = GDF.read(
    bioregion_shp_path
)

# Load manta tow ltmp reef level data
ltmp_reef_data = GDF.read(OUT_RME_MANTA)

# Order year columns in ascending order
ltmp_reef_years = parse.(Int64, names(ltmp_reef_data)[5:end])
ltmp_reef_perm = sortperm(ltmp_reef_years) .+ 4

ltmp_reef_data_names = names(ltmp_reef_data)
ltmp_reef_data_names[5:end] .= ltmp_reef_data_names[ltmp_reef_perm]

# Rorder columns
ltmp_reef_data = select!(ltmp_reef_data, ltmp_reef_data_names...)

# Rescale to be proportions
ltmp_reef_data[:, 5:end] ./= 100

ltmp_reef_data = ltmp_reef_data[(!).(ismissing.(ltmp_reef_data.RME_UNIQUE_ID)), :]

first_yr_idx = findfirst(x -> x == "2008", names(ltmp_reef_data))

ltmp_RME_UNIQUE_ID = ltmp_reef_data.RME_UNIQUE_ID
domain_ltmp_idx = [
    findfirst(x->x==id, dom.loc_data.UNIQUE_ID) for id in ltmp_RME_UNIQUE_ID
]

ltmp_reef_data.geometry = GDF.reproject(
    ltmp_reef_data.geometry,
    GI.crs(ltmp_reef_data),
    GI.crs(bioregions_df)
)

function plot_ltmp_locs(xs, ys, title::String, n_locs::Int64; dom=dom)::Figure

    opts::Dict{Symbol, Any} = Dict(
        :show_colorbar => false
    )
    axis_opts::Dict{Symbol, Any} = Dict(
        :title => title * "\n$(n_locs) Locations"
    )

    f = ADRIA.viz.map(dom; opts=opts, axis_opts=axis_opts)
    scatter!(xs, ys, color=(:blue, 0.4), overdraw=true)
    save(joinpath(OUTPUT_DIR, "ltmp_spatial", "$(title).png"), f)

    return f
end

mkpath(joinpath(OUTPUT_DIR, "ltmp_spatial"))

@showprogress desc="Plotting LTMP Data Availability" for i in 1:14
    obs_count_mask = [
        count((!).(ismissing.(collect(df_row[first_yr_idx:end])))) > i for df_row in eachrow(ltmp_reef_data)
    ]

    n_locs = count(obs_count_mask)

    local xs = dom.loc_data.X_COORD[domain_ltmp_idx][obs_count_mask]
    local ys = dom.loc_data.Y_COORD[domain_ltmp_idx][obs_count_mask]

    local f = plot_ltmp_locs(xs, ys, "More than $(i) observations", n_locs)
end
@info "Saved LTMP Data availability plots to $(joinpath(OUTPUT_DIR, "ltmp_spatial")) directory."

xs = dom.loc_data.X_COORD[domain_ltmp_idx]
ys = dom.loc_data.Y_COORD[domain_ltmp_idx]

f = plot_ltmp_locs(xs, ys, "All LTMP Locs", length(xs))

avail_masks = [
    [
        count((!).(ismissing.(collect(df_row[first_yr_idx:end])))) > i for df_row in eachrow(ltmp_reef_data)
    ] for i in 0:14
]

function get_ltmp_idx(ltmp_pt, bioregs)
    possible_idx = findfirst(x -> AG.contains(x, ltmp_pt), bioregs.geometry)
    if isnothing(possible_idx)
        dists = [AG.distance(ltmp_pt, bio_geom) for bio_geom in bioregs.geometry]
        possible_idx = argmin(dists)
        if dists[possible_idx] > 1.0
            return nothing
        end
    end
    return possible_idx
end

bioregion_ltmp_idx = [
    get_ltmp_idx(pt, bioregions_df) for pt in ltmp_reef_data.geometry
]

ltmp_bioregion = [
    isnothing(idx) ? -1 : bioregions_df[idx, :BIOREGION] for idx in bioregion_ltmp_idx
]

mkpath(joinpath(OUTPUT_DIR, "bioregions"))

@showprogress desc="Plotting LTMP Bioregion Alignment" for (idx, mask) in zip(0:14, avail_masks)

    ltmp_bioregion_atleast_idx_obs = ltmp_bioregion[ltmp_bioregion .!= -1 .&& mask]
    local f = Figure()
    Axis(
        f[1, 1],
        xlabel="Bioregion",
        ylabel="Number of Reefs",
        title="Number of LTMP reefs in each Bioregion with at least $(idx+1) observations",
        limits=(0, 40, 0, 45)
    )
    hist!(ltmp_bioregion_atleast_idx_obs, bins = 1:maximum(bioregions_df.BIOREGION))
    save(joinpath(OUTPUT_DIR, "bioregions", "ltmp_hist_atleast_$(idx+1).png"), f)
end
@info "Saved Bioregion LTMP alignment plots to $(joinpath(OUTPUT_DIR, "bioregions")) directory."

domain_bioregion_idxs = [
    findfirst(
        x->AG.intersects(dom_geom, x),
        bioregions_df.geometry
    ) for dom_geom in dom.loc_data.geom
]

for (dom_idx, bio_idx) in enumerate(domain_bioregion_idxs)
    if !isnothing(bio_idx)
        continue
    end

    dists = [
        AG.distance(dom.loc_data.geom[dom_idx], bio_geom)
        for bio_geom in bioregions_df.geometry
    ]
    domain_bioregion_idxs[dom_idx] = argmin(dists)
end

domain_bioregions = [bioregions_df.BIOREGION[j] for j in domain_bioregion_idxs]

opts::Dict{Symbol, Any} = Dict(
    :colorbar_label => "Bioregion",
    :color_map => :viridis
)
axis_opts::Dict{Symbol, Any} = Dict(
    :title => "GBRMPA Bioregions"
)

f = ADRIA.viz.map(dom, domain_bioregions; opts=opts, axis_opts=axis_opts)
save(joinpath(OUTPUT_DIR, "bioregions_map.png"), f)

function get_temporal_range(dfr::DataFrameRow)
    first_yr = findfirst(x -> !ismissing(x), dfr)
    last_yr  = findlast(x -> !ismissing(x), dfr)

    if isnothing(first_yr) || isnothing(last_yr)
        return 0
    end

    return parse(Int64, String(last_yr)) - parse(Int64, String(first_yr))
end

function count_data_points(dfr::DataFrameRow)::Int64
    return count((!).(ismissing.(collect(dfr))))
end

mkpath(joinpath(OUTPUT_DIR, "bioregion_dist"))

"""
    plot_bioregion_locs(loc_mask, title::String, n_locs::Int64; dom=dom)::Figure
"""
function plot_bioregion_locs(loc_mask, title::String; dom=dom)::Figure

    opts::Dict{Symbol, Any} = Dict(
        :show_colorbar => false
    )
    axis_opts::Dict{Symbol, Any} = Dict(
        :title => title
    )

    f = ADRIA.viz.map(dom, 1 .- loc_mask; opts=opts, axis_opts=axis_opts)
    save(joinpath(OUTPUT_DIR, "bioregion_dist", "$(title).png"), f)

    return f
end

# Plot the distribution of location for each bioregion
@showprogress desc="Plotting Bioregions" for bioreg in unique(domain_bioregions)
    target_bioreg_mask = domain_bioregions .== bioreg
    plot_bioregion_locs(target_bioreg_mask, "bioregion_$(bioreg)")
end
@info "Saved Bioregions Location distribution plots to $(joinpath(OUTPUT_DIR, "bioregion_dist")) directory."

function enough_data(df::DataFrame)
    return [count_data_points(dfr) >= 4 && get_temporal_range(dfr) >= 10 for dfr in eachrow(df)]
end

locs_enough_data = enough_data(ltmp_reef_data[:, first_yr_idx:end])

all_bioregs = unique(bioregions_df.BIOREGION)

ltmp_counts = zeros(Int64, length(all_bioregs))

for (idx, bioreg) in enumerate(all_bioregs)
    bioreg_mask =  ltmp_bioregion .== bioreg
    msk = bioreg_mask .&& locs_enough_data
    ltmp_counts[idx] = count(msk)
end

enough_data_mask = fill(false, 3806)
for (bioreg, cnt) in zip(all_bioregs, ltmp_counts)
    if cnt >= 4
        enough_data_mask .|= domain_bioregions .== bioreg
    end
end

plot_bioregion_locs(enough_data_mask, "enough_data")
plot_bioregion_locs(1 .- enough_data_mask, "not_enough_data")

function construct_distance_matrix(geoms1, geoms2)::Matrix{Float64}
    n1::Int64 = length(geoms1)
    n2::Int64 = length(geoms2)

    dists::Matrix{Float64} = Matrix{Float64}(undef, n1, n2)
    for i in 1:n1, j in 1:n2
        dists[i, j] = AG.distance(geoms1[i], geoms2[j])
    end

    return dists
end

"""
Quantify the distance between two bioregions by taking the 10% quantile of all distances
between locations. The distances are only meant to be used in a relative manner.
"""
function bioregion_distance(
    b1::Int64,
    b2::Int64;
    dom_bioregs=domain_bioregions,
    domain=dom
)::Float64
    bio1::DataFrame = domain.loc_data[dom_bioregs .== b1, :]
    bio2::DataFrame = domain.loc_data[dom_bioregs .== b2, :]

    dists_mats = construct_distance_matrix(bio1.geom, bio2.geom)

    return quantile(vec(dists_mats), 0.1)
end


bioregion_no_data = ltmp_counts .< 4

# Grouped bioregions guarenteeing enough data
bioregion_assignment = zeros(Int64, length(all_bioregs))

# bioregions with sufficient ltmp data are assignmed to themselves
bioregion_assignment[(!).(bioregion_no_data)] = all_bioregs[(!).(bioregion_no_data)]

@showprogress desc="Assigning Bioregions" for (idx, bioreg) in enumerate(all_bioregs)
    if bioregion_assignment[idx] != 0
        continue
    end
    # find the closest bioregion with enough data
    min_dist = 10000
    for bioreg_w_data in all_bioregs[(!).(bioregion_no_data)]
        bioreg_dist = bioregion_distance(bioreg, bioreg_w_data)
        if bioreg_dist < min_dist
            bioregion_assignment[idx] = bioreg_w_data
            min_dist = bioreg_dist
        end
    end
end

"""
    orginal_bio_to_assignmed_bio(original_bio::Int64; original_bio_regs=all_bioregs, assigned_bio_regs=bioregion_assignment)::Int64

Convert a bioregion to the reassigned bioregion.
"""
function orginal_bio_to_assignmed_bio(
    original_bio::Int64;
    original_bio_regs=all_bioregs,
    assigned_bio_regs=bioregion_assignment
)::Int64
    idx::Int64 = findfirst(x->x==original_bio, original_bio_regs)
    if isnothing(idx)
        throw(ArgumentError("Bioregion $(idx) not found in original bioregions"))
    end
    return assigned_bio_regs[idx]
end

domain_assigned_bios = orginal_bio_to_assignmed_bio.(domain_bioregions)

@showprogress desc="Plotting Assignmed Bioregions" for ass_bio in unique(bioregion_assignment)
    target_bioreg_mask = domain_assigned_bios .== ass_bio
    plot_bioregion_locs(target_bioreg_mask, "assigned_bioregion_$(ass_bio)")
end

ltmp_assigned_bioregions = orginal_bio_to_assignmed_bio.(ltmp_bioregion)
ltmp_assigned_bioregions = ltmp_assigned_bioregions[locs_enough_data]


f = Figure()
Axis(
    f[1, 1],
    xlabel="Assignmed Bioregion",
    ylabel="Number of Reefs",
    title="Number of LTMP reefs in each Assignmed Bioregion",
    limits=(0, 40, 0, 45)
)
hist!(ltmp_assigned_bioregions, bins = 1:maximum(bioregions_df.BIOREGION))
save(joinpath(OUTPUT_DIR, "final_ltmp_assigned_bio_counts.png"), f)

domain_gpkg_copy = copy(dom.loc_data)
domain_gpkg_copy[!, :ASSIGNED_BIOREGION] .= domain_assigned_bios

@info "Writing gpkg with bioregions to path: $(OUTPUT_DOMAIN_BIOREGION_GPKG)."
GDF.write(OUTPUT_DOMAIN_BIOREGION_GPKG, domain_gpkg_copy; geom_columns=(:geom,))
