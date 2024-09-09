using Infiltrator
using Revise
using ADRIA
using WGLMakie, GeoMakie, GraphMakie

using DimensionalData, NetCDF, YAXArrays
using ADRIA: AG, GDF

using CSV, DataFrames

using Statistics

if !isdefined(Main, :dom)
    dom = ADRIA.load_domain(RMEDomain, "C:\\Users\\dtan\\repos\\rme_ml_2024_01_08", "45")
end

ltmp_data_dir        = "C:\\Users\\dtan\\data\\ltmp_reefmod_format"

ltmp_hard_soft_fn    = joinpath(ltmp_data_dir, "ltmp_reef_hard_and_soft_coral.gpkg")
ltmp_hard_fn         = joinpath(ltmp_data_dir, "ltmp_reef_hard_cover_ts.gpkg")
manta_tow_fn         = joinpath(ltmp_data_dir, "manta_tow_data_reef_lvl.gpkg")

ltmp_hard_soft       = GDF.read(ltmp_hard_soft_fn)
ltmp_hard            = GDF.read(ltmp_hard_fn)
manta_tow            = GDF.read(manta_tow_fn)

# Sort columns by ascending year
manta_tow_years = parse.(Int64, names(manta_tow)[5:end])
ltmp_hard_years = parse.(Int64, names(ltmp_hard)[7:end])
manta_tow_perm = sortperm(manta_tow_years) .+ 4
ltmp_hard_perm = sortperm(ltmp_hard_years) .+ 6

manta_tow_names = names(manta_tow)
ltmp_hard_names = names(ltmp_hard)
manta_tow_names[5:end] .= manta_tow_names[manta_tow_perm]
ltmp_hard_names[7:end] .= ltmp_hard_names[ltmp_hard_perm]

# Rorder columns
manta_tow = select!(manta_tow, manta_tow_names...)
ltmp_hard = select!(ltmp_hard, ltmp_hard_names...)

# Rescale to be proportions
manta_tow[:, 5:end] ./= 100
ltmp_hard[:, 7:end] ./= 100

ltmp_hard_soft_geoms = ltmp_hard_soft.geometry
ltmp_hard_geoms      = ltmp_hard.geometry
manta_tow_geoms      = manta_tow.geometry

reefmod_historic_dir = "C:\\Users\\dtan\\data\\reefmod_domain_compressed"
reefmod_historic     = open_dataset(joinpath(reefmod_historic_dir, "ReefMod_RCP85.nc"))
reefmod_cover        = reefmod_historic.coral_cover_per_taxa[timestep=1:15]

loc_class = CSV.read("C:\\Users\\dtan\\repos\\ADRIA.jl\\sandbox\\ltmp_calibration\\spatial_data\\location_classification_MPA.csv", DataFrame)

function plot_reefmod_ltmp(geo_row)::Figure
    col_idx::Int64 = findfirst(x -> x[1] == '1', names(geo_row))
    id = geo_row.RME_UNIQUE_ID
    if ismissing(id)
        @warn "Skipping Missing ID"
        return Figure()
    end
    idx = findfirst(x->x==id, dom.site_data.UNIQUE_ID)

    cover = dropdims(sum(reefmod_cover[location=idx], dims=:group), dims=:group)
    uniq_cover = unique(cover[timestep=2])
    idxs = [findfirst(x->cover[timestep=2, scenario=x][1] == val_c, 1:size(cover, 2)) for val_c in uniq_cover]
    cover = cover[scenario=idxs]
    mn = dropdims(mean(cover, dims=:scenario), dims=:scenario)

    missing_mask = (!).(ismissing.(Vector(geo_row[col_idx:end])))
    xs = parse.(Int64, names(geo_row)[col_idx:end])[missing_mask]
    sorted_perm = sortperm(xs)
    xs = xs[sorted_perm]

    f = Figure(; size=(1200, 800))
    Axis(f[1, 1], xlabel="year", ylabel="relative cover", title="LTMP and ReefMod at $(id)")
    obs = scatter!(xs, Vector(geo_row[col_idx:end][missing_mask][sorted_perm]); color=:transparent, strokewidth=2, strokecolor=:black, markersize=15)
    sr = lines!(2008:2022, mn.data[:], color=:red, linewidth=5)
    Legend(
        f[1, 2],
        [obs, sr],
        ["LTMP", "ReefMod"]
    )
    series!(2008:2022, cover.data[:, :]', solid_color=(:red, 0.1))
    if col_idx > 5
        save("reefmod_ltmp/Figures/photo/photo_$(id).png", f)
    else
        save("reefmod_ltmp/Figures/manta_tow/manta_tow_$(id).png", f)
    end
    return f
end

function _get_classification(geo_row)::Int64
    id = geo_row.RME_UNIQUE_ID
    if ismissing(id)
        @warn "Skipping Missing ID"
        return -1
    end
    idx = findfirst(x->x==id, dom.site_data.UNIQUE_ID)

    return loc_class.consecutive_classification[idx]
end

manta_tow_classification = _get_classification.(eachrow(manta_tow))
ltmp_hard_classification = _get_classification.(eachrow(ltmp_hard))

function matrix_flatten_ignore(mat)
    return [el for el in mat if !ismissing(el)]
end

function subdf_mean(flattened)::Float64
    return mean(flattened)
end
function subdf_std(flattened)::Float64
    penalty = 1.0
    n_points = length(flattened)
    stdev = std(flattened)
    if n_points == 1
        return 0.2 # defaults std
    elseif n_points <= 6
        penalty += 5 / n_points
    end
    stdev = stdev < 0.025 ? 0.1 : stdev
    return stdev * penalty
end
function mean_std_series(df)::Tuple{Vector{Union{Missing, Float64}}, Vector{Union{Missing, Float64}}}
    start_index::Int64 = findfirst(x -> x == "2008", names(df[1, :]))
    end_index::Int64 = size(df, 2)

    mean_vec::Vector{Union{Missing, Float64}} = Vector{Union{Missing, Float64}}(
        missing,
        end_index - start_index + 1
    )
    std_vec::Vector{Union{Missing, Float64}} = Vector{Union{Missing, Float64}}(
        missing,
        end_index - start_index + 1
    )
    vec_idx::Int64 = 1
    for yr in start_index:end_index
        time_lb = max(start_index, yr - 1)
        time_ub = min(end_index, yr + 1)
        flattened = matrix_flatten_ignore(Matrix(df[:, time_lb:time_ub]))
        window_penalty = 1 # default to no penalty
        if length(flattened) == 0
            time_lb = max(start_index, yr - 1)
            time_ub = min(end_index, yr + 1)
            flattened = matrix_flatten_ignore(Matrix(df[:, time_lb:time_ub]))
            window_penalty = 3
        end
        if length(flattened) == 0
            time_lb = max(start_index, yr - 2)
            time_ub = min(end_index, yr + 2)
            flattened = matrix_flatten_ignore(Matrix(df[:, time_lb:time_ub]))
            window_penalty = 5
        end
        if length(flattened) == 0
            vec_idx += 1
            continue
        end
        mean_vec[vec_idx] = subdf_mean(flattened)
        std_vec[vec_idx] = subdf_std(flattened) * window_penalty
        vec_idx += 1
    end
    return mean_vec, std_vec
end

function plot_line_band(subdf)::Nothing
    mean_vec, std_vec = mean_std_series(subdf)
    non_missing_mask = (!).(ismissing.(mean_vec))
    xs = 2008:2022
    if length(non_missing_mask) == 15
        xs = (2008:2022)[non_missing_mask]
    else
        xs = (2008:2021)[non_missing_mask]
    end
    if count(non_missing_mask) == 0
        @warn "No years"
        return nothing
    end
    lines!(
        xs,
        mean_vec[non_missing_mask],
        color = :black
    )
    band!(
        xs,
        mean_vec[non_missing_mask] - std_vec[non_missing_mask],
        mean_vec[non_missing_mask] + std_vec[non_missing_mask],
        color = (:black, 0.2)
    )
    return nothing
end
function plot_loc(geo_row)::Nothing
    col_idx::Int64 = findfirst(x -> x[1] == '1', names(geo_row))

    missing_mask = (!).(ismissing.(Vector(geo_row[col_idx:end])))
    xs = parse.(Int64, names(geo_row)[col_idx:end])[missing_mask]
    scatter!(xs, Vector(geo_row[col_idx:end][missing_mask]); color=:transparent, strokewidth=2, strokecolor=:black, markersize=10)

    return nothing
end
function plot_classification(df, classes, class; plot_line=false)::Figure
    @info "Class: $(class)"
    class_mask = classes .== class
    @info "No. Sites: $(count(class_mask))"
    locs = df[class_mask, :]
    f = Figure(; size=(1200, 900))
    Axis(f[1, 1]; xlabel="Year", ylabel="cover", title="Classification: $(class)")
    for l in eachrow(locs)
        plot_loc(l)
    end
    subdir_app = ""
    if plot_line
        subdur_app = "_line"
        plot_line_band(locs)
    end
    col_idx::Int64 = findfirst(x -> x[1] == '1', names(df[1, :]))
    if col_idx > 5
        save("reefmod_ltmp/Figures/loc_classification/photo$(subdir_app)/photo_$(class).png", f)
    else
        save("reefmod_ltmp/Figures/loc_classification/manta_tow$(subdir_app)/manta_tow_$(class).png", f)
    end
    return f
end

manta_tow_classes = unique(manta_tow_classification)
ltmp_hard_classes = unique(ltmp_hard_classification)

manta_tow_yrs = 2008:2022
ltmp_hard_yrs = 2008:2021

manta_tow_mean::Matrix{Union{Missing, Float64}} = Matrix{Union{Missing, Float64}}(
    missing,
    length(manta_tow_classes),
    length(manta_tow_yrs)
)
manta_tow_std::Matrix{Union{Missing, Float64}} = Matrix{Union{Missing, Float64}}(
    missing,
    length(manta_tow_classes),
    length(manta_tow_yrs)
)
ltmp_hard_mean::Matrix{Union{Missing, Float64}} = Matrix{Union{Missing, Float64}}(
    missing,
    length(ltmp_hard_classes),
    length(ltmp_hard_yrs)
)
ltmp_hard_std::Matrix{Union{Missing, Float64}} = Matrix{Union{Missing, Float64}}(
    missing,
    length(ltmp_hard_classes),
    length(ltmp_hard_yrs)
)

manta_tow_axlist = (
    Dim{:class}(manta_tow_classes),
    Dim{:timesteps}(manta_tow_yrs)
)
ltmp_hard_axlist = (
    Dim{:class}(ltmp_hard_classes),
    Dim{:timesteps}(ltmp_hard_yrs)
)

using ProgressMeter

@showprogress for (idx, cls) in enumerate(manta_tow_classes)
    class_mask = manta_tow_classification .== cls
    subdf = manta_tow[class_mask, :]
    mn, stdev = mean_std_series(subdf)
    manta_tow_mean[idx, :] .= mn
    manta_tow_std[idx, :] .= stdev
end

@showprogress for (idx, cls) in enumerate(ltmp_hard_classes)
    class_mask = ltmp_hard_classification .== cls
    subdf = ltmp_hard[class_mask, :]
    mn, stdev = mean_std_series(subdf)
    ltmp_hard_mean[idx, :] .= mn
    ltmp_hard_std[idx, :] .= stdev
end

manta_tow_yax_mean = YAXArray(manta_tow_axlist, manta_tow_mean)
manta_tow_yax_std = YAXArray(manta_tow_axlist, manta_tow_std)

ltmp_hard_yax_mean = YAXArray(ltmp_hard_axlist, ltmp_hard_mean)
ltmp_hard_yax_std = YAXArray(ltmp_hard_axlist, ltmp_hard_std)

manta_tow_dataset = Dataset(; :mean => manta_tow_yax_mean, :std => manta_tow_yax_std)
ltmp_hard_dataset = Dataset(; :mean => ltmp_hard_yax_mean, :std => ltmp_hard_yax_std)

savedataset(manta_tow_dataset, path="reefmod_ltmp/manta_tow_mean_std.nc", backend=:netcdf, overwrite=true)
savedataset(ltmp_hard_dataset, path="reefmod_ltmp/ltmp_hard_mean_std.nc", backend=:netcdf, overwrite=true)
