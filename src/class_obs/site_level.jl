"""
Get LTMP observation for location classifications. Calculate the standard deviation as well.
"""

include("common.jl")

using ADRIA
using ADRIA: AG, GDF

using CSV,
    DataFrames,
    DimensionalData,
    NetCDF,
    Statistics,
    YAXArrays

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

if !isdefined(Main, :dom)
    dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_DIR, "45")
end

ltmp_hard            = GDF.read(OUT_RME_PHOTO_HC)
manta_tow            = GDF.read(OUT_RME_MANTA)

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

loc_class = CSV.read(OUTPUT_CSV, DataFrame)
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

for (idx, cls) in enumerate(manta_tow_classes)
    class_mask = manta_tow_classification .== cls
    subdf = manta_tow[class_mask, :]
    mn, stdev = mean_std_series(subdf)
    manta_tow_mean[idx, :] .= mn
    manta_tow_std[idx, :] .= stdev
end

for (idx, cls) in enumerate(ltmp_hard_classes)
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

@info "Writing LTMP Manta Tow Observations for Location Classes to $(OUT_CLASS_MANTA)"
savedataset(manta_tow_dataset, path=OUT_CLASS_MANTA, backend=:netcdf, overwrite=true)
@info "Writing LTMP PHOTO TS Observations for Location Classes to $(OUT_CLASS_MANTA)"
savedataset(ltmp_hard_dataset, path=OUT_CLASS_PHOTO_HC, backend=:netcdf, overwrite=true)
