include("utils.jl")
include("common.jl")
using CSV,
    Dates,
    DataFrames,
    DimensionalData,
    YAXArrays,
    NetCDF


ltmp_site_data = CSV.read(LTMP_PHOTO_FN, DataFrame)
manta_tow_data = CSV.read(LTMP_MANTA_FN, DataFrame)

# Add years
ltmp_site_data[!, :REPORT_YEAR] = Dates.year.(ltmp_site_data.SAMPLE_DATE)

hard_coral_mask = ltmp_site_data.GROUP_CODE .== "Hard Coral"
soft_coral_mask = ltmp_site_data.GROUP_CODE .== "Soft Coral"
algae_mask      = ltmp_site_data.GROUP_CODE .== "Algae"
other_mask      = ltmp_site_data.GROUP_CODE .== "Other"

ltmp_hard_coral = ltmp_site_data[hard_coral_mask, :]
ltmp_soft_coral = ltmp_site_data[soft_coral_mask, :]
ltmp_algae      = ltmp_site_data[algae_mask, :]
ltmp_other      = ltmp_site_data[other_mask, :]

manta_tow_data[!, :MEDIAN_LIVE_CORAL_LB] = lower_bound.(manta_tow_data.MEDIAN_LIVE_CORAL)
manta_tow_data[!, :MEDIAN_SOFT_CORAL_LB] = lower_bound.(manta_tow_data.MEDIAN_SOFT_CORAL)
manta_tow_data[!, :MEDIAN_DEAD_CORAL_LB] = lower_bound.(manta_tow_data.MEDIAN_DEAD_CORAL)
manta_tow_data[!, :MEDIAN_LIVE_CORAL_UB] = upper_bound.(manta_tow_data.MEDIAN_LIVE_CORAL)
manta_tow_data[!, :MEDIAN_SOFT_CORAL_UB] = upper_bound.(manta_tow_data.MEDIAN_SOFT_CORAL)
manta_tow_data[!, :MEDIAN_DEAD_CORAL_UB] = upper_bound.(manta_tow_data.MEDIAN_DEAD_CORAL)

loc_ids = unique(manta_tow_data.REEF_ID)
years   = sort(unique(manta_tow_data.REPORT_YEAR))

vars = [
    :MEDIAN_LIVE_CORAL_LB,
    :MEDIAN_SOFT_CORAL_LB,
    :MEDIAN_DEAD_CORAL_LB,
    :MEDIAN_LIVE_CORAL_UB,
    :MEDIAN_SOFT_CORAL_UB,
    :MEDIAN_DEAD_CORAL_UB,
    :MEAN_LIVE_CORAL,
    :MEAN_SOFT_CORAL,
    :MEAN_DEAD_CORAL,
]

function create_var(df, var_name; T=Float64)
    tsteps = sort(unique(df.REPORT_YEAR))
    location_ids = unique(df.REEF_ID)

    axlist = (
        Dim{:timesteps}(tsteps),
        Dim{:locs}(String.(location_ids))
    )

    preallocation = Matrix{Union{Missing, T}}(
        missing,
        length(tsteps),
        length(location_ids)
    )

    for (idx, loc_id) in enumerate(location_ids)
        sub_df   = df[df.REEF_ID .== loc_id, :]
        idx_mask = [findfirst(x -> x == l, tsteps) for l in sub_df.REPORT_YEAR]
        preallocation[idx_mask, idx] .= sub_df[:, var_name]
    end

    return var_name => YAXArray(axlist, preallocation)
end


var_creater = x -> create_var(manta_tow_data, x)

manta_tow = Dataset(; var_creater.(vars)...)
@info "Writing manta tow observations for reefmod location to $(OUT_MANTA_NETCDF)"
savedataset(manta_tow; path=OUT_MANTA_NETCDF, driver=:netcdf, overwrite=true)

hard_coral = create_var(ltmp_hard_coral, :COVER)
soft_coral = create_var(ltmp_soft_coral, :COVER)
algae = create_var(ltmp_algae, :COVER)
other = create_var(ltmp_other, :COVER)

ltmp_hc = Dataset(
    ;
    :HARD_CORAL => hard_coral[2],
    :SOFT_CORAL => soft_coral[2],
    :ALGAE => algae[2],
    :OTHER => other[2]
)

@info "Writing manta tow observations for reefmod location to $(OUT_MANTA_NETCDF)"
savedataset(ltmp_hc; path=OUT_PHOTO_NETCDF, driver=:netcdf, overwrite=true)
