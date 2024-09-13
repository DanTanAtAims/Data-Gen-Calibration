"""
Calculate the mean and standard deviation for wave statistics (ub, hs, tp), bathymetry,
and turbidity.
"""


include("common.jl")

using ADRIA, Rasters

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using CSV, DataFrames, Statistics

"""
    clean_and_combine_vectors(clean, unclean)

Remove Nan and Missing Values from unclean vector and add to clean vector.
"""
function clean_and_combine_vectors(clean, unclean)
    unclean[ismissing.(unclean)] .= 0
    unclean[isnan.(unclean)] .= 0
    clean += unclean
    return clean
end

"""
    calc_stats(slps, flts, habitable, var_rast, dom)

Calculate the mean and standard deviation of a given variable at each polygon in the domain.
Furthermore, calculate the mean and standard deviation for the slope and flats.
Geomorphologies are masked by habitability.
"""
function calc_stats(slps, flts, habitable, var_rast, dom)
    var_habitable  = Rasters.mask(var_rast; with=habitable, boundary=:touches, progress=false, verbose=false)
    var_slopes     = Rasters.mask(var_habitable; with=slps, boundary=:touches, progress=false, verbose=false)
    var_flats      = Rasters.mask(var_habitable; with=flts, boundary=:touches, progress=false, verbose=false)

    var_mean       = Rasters.zonal(mean, var_habitable; of=dom.site_data.geom, progress=false, verbose=false)
    var_std        = Rasters.zonal(std,  var_habitable; of=dom.site_data.geom, progress=false, verbose=false)

    var_slope_mean = Rasters.zonal(mean, var_slopes; of=dom.site_data.geom, progress=false, verbose=false)
    var_slope_std  = Rasters.zonal(std,  var_slopes; of=dom.site_data.geom, progress=false, verbose=false)

    var_flats_mean = Rasters.zonal(mean, var_flats; of=dom.site_data.geom, progress=false, verbose=false)
    var_flats_std  = Rasters.zonal(std,  var_flats; of=dom.site_data.geom, progress=false, verbose=false)

    return var_flats_mean, var_flats_std, var_slope_mean, var_slope_std, var_mean, var_std
end

"""
    total_aggregation(x, y)

Add two vectors togethor, setting missing values to 0.
"""
function total_aggregation(x, y)
    if ismissing(x) && ismissing(y)
        return 1.0
    elseif ismissing(x)
        return y
    elseif ismissing(y)
        return x
    else
        return x + y
    end
end

"""
    slope_flat_prop(slps, flts, habitable, one_rast, dom)

Calculate the proportion of polygon covered by habitable slope and flats.
"""
function slope_flat_prop(slps, flts, habitable, one_rast, dom)
    one_habitable = Rasters.mask(one_rast; with=habitable, boundary=:touches, progress=false, verbose=false)
    one_slopes = Rasters.mask(one_habitable; with=slps, bounary=:touches, progress=false, verbose=false)
    one_flats = Rasters.mask(one_habitable; with=flts, boundary=:touches, progress=false, verbose=false)

    slp_count = Rasters.zonal(sum, one_slopes; of=dom.site_data.geom, progress=false, verbose=false)
    flt_count = Rasters.zonal(sum, one_flats; of=dom.site_data.geom, progress=false, verbose=false)

    total = total_aggregation.(slp_count, flt_count)
    return flt_count ./ total, slp_count ./ total
end

dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_DIR, "45")

var_names = [
    "ub_", # wave
    "hs_", # wave
    "tp_", # wave
    "bathy_",
    "turbid_"
]

geomorphic_names = [
    "",
    "slp_", # slope
    "flt_" # flat
]

stat_names = [
    "mean",
    "std"
]

df_columns = []
for v in var_names, g in geomorphic_names, s in stat_names
    push!(df_columns, Symbol(g * v * s))
end
push!(df_columns, :slp_prop)
push!(df_columns, :flt_prop)

res = DataFrame(zeros(Float64, 3806, 32), df_columns)
# wave UB
for (reg_idx, reg) in enumerate(REGIONS)
    @info "Processing $(reg)"
    slopes_fn  = SLOPE_FNS[reg_idx]
    flats_fn   = FLAT_FNS[reg_idx]
    benthic_fn = BENTHIC_FNS[reg_idx]

    slopes  = GDF.read(slopes_fn)
    flats   = GDF.read(flats_fn)
    benthic = GDF.read(benthic_fn)

    habitable_mask = benthic[benthic.class .!= "Sand", :]

    waves_ub_fn = WAVES_UB_FNS[reg_idx]
    waves_hs_fn = WAVES_HS_FNS[reg_idx]
    waves_tp_fn = WAVES_TP_FNS[reg_idx]
    bathy_fn    = BATHY_FNS[reg_idx]
    turbid_fn   = TURBID_FNS[reg_idx]

    waves_ub = Raster(waves_ub_fn; lazy=true)
    raster_ones = copy(waves_ub)
    raster_ones .= 1

    @info "Processing $(reg): proportions"
    f_mean, s_mean = slope_flat_prop(slopes, flats, habitable_mask, raster_ones, dom)
    res.slp_prop .= clean_and_combine_vectors(res.slp_prop, s_mean)
    res.flt_prop .= clean_and_combine_vectors(res.flt_prop, f_mean)
    raster_ones = nothing
    GC.gc()

    @info "Processing $(reg): waves_ub"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_ub, dom)
    res.ub_mean     .= clean_and_combine_vectors(res.ub_mean, mn)
    res.ub_std      .= clean_and_combine_vectors(res.ub_std,  st)
    res.slp_ub_mean .= clean_and_combine_vectors(res.slp_ub_mean, s_mean)
    res.slp_ub_std  .= clean_and_combine_vectors(res.slp_ub_std,  s_std)
    res.flt_ub_mean .= clean_and_combine_vectors(res.flt_ub_mean, f_mean)
    res.flt_ub_std  .= clean_and_combine_vectors(res.flt_ub_std,  f_std)
    waves_ub = nothing
    GC.gc()

    @info "Processing $(reg): waves_hs"
    waves_hs = Raster(waves_hs_fn; lazy=true)
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_hs, dom)
    res.hs_mean     .= clean_and_combine_vectors(res.hs_mean, mn)
    res.hs_std      .= clean_and_combine_vectors(res.hs_std,  st)
    res.slp_hs_mean .= clean_and_combine_vectors(res.slp_hs_mean, s_mean)
    res.slp_hs_std  .= clean_and_combine_vectors(res.slp_hs_std,  s_std)
    res.flt_hs_mean .= clean_and_combine_vectors(res.flt_hs_mean, f_mean)
    res.flt_hs_std  .= clean_and_combine_vectors(res.flt_hs_std,  f_std)
    waves_hs = nothing
    GC.gc()

    @info "Processing $(reg): waves_tp"
    waves_tp = Raster(waves_tp_fn; lazy=true)
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_tp, dom)
    res.tp_mean .= clean_and_combine_vectors(res.tp_mean, mn)
    res.tp_std  .= clean_and_combine_vectors(res.tp_std,  st)
    res.slp_tp_mean .= clean_and_combine_vectors(res.slp_tp_mean, s_mean)
    res.slp_tp_std  .= clean_and_combine_vectors(res.slp_tp_std,  s_std)
    res.flt_tp_mean .= clean_and_combine_vectors(res.flt_tp_mean, f_mean)
    res.flt_tp_std  .= clean_and_combine_vectors(res.flt_tp_std,  f_std)
    waves_tp = nothing
    GC.gc()

    @info "Processing $(reg): bathy"
    bathy    = Raster(bathy_fn;    lazy=true)
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, bathy, dom)
    res.bathy_mean .= clean_and_combine_vectors(res.bathy_mean, mn)
    res.bathy_std  .= clean_and_combine_vectors(res.bathy_std,  st)
    res.slp_bathy_mean .= clean_and_combine_vectors(res.slp_bathy_mean, s_mean)
    res.slp_bathy_std  .= clean_and_combine_vectors(res.slp_bathy_std,  s_std)
    res.flt_bathy_mean .= clean_and_combine_vectors(res.flt_bathy_mean, f_mean)
    res.flt_bathy_std  .= clean_and_combine_vectors(res.flt_bathy_std,  f_std)
    bathy = nothing
    GC.gc()

    @info "Processing $(reg): turbid"
    turbid   = Raster(turbid_fn;   lazy=true)
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, turbid, dom)
    res.turbid_mean .= clean_and_combine_vectors(res.turbid_mean, mn)
    res.turbid_std  .= clean_and_combine_vectors(res.turbid_std,  st)
    res.slp_turbid_mean .= clean_and_combine_vectors(res.slp_turbid_mean, s_mean)
    res.slp_turbid_std  .= clean_and_combine_vectors(res.slp_turbid_std,  s_std)
    res.flt_turbid_mean .= clean_and_combine_vectors(res.flt_turbid_mean, f_mean)
    res.flt_turbid_std  .= clean_and_combine_vectors(res.flt_turbid_std,  f_std)
    turbid = nothing
    GC.gc()
end

CSV.write(OUTPUT_CSV, res)
