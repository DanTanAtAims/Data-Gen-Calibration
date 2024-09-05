include("../common.jl")

using ADRIA, Rasters

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using Statistics

dom = ADRIA.load_domain(RMEDomain, RME_DOMAIN_DIR, "45")

function calc_stats(slps, flts, habitable, var_rast, dom)
    var_habitable  = Rasters.mask(var_rast; with=habitable, boundary=:touches, progress=false)
    var_slopes     = Rasters.mask(var_habitable; with=slps, boundary=:touches, progress=false)
    var_flats      = Rasters.mask(var_habitable; with=flts, boundary=:touches, progress=false)

    var_mean       = Rasters.zonal(mean, var_habitable; of=dom.site_data.geom, progress=false)
    var_std        = Rasters.zonal(std,  var_habitable; of=dom.site_data.geom, progress=false)

    var_slope_mean = Rasters.zonal(mean, var_slopes; of=dom.site_data.geom, progress=false)
    var_slope_std  = Rasters.zonal(std,  var_slopes; of=dom.site_data.geom, progress=false)

    var_flats_mean = Rasters.zonal(mean, var_flats; of=dom.site_data.geom, progress=false)
    var_flats_std  = Rasters.zonal(std,  var_flats; of=dom.site_data.geom, progress=false)

    return var_flats_mean, var_flats_std, var_slope_mean, var_slope_std, var_mean, var_std
end

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

function slope_flat_prop(slps, flts, habitable, one_rast, dom)
    one_habitable = Rasters.mask(one_rast; with=habitable, boundary=:touches, progress=false)
    one_slopes = Rasters.mask(one_habitable; with=slps, bounary=:touches, progress=false)
    one_flats = Rasters.mask(one_habitable; with=flts, boundary=:touches, progress=false)

    slp_count = Rasters.zonal(sum, one_slopes; of=dom.site_data.geom, progress=false)
    flt_count = Rasters.zonal(sum, one_flats; of=dom.site_data.geom, progress=false)

    total = total_aggregation.(slp_count, flt_count)
    return flt_count ./ total, slp_count ./ total
end

# wave UB
ub_mean::Vector{Float64} = zeros(Float64, 3806)
ub_std::Vector{Float64}  = zeros(Float64, 3806)

slp_ub_mean::Vector{Float64} = zeros(Float64, 3806)
slp_ub_std::Vector{Float64}  = zeros(Float64, 3806)

flt_ub_mean::Vector{Float64} = zeros(Float64, 3806)
flt_ub_std::Vector{Float64}  = zeros(Float64, 3806)

# wave HS
hs_mean::Vector{Float64} = zeros(Float64, 3806)
hs_std::Vector{Float64}  = zeros(Float64, 3806)

slp_hs_mean::Vector{Float64} = zeros(Float64, 3806)
slp_hs_std::Vector{Float64}  = zeros(Float64, 3806)

flt_hs_mean::Vector{Float64} = zeros(Float64, 3806)
flt_hs_std::Vector{Float64}  = zeros(Float64, 3806)

# wave Tp
tp_mean::Vector{Float64} = zeros(Float64, 3806)
tp_std::Vector{Float64}  = zeros(Float64, 3806)

slp_tp_mean::Vector{Float64} = zeros(Float64, 3806)
slp_tp_std::Vector{Float64}  = zeros(Float64, 3806)

flt_tp_mean::Vector{Float64} = zeros(Float64, 3806)
flt_tp_std::Vector{Float64}  = zeros(Float64, 3806)

# Bathymetry
bathy_mean::Vector{Float64} = zeros(Float64, 3806)
bathy_std::Vector{Float64}  = zeros(Float64, 3806)

slp_bathy_mean::Vector{Float64} = zeros(Float64, 3806)
slp_bathy_std::Vector{Float64}  = zeros(Float64, 3806)

flt_bathy_mean::Vector{Float64} = zeros(Float64, 3806)
flt_bathy_std::Vector{Float64}  = zeros(Float64, 3806)

# Turbidity
turbid_mean::Vector{Float64} = zeros(Float64, 3806)
turbid_std::Vector{Float64}  = zeros(Float64, 3806)

slp_turbid_mean::Vector{Float64} = zeros(Float64, 3806)
slp_turbid_std::Vector{Float64}  = zeros(Float64, 3806)

flt_turbid_mean::Vector{Float64} = zeros(Float64, 3806)
flt_turbid_std::Vector{Float64}  = zeros(Float64, 3806)

slp_prop::Vector{Float64} = zeros(Float64, 3806)
flt_prop::Vector{Float64} = zeros(Float64, 3806)

for reg_idx in 1:length(REGIONS)
    @info "Processing $(reg)"
    slopes_fn  = SLOPE_FNS[reg_idx]
    flats_fn   = FLAT_FNS[reg_idx]
    benthic_fn = BENTHIC_FNS[reg_idx]

    slopes  = GDF.read(slopes_fn)
    flats   = GDF.read(flats_fn)
    benthic = GDF.read(benthic_fn)

    habitable_mask = benthic[benthic.class .!= "Sand", :]

    waves_ub_fn = WAVES_UB_FNS[reg_idx]
    bathy_fn    = BATHY_FNS[reg_idx]
    turbid_fn   = TURBID_FNS[reg_idx]

    waves_ub = Raster(waves_ub_fn; lazy=true)
    waves_hs = Raster(waves_hs_fn; lazy=true)
    waves_tp = Raster(waves_tp_fn; lazy=true)
    bathy    = Raster(bathy_fn;    lazy=true)
    turbid   = Raster(turbid_fn;   lazy=true)

    raster_ones = copy(waves_hs)
    raster_ones .= 1

    @info "Processing $(reg): proportions"
    f_mean, s_mean = slope_flat_prop(slopes, flats, habitable_mask, raster_ones, dom)
    slp_prop .= clean_and_combine_vectors(slp_prop, s_mean)
    flt_prop .= clean_and_combine_vectors(flt_prop, f_mean)

    @info "Processing $(reg): waves_ub"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_ub, dom)
    ub_mean     .= clean_and_combine_vectors(ub_mean, mn)
    ub_std      .= clean_and_combine_vectors(ub_std,  st)
    slp_ub_mean .= clean_and_combine_vectors(slp_ub_mean, s_mean)
    slp_ub_std  .= clean_and_combine_vectors(slp_ub_std,  s_std)
    flt_ub_mean .= clean_and_combine_vectors(flt_ub_mean, f_mean)
    flt_ub_std  .= clean_and_combine_vectors(flt_ub_std,  f_std)

    @info "Processing $(reg): waves_hs"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_hs, dom)
    hs_mean     .= clean_and_combine_vectors(hs_mean, mn)
    hs_std      .= clean_and_combine_vectors(hs_std,  st)
    slp_hs_mean .= clean_and_combine_vectors(slp_hs_mean, s_mean)
    slp_hs_std  .= clean_and_combine_vectors(slp_hs_std,  s_std)
    flt_hs_mean .= clean_and_combine_vectors(flt_hs_mean, f_mean)
    flt_hs_std  .= clean_and_combine_vectors(flt_hs_std,  f_std)

    @info "Processing $(reg): waves_tp"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, waves_tp, dom)
    tp_mean .= clean_and_combine_vectors(tp_mean, mn)
    tp_std  .= clean_and_combine_vectors(tp_std,  st)
    slp_tp_mean .= clean_and_combine_vectors(slp_tp_mean, s_mean)
    slp_tp_std  .= clean_and_combine_vectors(slp_tp_std,  s_std)
    flt_tp_mean .= clean_and_combine_vectors(flt_tp_mean, f_mean)
    flt_tp_std  .= clean_and_combine_vectors(flt_tp_std,  f_std)

    @info "Processing $(reg): bathy"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, bathy, dom)
    bathy_mean .= clean_and_combine_vectors(bathy_mean, mn)
    bathy_std  .= clean_and_combine_vectors(bathy_std,  st)
    slp_bathy_mean .= clean_and_combine_vectors(slp_bathy_mean, s_mean)
    slp_bathy_std  .= clean_and_combine_vectors(slp_bathy_std,  s_std)
    flt_bathy_mean .= clean_and_combine_vectors(flt_bathy_mean, f_mean)
    flt_bathy_std  .= clean_and_combine_vectors(flt_bathy_std,  f_std)

    @info "Processing $(reg): turbid"
    f_mean, f_std, s_mean, s_std, mn, st = calc_stats(slopes, flats, habitable_mask, turbid, dom)
    turbid_mean .= clean_and_combine_vectors(turbid_mean, mn)
    turbid_std  .= clean_and_combine_vectors(turbid_std,  st)
    slp_turbid_mean .= clean_and_combine_vectors(slp_turbid_mean, s_mean)
    slp_turbid_std  .= clean_and_combine_vectors(slp_turbid_std,  s_std)
    flt_turbid_mean .= clean_and_combine_vectors(flt_turbid_mean, f_mean)
    flt_turbid_std  .= clean_and_combine_vectors(flt_turbid_std,  f_std)
end

# Visualise Distributions
using WGLMakie, GeoMakie, GraphMakie

f = Figure(; size=(1600, 1600))

# wave hs
var_name = "wave Ub"
col = :blue
Axis(f[1, 1]; xlabel=var_name, ylabel="location count", title="slope $(var_name) mean")
hist!(slp_ub_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[1, 2]; xlabel=var_name, ylabel="location count", title="flats $(var_name) mean")
hist!(flt_ub_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[2, 1]; xlabel=var_name, ylabel="location count", title="slope $(var_name) standard deviation")
hist!(slp_ub_std; strokewidth=1, strokecolor=:black, color=col)
Axis(f[2, 2]; xlabel=var_name, ylabel="location count", title="flats $(var_name) standard deviation")
hist!(flt_ub_std; strokewidth=1, strokecolor=:black, color=col)

# wave tp
var_name = "wave Tp"
col = :purple
Axis(f[3, 1]; xlabel=var_name, ylabel="location count", title="slope $(var_name) mean")
hist!(slp_tp_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[3, 2]; xlabel=var_name, ylabel="location count", title="flats $(var_name) mean")
hist!(flt_tp_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[4, 1]; xlabel=var_name, ylabel="location count", title="slope $(var_name) standard deviation")
hist!(slp_tp_std; strokewidth=1, strokecolor=:black, color=col)
Axis(f[4, 2]; xlabel=var_name, ylabel="location count", title="flats $(var_name) standard deviation")
hist!(flt_tp_std; strokewidth=1, strokecolor=:black, color=col)

# Bathy
var_name = "Bathy"
col = :red
Axis(f[1, 3]; xlabel=var_name, ylabel="location count", title="slope $(var_name) mean")
hist!(slp_bathy_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[1, 4]; xlabel=var_name, ylabel="location count", title="flats $(var_name) mean")
hist!(flt_bathy_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[2, 3]; xlabel=var_name, ylabel="location count", title="slope $(var_name) standard deviation")
hist!(slp_bathy_std; strokewidth=1, strokecolor=:black, color=col)
Axis(f[2, 4]; xlabel=var_name, ylabel="location count", title="flats $(var_name) standard deviation")
hist!(flt_bathy_std; strokewidth=1, strokecolor=:black, color=col)

# Turbid
var_name = "Turbid"
col = :orange
Axis(f[3, 3]; xlabel=var_name, ylabel="location count", title="slope $(var_name) mean")
hist!(slp_turbid_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[3, 4]; xlabel=var_name, ylabel="location count", title="flats $(var_name) mean")
hist!(flt_turbid_mean; strokewidth=1, strokecolor=:black, color=col)
Axis(f[4, 3]; xlabel=var_name, ylabel="location count", title="slope $(var_name) standard deviation")
hist!(slp_turbid_std; strokewidth=1, strokecolor=:black, color=col)
Axis(f[4, 4]; xlabel=var_name, ylabel="location count", title="flats $(var_name) standard deviation")
hist!(flt_turbid_std; strokewidth=1, strokecolor=:black, color=col)

using CSV, DataFrames

res = DataFrame(
    ub_mean=ub_mean,
    ub_std=ub_std,
    slp_ub_mean=slp_ub_mean,
    slp_ub_std=slp_ub_std,
    flt_ub_mean=flt_ub_mean,
    flt_ub_std=flt_ub_std,
    hs_mean=hs_mean,
    hs_std=hs_std,
    slp_hs_mean=slp_hs_mean,
    slp_hs_std=slp_hs_std,
    flt_hs_mean=flt_hs_mean,
    flt_hs_std=flt_hs_std,
    tp_mean=tp_mean,
    tp_std=tp_std,
    slp_tp_mean=slp_tp_mean,
    slp_tp_std=slp_tp_std,
    flt_tp_mean=flt_tp_mean,
    flt_tp_std=flt_tp_std,
    bathy_mean=bathy_mean,
    bathy_std=bathy_std,
    slp_bathy_mean=slp_bathy_mean,
    slp_bathy_std=slp_bathy_std,
    flt_bathy_mean=flt_bathy_mean,
    flt_bathy_std=flt_bathy_std,
    turbid_mean=turbid_mean,
    turbid_std=turbid_std,
    slp_turbid_mean=slp_turbid_mean,
    slp_turbid_std=slp_turbid_std,
    flt_turbid_mean=flt_turbid_mean,
    flt_turbid_std=flt_turbid_std,
    slp_prop=slp_prop,
    flt_prop=flt_prop
)
