include("data_matching.jl")

# Process manta tow data.

manta_data = open_dataset(OUT_MANTA_NETCDF)
manta_live_coral = manta_data.MEAN_LIVE_CORAL

manta_LC_df = DataFrame(manta_live_coral.data[:,:], collect(getAxis("locs", manta_live_coral).val))
manta_LC_df.timesteps = collect(getAxis("timesteps", manta_live_coral).val)
select!(manta_LC_df, :timesteps, Not(:timesteps))
manta_LC_df.timesteps = string.(manta_LC_df.timesteps)
manta_LC_df = permutedims(manta_LC_df, "timesteps", "REEF_ID")

manta_tow_reefs = CSV.read(LTMP_MANTA_FN, DataFrame)
manta_tow_reefs = unique(manta_tow_reefs[:, [:REEF_ID, :LATITUDE, :LONGITUDE]])
manta_tow_reefs.geometry = Vector{Any}(missing, size(manta_tow_reefs, 1))
for row in eachrow(manta_tow_reefs)
    row.geometry = ArchGDAL.createpoint()
    ArchGDAL.addpoint!(row.geometry, row.LONGITUDE, row.LATITUDE)
end

matching_labels = find_intersections(manta_tow_reefs, rme_reefs, :REEF_ID, :GBRMPA_ID, "km", 1; nearest=true)
rename!(matching_labels, :area_ID => :GBRMPA_ID)

manta_tow_reefs = leftjoin(manta_tow_reefs, matching_labels; on=:REEF_ID, order=:left)
manta_LC_df = leftjoin(manta_LC_df, manta_tow_reefs[:, [:REEF_ID, :GBRMPA_ID, :geometry]], on=:REEF_ID, order=:left)
manta_LC_df = leftjoin(manta_LC_df, rme_reefs[:, [:GBRMPA_ID, :RME_UNIQUE_ID]], on=:GBRMPA_ID, order=:left, matchmissing=:notequal)
select!(manta_LC_df, [:REEF_ID, :GBRMPA_ID, :RME_UNIQUE_ID], Not([:REEF_ID, :GBRMPA_ID, :RME_UNIQUE_ID]))
manta_LC_df.GBRMPA_ID = Vector{Union{Missing, String}}(manta_LC_df.GBRMPA_ID)

@info "Writing LTMP Manta Tow Coral Observations At ReefMod Locations"
GDF.write(OUT_RME_MANTA, manta_LC_df; crs=GFT.EPSG(4326))
