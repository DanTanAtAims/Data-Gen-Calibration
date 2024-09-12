include("data_matching.jl")

# Process ltmp data

ltmp_data = CSV.read(LTMP_PHOTO_FN, DataFrame)
gd = DataFrames.groupby(ltmp_data, [:REEF_ID, :YEAR_CODE, :GROUP_CODE])
ltmp_reef_level = combine(gd, [:LATITUDE, :LONGITUDE, :COVER] .=> mean)

# Particular years have incorrect lon/lat at some REEF_IDs. Need to standardise them for later
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "16028S", :LATITUDE_mean] .= -16.3867
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "16028S", :LONGITUDE_mean] .= 145.572
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "14139S", :LONGITUDE_mean] .= 145.649
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "14139S", :LATITUDE_mean] .= -14.6297
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "22088S", :LATITUDE_mean] .= -22.0273
ltmp_reef_level[ltmp_reef_level.REEF_ID .== "22088S", :LONGITUDE_mean] .= 152.193

ltmp_reefs = unique(ltmp_reef_level[:, [:REEF_ID, :LATITUDE_mean, :LONGITUDE_mean]])
ltmp_reefs.geometry = Vector{Union{Missing, ArchGDAL.IGeometry{ArchGDAL.wkbPoint}}}(missing, size(ltmp_reefs, 1))
for row in eachrow(ltmp_reefs)
    row.geometry = ArchGDAL.createpoint()
    ArchGDAL.addpoint!(row.geometry, row.LONGITUDE_mean, row.LATITUDE_mean)
end

matching_labels = find_intersections(ltmp_reefs, rme_reefs, :REEF_ID, :GBRMPA_ID, "km", 1; nearest=true)
rename!(matching_labels, :area_ID => :GBRMPA_ID)
ltmp_reefs = leftjoin(ltmp_reefs, matching_labels, on=:REEF_ID, order=:left)
ltmp_reef_level = leftjoin(ltmp_reef_level, ltmp_reefs[:, [:REEF_ID, :geometry, :GBRMPA_ID]]; on=:REEF_ID, order=:left)
ltmp_reef_level = leftjoin(ltmp_reef_level, rme_reefs[:, [:GBRMPA_ID, :RME_UNIQUE_ID]], on=:GBRMPA_ID, order=:left, matchmissing=:notequal)

ltmp_reef_level.GBRMPA_ID = Vector{Union{Missing, String}}(ltmp_reef_level.GBRMPA_ID)
ltmp_reef_level.REEF_ID = String.(ltmp_reef_level.REEF_ID)
ltmp_reef_level.GROUP_CODE = String.(ltmp_reef_level.GROUP_CODE)

GDF.write(OUT_RME_PHOTO, ltmp_reef_level; crs=GFT.EPSG(4326))

# Process ltmp-hard-coral cover data

ltmp_reef_hc = ltmp_reef_level[ltmp_reef_level.GROUP_CODE .== "Hard Coral", :]
ltmp_reef_hc.YEAR .= ""
ltmp_reef_hc.YEAR_CODE = string.(ltmp_reef_hc.YEAR_CODE)
for row in eachrow(ltmp_reef_hc)
    year_code = row.YEAR_CODE
    row.YEAR = row.YEAR_CODE[1:4]
end

select!(ltmp_reef_hc, [:REEF_ID, :LATITUDE_mean, :LONGITUDE_mean, :COVER_mean, :YEAR, :geometry, :GBRMPA_ID, :RME_UNIQUE_ID])
ltmp_reef_hc = unstack(ltmp_reef_hc, :YEAR, :COVER_mean)

@info "Writing LTMP Photo Hard coral Observations At ReefMod Locations"
GDF.write(OUT_RME_PHOTO_HC, ltmp_reef_hc; crs=GFT.EPSG(4326))


# Process ltmp-total-cover data
ltmp_hc_sc = ltmp_reef_level[ltmp_reef_level.GROUP_CODE .∈ [["Hard Coral", "Soft Coral"]], :]
gdf = DataFrames.groupby(ltmp_hc_sc, [:REEF_ID, :YEAR_CODE, :LATITUDE_mean, :LONGITUDE_mean, :geometry, :GBRMPA_ID, :RME_UNIQUE_ID])
ltmp_hc_sc = combine(gdf, :COVER_mean .=> sum)

ltmp_hc_sc.YEAR .= ""
ltmp_hc_sc.YEAR_CODE = string.(ltmp_hc_sc.YEAR_CODE)
for row in eachrow(ltmp_hc_sc)
    year_code = row.YEAR_CODE
    row.YEAR = row.YEAR_CODE[1:4]
end

select!(ltmp_hc_sc, [:REEF_ID, :LATITUDE_mean, :LONGITUDE_mean, :COVER_mean_sum, :YEAR, :geometry, :GBRMPA_ID, :RME_UNIQUE_ID])
ltmp_hc_sc = unstack(ltmp_hc_sc, :YEAR, :COVER_mean_sum)

@info "Writing LTMP Photo Hard and Soft coral Observations At ReefMod Locations"
GDF.write(OUT_RME_PHOTO_HC_SC, ltmp_hc_sc; crs=GFT.EPSG(4326))
