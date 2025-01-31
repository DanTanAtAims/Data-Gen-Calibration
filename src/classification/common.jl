include("../common.jl")

if !haskey(config_toml, "Location_Classification")
    throw(ArgumentError(
        "Unable to find \"Location_Classification\" section in config.toml"
    ))
end

# Location Classification Paths
LOC_CLASS_CONFIG = config_toml["Location_Classification"]

assert_config_key_exists(LOC_CLASS_CONFIG, "aca_geospatial_dir")
assert_config_key_exists(LOC_CLASS_CONFIG, "mpa_geospatial_dir")

# ACA data is only used for slope and flat categorisation
ACA_GEOSPATIAL_DIR = LOC_CLASS_CONFIG["aca_geospatial_dir"]
# MPA data is used for all other geospatial data
MPA_GEOSPATIAL_DIR = LOC_CLASS_CONFIG["mpa_geospatial_dir"]
# Reprojected Waves Ubed90 data from the Callaghan dataset
WAVES_UBED90_DIR = LOC_CLASS_CONFIG["waves_ubed90_dir"]

@assert isdir(ACA_GEOSPATIAL_DIR) "ACA geospatial directory not found."
@assert isdir(MPA_GEOSPATIAL_DIR) "MPA geospatial directory not found."

REGIONS = [
    "Cairns-Cooktown",
    "FarNorthern",
    "Mackay-Capricorn",
    "Townsville-Whitsunday"
]

SLOPE_FNS = [
    joinpath(ACA_GEOSPATIAL_DIR, "aca_target_slopes_"  * reg * ".gpkg") for reg in REGIONS
]
FLAT_FNS = [
    joinpath(ACA_GEOSPATIAL_DIR, "aca_target_flats_"   * reg * ".gpkg") for reg in REGIONS
]
BENTHIC_FNS = [
    joinpath(ACA_GEOSPATIAL_DIR, "aca_benthic_" * reg * ".gpkg") for reg in REGIONS
]

WAVES_UB_FNS = [
    joinpath(WAVES_UBED90_DIR, reg * "_waves_Ub.tif") for reg in REGIONS
]
WAVES_HS_FNS = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_waves_Hs.tif") for reg in REGIONS
]
WAVES_TP_FNS = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_waves_Tp.tif") for reg in REGIONS
]
BATHY_FNS    = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_bathy.tif")    for reg in REGIONS
]
TURBID_FNS   = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_turbid.tif")   for reg in REGIONS
]

assert_file_exist.(SLOPE_FNS)
assert_file_exist.(FLAT_FNS)
assert_file_exist.(BENTHIC_FNS)
assert_file_exist.(WAVES_UB_FNS)
assert_file_exist.(BATHY_FNS)
assert_file_exist.(TURBID_FNS)

assert_config_key_exists(LOC_CLASS_CONFIG, "outer_shelf_gpkg")
assert_config_key_exists(LOC_CLASS_CONFIG, "inner_shelf_gpkg")

OUTER_SHELF_GPKG = LOC_CLASS_CONFIG["outer_shelf_gpkg"]
INNER_SHELF_GPKG = LOC_CLASS_CONFIG["inner_shelf_gpkg"]
