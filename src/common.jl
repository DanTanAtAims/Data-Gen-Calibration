using TOML

config_fn = joinpath(@__DIR__, "../config.toml")
if !isfile(config_fn)
    throw(LoadError(
        "Configuration file not found. File \'config.toml\' not found in repository root."
    ))
end

config_toml = TOML.parsefile(config_fn)

COMMON_CONFIG    = config_toml["Common"]
RME_DOMAIN_DIR   = COMMON_CONFIG["rme_domain_path"]

# Location Classification Paths
LOC_CLASS_CONFIG = config_toml["location_classification"]

# ACA data is only used for slope and flat categorisation
ACA_GEOSPATIAL_DIR = LOC_CLASS_CONFIG["aca_geospatial_dir"]
# MPA data is used for all other geospatial data
MPA_GEOSPATIAL_DIR = LOC_CLASS_CONFIG["mpa_geospatial_dir"]

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
    joinpath(ACA_GEOSPATIAL_DIR, "aca_target_benthic_" * reg * ".gpkg") for reg in REGIONS
]

WAVES_UB_FNS = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_waves_Ub.tif") for reg in REGIONS
]
BATHY_FNS    = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_bathy.tif")    for reg in REGIONS
]
TURBID_FNS   = [
    joinpath(MPA_GEOSPATIAL_DIR, reg * "_turbid.tif")   for reg in REGIONS
]

# Validating Loccation Classification File Names
function assert_file_exist(filepath::String)::Nothing
    @assert isfile(filepath) "File: $(filepath) does not exist."
    return nothing
end

assert_file_exist.(SLOPE_FNS)
assert_file_exist.(FLAT_FNS)
assert_file_exist.(BENTHIC_FNS)
assert_file_exist.(WAVES_UB_FNS)
assert_file_exist.(BATHY_FNS)
assert_file_exist.(TURBID_FNS)

# ReefMod Location Observation Path Construction
LOC_OBS_CONFIG   = config_toml["reefmod_loc_obs"]

CANONICAL_GPKG_FN = LOC_OBS_CONFIG["canonical_gpkg_path"]
MANTA_TOW_FN      = LOC_OBS_CONFIG["ltmp_manta_tow_path"]

assert_file_exist(CANONICAL_GPKG_FN)
assert_file_exist(MANTA_TOW_FN)

OUTPUT_CONFIG = config_toml["Output"]
OUTPUT_DIR    = OUTPUT_CONFIG["output_dir"]

@assert isdir(OUTPUT_DIR) "Output directory, $(OUTPUT_DIR), does not exist."
