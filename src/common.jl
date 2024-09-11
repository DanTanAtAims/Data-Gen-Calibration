using TOML

function assert_file_exist(filepath::String)::Nothing
    @assert isfile(filepath) "File: $(filepath) does not exist."
    return nothing
end

function assert_config_key_exists(config_dict, key)::Nothing
    msg  = "Configuration key: $(key) does not exist."
    msg *= " Check that all required file paths are provided."
    @assert haskey(config_dict, key)
    return nothing
end

config_fn = joinpath(@__DIR__, "../config.toml")
if !isfile(config_fn)
    throw(ArgumentError(
        "Configuration file not found. File \'config.toml\' not found in repository root."
    ))
end

config_toml = TOML.parsefile(config_fn)

COMMON_CONFIG    = config_toml["Common"]
RME_DOMAIN_DIR   = COMMON_CONFIG["rme_domain_path"]
# LOC_OBS_CONFIG   = config_toml["reefmod_loc_obs"]
# 
# CANONICAL_GPKG_FN = LOC_OBS_CONFIG["canonical_gpkg_path"]
# MANTA_TOW_FN      = LOC_OBS_CONFIG["ltmp_manta_tow_path"]
# 
# assert_file_exist(CANONICAL_GPKG_FN)
# assert_file_exist(MANTA_TOW_FN)

OUTPUT_DIR = joinpath(@__DIR__, "../Outputs")
@assert isdir(OUTPUT_DIR) "Output directory, $(OUTPUT_DIR), does not exist."

OUTPUT_CSV = joinpath(OUTPUT_DIR, "location_classification.csv")
