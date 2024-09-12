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

OUTPUT_DIR = joinpath(@__DIR__, "..", "Outputs")
@assert isdir(OUTPUT_DIR) "Output directory, $(OUTPUT_DIR), does not exist."

OUTPUT_CSV   = joinpath(OUTPUT_DIR, "location_classification.csv")

# LTMP ReefMod Locations
OUT_PHOTO_NETCDF = joinpath(OUTPUT_DIR, "ltmp_vid_photo.nc")
OUT_MANTA_NETCDF = joinpath(OUTPUT_DIR, "ltmp_manta_tow.nc")

OUT_RME_MANTA       = joinpath(OUTPUT_DIR, "ltmp_reefmod_manta.gpkg")
OUT_RME_PHOTO       = joinpath(OUTPUT_DIR, "ltmp_reefmod_photo.gpkg")
OUT_RME_PHOTO_HC    = joinpath(OUTPUT_DIR, "ltmp_reefmod_photo_hc.gpkg")
OUT_RME_PHOTO_HC_SC = joinpath(OUTPUT_DIR, "ltmp_reefmod_photo_hc_sc.gpkg")

OUT_CLASS_MANTA    = joinpath(OUTPUT_DIR, "manta_tow_mean_std.nc")
OUT_CLASS_PHOTO_HC = joinpath(OUTPUT_DIR, "ltmp_hard_mean_std.nc")
