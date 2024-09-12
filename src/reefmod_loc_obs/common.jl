include("../common.jl")

if !haskey(config_toml, "ReefMod_Loc_Obs")
    throw(ArgumentError(
        "Unable to find \"ReefMod_Loc_Obs\" section in config.toml"
    ))
end

REEFMOD_LOC_CONFIG = config_toml["ReefMod_Loc_Obs"]

assert_config_key_exists(REEFMOD_LOC_CONFIG, "ltmp_vid_photo_path")
assert_config_key_exists(REEFMOD_LOC_CONFIG, "ltmp_manta_tow_path")

LTMP_PHOTO_FN = REEFMOD_LOC_CONFIG["ltmp_vid_photo_path"]
LTMP_MANTA_FN = REEFMOD_LOC_CONFIG["ltmp_manta_tow_path"]

assert_file_exist(LTMP_PHOTO_FN)
assert_file_exist(LTMP_MANTA_FN)
