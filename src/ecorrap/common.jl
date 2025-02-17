include("../common.jl")

ECORRAP_CONFIG =  config_toml["EcoRRAP"]

# Input data paths
CORAL_OBS_PATH            = ECORRAP_CONFIG["coral_observation_path"]
CORAL_CLASSIFICATION_PATH = ECORRAP_CONFIG["coral_classification_path"]

# Output data paths
OUT_PLOT_DIR = joinpath(@__DIR__, "..", "..", "Outputs", "Plots")

ECORRAP_SAVE_PATH = joinpath(OUTPUT_DIR, "ecorrap_params.nc")
