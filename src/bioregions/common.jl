include("../common.jl")

BIOREGION_CONFIG = config_toml["Bioregions"]

bioregion_shp_path = BIOREGION_CONFIG["bioregion_shp_path"]

OUTPUT_DOMAIN_BIOREGION_GPKG = joinpath(OUTPUT_DIR, "domain_bioregions.gpkg")
