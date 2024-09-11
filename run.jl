@info "Generating Location Classifications"
include("src/classification/run.jl")
@info "Generating LTMP Observations for Reefmod Locations (empty)"
include("src/reefmod_loc_obs/run.jl")
@info "Generating LTMP Observations for Location Classifications (empty)"
include("src/class_obs/run.jl")
