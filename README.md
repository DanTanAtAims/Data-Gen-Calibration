# Calibration Data Generation

This repository contains scripts used to generate the data for the calibration of
[ADRIAmod](https://github.com/open-AIMS/ADRIA.jl)/[CoralBlox](https://github.com/open-AIMS/CoralBlox.jl).

## Setup

### Required Data and Path Configuration

Before executing any scripts, a `config.toml` file must be created in the repository root
directory. See the `example_config.toml` for the expected TOML format. Given the large range
of data required each section is able to run seperately if data requirements are not met for
all sections.

### Script Execution

The scripts are assumed to be executed from the repository root directory.

```julia-repl
julia> ]instantiate

julia> include("run.jl")
```

Scripts need not be executed at the same time using the top level `run.jl` script. Each
section has its own run and execution order. See below for more information.

## Location Classification

## Location Level Ground Truth

[comment]: <> (The following was written by Ben Grier)

Matching LTMP locations to ReefModEngine reef locations, where multiple sites are found in a
reef, these locations are aggregated with mean (including lat/lon). GBRMPA_ID and
RME_UNIQUE_ID are attached to output dataframes, as well as geometries for each site's point.

Lat/Lon mismatch at certain years/sites in the ltmp csv data. We have manually adjusted
these values to match the rest of the data.

Where an ltmp site does not overlap with a RME reef outline the site is joined to the
closest reef, if there is a reef within 1km distance.

Have to input the canonical-reefs geopackage dataset and the path to ltmp data directory.

It is assumed that the Lat/Lon values are recorded with a GPS device or similar and are in
EPSG:4326 CRS.

[comment]: <> (Contribution by Ben Grier end here)

### Script Execution

The Location classification script must be executed before this section.

```julia-repl
julia> include("src/reefmod_loc_obs/run.jl")
[ Info: Writing manta tow observations for reefmod location to <OUTPUT_PATH>
[ Info: Writing photo and video observations for reefmod location to <OUTPUT_PATH>
[ Info: Writing LTMP Photo Hard coral Observations At ReefMod Locations <OUTPUT_PATH>
[ Info: Writing LTMP Photo Hard and Soft coral Observations At ReefMod Locations <OUTPUT_PATH>
[ Info: Writing LTMP Manta Tow Coral Observations At ReefMod Locations to <OUTPUT_PATH>
```

## Classification Level Ground Truth


