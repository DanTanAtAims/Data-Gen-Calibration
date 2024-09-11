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

```julia
julia> ]instantiate

julia> include("run.jl")
```

Scripts need not be executed at the same time using the top level `run.jl` script. Each
section has its own run and execution order. See below for more information.

## Location Classification

## Location Level Ground Truth

## Classification Level Ground Truth


