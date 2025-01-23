"""
Define functions that extract data rows to be used for growht and mortality. These functions
should also convert the type of strings variables to usable types such as booleans, ints and
floats where possible.
"""

using DataFrames

"""
    _area_to_diam(area::Float64)::Float64

Calculate the diameter of a coral given its area by assuming it is a circle.
"""
function _area_to_diam(area::Float64)::Float64
    return sqrt(4 * area / pi)
end

"""
    get_growth_entries(raw_data::DataFrame)::DataFrame

Given the csv from containing all entries of the coral demograph data, remove rows not
related to the calculation of growth statistics and add diameter, log diameter and linear
extension columns.
"""
function get_growth_entries(raw_data::DataFrame)::DataFrame
    # Maintain compatanility with previous yearly data versions
    if "GROWTH_USE" in names(raw_data)
        rename_map = [
            "GROWTH_USE" => "to_use_for_growth",
            "SURVIVAL" => "surv",
            "AREA_T1_SQCM" => "size",
            "AREA_T2_SQCM" => "sizeNext",
            "CLUSTER" => "Cluster",
            "REEF" => "Reef",
            "SITE" => "Site_UID",
            "TAXON" => "Taxa"
        ]
        rename!(raw_data, rename_map...)

        surv_na_mask = String.(raw_data.surv) .!= "NA"
        growth_na_mask = String.(raw_data.size) .!= "NA"
        raw_data = raw_data[surv_na_mask .&& growth_na_mask, :]

        raw_data[!, :surv] .= String.(raw_data.surv)
        raw_data[!, :surv] .= parse.(Int64, raw_data.surv)

        # Cast sizeNext column to Float64
        raw_data[!, :size] .= String.(raw_data.size)
        raw_data[!, :size] .= parse.(Float64, raw_data.size)
    end
    # Construct masks to remove unused and missing data
    growth_mask = raw_data.to_use_for_growth.== "yes"
    survived_mask = raw_data.surv .== 1
    non_missing_size_mask = raw_data.size .!= "NA"

    # Remove missing and unused data
    growth_data::DataFrame = raw_data[
        growth_mask .&& survived_mask .&& non_missing_size_mask, :
    ]

    # Cast sizeNext column to Float64
    growth_data[!, :sizeNext] .= String.(growth_data.sizeNext)
    growth_data[!, :sizeNext] .= parse.(Float64, growth_data.sizeNext)

    # Add diameter column and diameter Next column
    growth_data[!, :diam]     .= _area_to_diam.(growth_data.size)
    growth_data[!, :diamNext] .= _area_to_diam.(growth_data.sizeNext)

    # Add log diameter column
    growth_data[!, :logdiam]  .= log.(2, growth_data.diam)

    # Add growth and linear extension entries into data frame
    growth_data[!, :growth]  .= growth_data.sizeNext .- growth_data.size
    growth_data[!, :lin_ext] .= growth_data.diamNext .- growth_data.diam

    # Cast taxa String15 type to string type
    growth_data[!, :Taxa] .= String.(growth_data.Taxa)

    no_partial_mask = growth_data.lin_ext .> 0.0
    growth_data = growth_data[no_partial_mask, :]

    return growth_data
end

"""
    get_survival_entries(raw_data:;DataFrame)::DataFrame

Given the csv from containing all entries of the coral demograph data, remove rows not
related to the calculation of survival statistics and add diameter and log diameter columns.
"""
function get_survival_entries(raw_data::DataFrame)::DataFrame
    # Construct masks to remove unused and missing data
    for_survival = raw_data[:, Symbol("To_use_for_survival_21.22")] .== "yes"
    non_missing_size_mask = raw_data.size .!= "NA"

    # Remove missing and unused data
    survival_data::DataFrame = raw_data[
        for_survival .&& non_missing_size_mask, :
    ]

    # insert diameter column
    survival_data[!, :diam] .= _area_to_diam.(survival_data.size)

    # Add log diameter column
    survival_data[!, :logdiam]  .= log.(2, survival_data.diam)

    # Cast taxa String15 type to string type
    survival_data[!, :Taxa] .= String.(survival_data.Taxa)

    return survival_data
end
