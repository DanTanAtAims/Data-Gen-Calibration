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
related to the calculation of growth statistics and add diameter and linear extension
columns.
"""
function get_growth_entries(raw_data::DataFrame)::DataFrame
    # Construct masks to remove unused and missing data
    growth_mask = raw_data.to_use_for_growth .== "yes"
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

    # Add growth and linear extension entries into data frame
    growth_data[!, :growth]  .= growth_data.sizeNext .- growth_data.size
    growth_data[!, :lin_ext] .= growth_data.diamNext .= growth_data.diam

    return growth_data
end

"""
    get_survival_entries(raw_data:;DataFrame)::DataFrame

Given the csv from containing all entries of the coral demograph data, remove rows not
related to the calculation of survival statistics and add diameter column.
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

    return survival_data
end
