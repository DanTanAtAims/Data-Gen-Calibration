"""
Define functions to classify ecorrap species classification to C~scape/ADRIAmod functional groups.

Defines
- is_tablular_Acroproa
- is_corymbose_Acropora
- is_corymbose_non_Acropora
- is_small_massive
- is_large_massive
"""

include("common.jl")

using CSV, DataFrames

const CScape_to_ADRIA::Dict{String, String} = Dict(
    :acro_table => :tabular_Acropora,
    :corymbose_Acropora => :corymbose_Acropora,
    :corym_non_acro => :corymbose_non_Acropora,
    :small_massive => :small_massives,
    :large_massives => :large_massives
)

global classification_csv = CSV.read(CORAL_CLASSIFICATION_PATH, DataFrame)

"""
    species_code_to_cscape(code::String)::String

Convert a ecorrap species code to a cscape functional group name.
"""
function _species_code_to_cscape(code::String)::String
    code_matches = classification_csv.Code .== code

    n_matches = count(code_matches)
    if n_matches > 1
        msg  = "Species code: $(code)"
        msg *= ", matches multiple codes: $(classification_csv.Code[code_matches]). "
        msg *= "Returning first match."
        @warn msg
    elseif n_matches == 0
        msg  = "Given species code, $(code), "
        msg *= "does not match any codes found in classification DataFrame."
        throw(ArgumentError(msg))
    end

    return classification_csv.Cscape_group[code_matches][1]
end

function _is_functional_group(code::String, target_cscape_group::String)::Bool
    cscape_group = species_code_to_cscape(code)
    return cscape_group == target_cscape_group
end

# Classification Functions
"""
    is_tabular_Acropora(code::String)::Bool
"""
function is_tabular_Acropora(code::String)::Bool
    # juv_quadrats do not distinguish between tabular and corymbose Acropora
    if code == "Acropora"
        return true
    end

    return _is_functional_group(code, "acro_table")
end
"""
    is_corymbose_Acropora(code::String)::Bool
"""
function is_corymbose_Acropora(code::String)::Bool
    # juv_quadrats do not distinguish between tabular and corymbose Acropora
    if code == "Acropora"
        return true
    end

    return _is_functional_group(code, "acro_corym")
end
"""
    is_corymbose_non_Acropora(code::String)::Bool
"""
function is_corymbose_non_Acropora(code::String)::Bool
    return _is_functional_group(code, "corym_non_acro")
end
"""
    is_small_massive(code::String)::Bool
"""
function is_small_massive(code::String)::Bool
    return _is_functional_group(code, "small_massive")
end
"""
    is_large_massive(code::String)::Bool
"""
function is_large_massive(code::String)::Bool
    return _is_functional_group(code, "large_massive")
end
