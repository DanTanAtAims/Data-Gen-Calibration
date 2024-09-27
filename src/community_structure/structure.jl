import GeoDataFrames as GDF
import ArchGDAL as AG
using CSV, DataFrames, Glob, Statistics

include("../common.jl")

loc_rme_labels = "C:/Users/dtan/Downloads/coral_composition_locations.gpkg"
composition_files = readdir("C:/Users/dtan/repos/Reef-Monitoring-Endpoints/Outputs/")

loc_rme_gpkg = GDF.read(loc_rme_labels)

function _has_composition(location_name::String; composition_files=composition_files)::Bool
    for fl in composition_files
        contains(fl, location_name) && return true
    end
    return false
end

has_data_mask = _has_composition.(loc_rme_gpkg.reef_names)
usable_locs = loc_rme_gpkg[has_data_mask, :]

function _get_matching_files(reef_name; composition_files=composition_files)::Vector{String}
    fls = []
    for fl in composition_files
        if contains(fl, reef_name)
            push!(fls, fl)
        end
    end
    return fls
end

taxa_names = ["name"]

function reefmon_to_index(taxa::String)::Int64
    if taxa == "Acropora"
        return -1 # Split acro table corym
    end
    if taxa == "Pocilloporidae"
        return 3 # Corym non Acro
    end
    if taxa == "Isopora"
        return 2  # Acro Corym
    end
    if taxa == "Merulinidae"
        return 4 # Small massive
    end
    if taxa == "Other"
        return 0
    end
    if taxa == "Porites"
        return 5 # Large Massive
    end
    if taxa == "Unidentified"
        return 0
    end
    if taxa == "Lobophylliidae"
        return 4 # Small Massive
    end
    if taxa == "Montipora"
        return 2 # Acro Corym
    end
    if taxa == "Euphylliidae"
        return 4 # Small Massive
    end
    if taxa == "Dendrophylliidae"
        return 2 # Acro Corym
    end
    if taxa == "Fungiidae"
        return 4 # Small Massive
    end
    if taxa == "Goniopora Alveopora"
        return 5 # Large Massive
    end
    if taxa == "Pachyseris"
        return 4 # Small Massive
    end
    if taxa == "Rare groups"
        return 0
    end
    if taxa == "Leptastrea"
        return 4 # Small Massive
    end
    if taxa == "Psammocora"
        return 4 # Small Massive
    end
    if taxa == "Agariciidae"
        return 4 # Small Massive
    end
    return 0
end

function comp_to_adria_comp(nms::Vector{String}, comp_avs::Vector{Float64})::Vector{Float64}
    inds::Vector{Int64} = reefmon_to_index.(nms)
    output::Vector{Float64} = zeros(Float64, 5)
    for (idx, val) in zip(inds, comp_avs)
        if idx == 0
            output .+= val / 5
        elseif idx == -1
            output[1:2] .+= val / 2
        else
            output[idx] += val
        end
    end
    return output
end

function extract_mean_composition(
    locs::DataFrameRow;
    composition_files=composition_files
)::Vector{Float64}
    rf_name::String = locs.reef_names
    fns = _get_matching_files(rf_name; composition_files=composition_files)
    dfs = CSV.read.("C:/Users/dtan/repos/Reef-Monitoring-Endpoints/Outputs/" .* fns, DataFrame)
    yrs_range = 2006:2008
    comp = fill(0.0, 5)
    counts = 0
    if length(dfs) == 0
        return comp
    end
    for df in dfs
        append!(taxa_names, String.(names(df[:, 2:end])))
        mask = [y in yrs_range for y in df.Year]
        if count(mask) == 0
            continue
        end
        counts += 1
        comp .+= comp_to_adria_comp(String.(names(df[:, 2:end])), dropdims(mean(Matrix(df[mask, 2:end]), dims=1), dims=1))
    end
    if counts == 0
        return comp
    end
    comp = comp ./ counts
    return comp ./ sum(comp)
end

prealloc = zeros(Float64, 5, size(usable_locs, 1))

for (idx, rw) in enumerate(eachrow(usable_locs))
    prealloc[:, idx] .= extract_mean_composition(rw)
end

res = DataFrame(prealloc', [:tabular_acropora, :corymbose_acropora, :corymbose_non_acropora, :small_massives, :large_massives])
res[!, :geometry] .= usable_locs.geometry
res[!, :reef_name] .= usable_locs.reef_names
res[!, :RME_UNIQUE_ID] .= usable_locs.RME_UNIQUE_ID
res[!, :reef_a_sector] .= usable_locs.reef_a_sector
res[!, :reef_p_code] .= usable_locs.reef_p_code

function is_empty(comp::Vector{Float64})
    return all(comp .== 0.0)
end

empty_mask = is_empty.(Vector.(eachrow(res[:, 1:5])))
res = res[(!).(empty_mask), :]

classification_csv = CSV.read(OUTPUT_CSV, DataFrame)

res_idxs = [findfirst(x -> string(x) == id, classification_csv.RME_UNIQUE_ID) for id in res.RME_UNIQUE_ID]
res[!, :consecutive_classification] .= classification_csv[res_idxs, :consecutive_classification]

comp_classes = unique(res.consecutive_classification)
av_res = zeros(Float64, length(comp_classes), 5)
for (idx, cls) in enumerate(comp_classes)
    mask = res.consecutive_classification .== cls
    av_res[idx, :] .= dropdims(mean(Matrix(res[mask, 1:5]), dims=1), dims=1)
end

class_res = DataFrame(av_res, Symbol.(names(res)[1:5]))
class_res[!, :consecutive_classification] .= comp_classes

CSV.write("Outputs/class_community_structure.csv", class_res)
