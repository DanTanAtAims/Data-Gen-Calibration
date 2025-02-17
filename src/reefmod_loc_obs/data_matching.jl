"""
Matching RME reef locations to the manta_tow_data survey results.

Written By Ben Grier and modified by Daniel Tan.
"""

include("common.jl")

using NetCDF, YAXArrays, CSV, DataFrames, ArchGDAL, Statistics
import GeoDataFrames as GDF
import GeoFormatTypes as GFT

function format_canonical_ltmp_id(canonical_ltmp_id::String)::String
    return uppercase(replace(canonical_ltmp_id, "-" => ""))
end

function is_same_reef(canonical_id::String, ltmp_id::AbstractString)::Bool
    formatted_canonical::String = format_canonical_ltmp_id(canonical_id)
    if formatted_canonical == ltmp_id
        return true
    end
    # find the beginning of the sub-reef indicator
    canonical_sub_ind = findfirst(isletter, formatted_canonical)
    canonical_sub_ind = isnothing(
        canonical_sub_ind
    ) ? length(formatted_canonical) : canonical_sub_ind
    ltmp_sub_ind = findfirst(isletter, ltmp_id)
    ltmp_sub_ind = isnothing(
        ltmp_sub_ind
    ) ? length(ltmp_id) : ltmp_sub_ind
    return formatted_canonical[1:canonical_sub_ind-1] == ltmp_id[1:ltmp_sub_ind-1]
end

function is_same_sub_reef(canonical_id::String, raw_ltmp_id::AbstractString)::Bool
    formatted_canonical::String = format_canonical_ltmp_id(canonical_id)
    if formatted_canonical == raw_ltmp_id
        return true
    end

    return isletter(formatted_canonical[end]) && formatted_canonical == raw_ltmp_id[1:end-1]
end

function is_same_id(canonical_id::String, ltmp_id::AbstractString)::Bool
    formatted_canonical::String = format_canonical_ltmp_id(canonical_id)
    return formatted_canonical == ltmp_id
end

function find_intersection(ltmp::DataFrame, canonical::DataFrame)::DataFrame
    matched_locs = DataFrame(
        [Vector{Any}(missing, size(ltmp, 1)) for _ in 1:3],
        [:REEF_ID, :area_ID, :match_reason])
    for (x_i, ltmp_row) in enumerate(eachrow(ltmp))
        for (y_i, canonical_row) in enumerate(eachrow(canonical))
            if is_same_id(canonical_row.LTMP_ID, ltmp_row.REEF_ID)
                matched_locs[x_i, :] .= [ltmp_row.REEF_ID, canonical_row.GBRMPA_ID, "Same ID"]
                break
            end
        end
        # If matched location, stop searching
        !ismissing(matched_locs[x_i, :REEF_ID]) && continue
        for (y_i, canonical_row) in enumerate(eachrow(canonical))
            if ismissing(matched_locs[x_i, 1]) &&
                is_same_sub_reef(canonical_row.LTMP_ID, ltmp_row.REEF_ID)
                matched_locs[x_i, :] .= [ltmp_row.REEF_ID, canonical_row.GBRMPA_ID, "Same Sub Reef"]
                break
            end
        end
        # If matched location, stop searching
        !ismissing(matched_locs[x_i, :REEF_ID]) && continue
        for (y_i, canonical_row) in enumerate(eachrow(canonical))
            if ismissing(matched_locs[x_i, 1]) &&
                is_same_reef(canonical_row.LTMP_ID, ltmp_row.REEF_ID)
                matched_locs[x_i, :] .= [ltmp_row.REEF_ID, canonical_row.GBRMPA_ID, "Same Reef"]
                break
            end
        end
        # If matched location, stop searching
        !ismissing(matched_locs[x_i, :REEF_ID]) && continue
        dists = [
            ArchGDAL.distance(ltmp_row.geometry, canon_geom)
            for canon_geom in canonical.geometry
        ]
        closest = argmin(dists)
        # Check closest location is within 1 km
        if dists[closest] * 111.1 < 1.0
            matched_locs[x_i, :] = [ltmp_row.REEF_ID, canonical.GBRMPA_ID[closest], "Closest Reef"]
        end
    end
    return matched_locs
end
"""
    find_intersections(
        x::DataFrame,
        y::DataFrame,
        x_id::Symbol,
        y_id::Symbol,
        units::String,
        threshold,
        y_geom_col::Symbol=:geometry;
        proportion::Bool=false,
        nearest::Bool=false
    )::DataFrame

Find the areas of `y` that intersect with each polygon in `x`.
`rel_areas` contains corresponding `y_id` for each intersecting polygon in x (can then be
joined to `x`).

If `proportion = true`: polygons of `y` are only chosen if the intersection with `x` is >
50% the area of `x`.

# Arguments
- `x` : The target GeoDataFrame to compare with
- `y` : GeoDataFrame containing polygons to match against `x`
- `xid` : Column name holding unique IDs for x geometries (referred to as GBRMPA_ID in rel_areas)
- `yid` : Column name holding variable of interest for y geometries
- `units` : String indicating the units to use for the nearest distance threshold
- `threshold` : Value to use for threshold for joining to nearest object when no overlap is found.
- `y_geom_col` : Column name holding geometries in y
- `proportion` : Only select y polygons if the intersection with x polygon is > 50% of x polygon area
                 (default: `false`).
- `nearest` : When there is no overlap join to nearest geometry instead of returning missing
                (within threshold), only specify units and threshold when nearest=true.
"""
function find_intersections(
    x::DataFrame,
    y::DataFrame,
    x_id::Symbol,
    y_id::Symbol,
    units::String,
    threshold,
    y_geom_col::Symbol=:geometry;
    proportion::Bool=false,
    nearest::Bool=false
)::DataFrame
    rel_areas = DataFrame(
        [Vector{Any}(missing, size(x, 1)) for _ in 1:2],
        [x_id, :area_ID]
    )
        for (x_i, reef_poly) in enumerate(eachrow(x))
            intersecting = DataFrame(
                [Vector{Any}(missing, size(y, 1)) for _ in 1:3],
                [x_id, :area_ID, :inter_area]
            )

            for (y_i, interest_area) in enumerate(eachrow(y))
                same_id::Bool = is_same_id(interest_area[y_id], String(reef_poly[x_id]))
                if same_id
                    data = [reef_poly[x_id], interest_area[y_id], 0.0]
                elseif ArchGDAL.intersects(reef_poly.geometry, interest_area[y_geom_col])
                    inter_area = ArchGDAL.intersection(
                        reef_poly.geometry, interest_area[y_geom_col]
                    )

                    inter_area = ArchGDAL.geomarea(inter_area)
                    if proportion
                        prop_area = inter_area / ArchGDAL.geomarea(reef_poly.geometry)

                        if prop_area >= 0.5
                            data = [reef_poly[x_id], interest_area[y_id], inter_area]

                        else
                            data = [missing, missing, missing]
                        end
                    else
                        data = [reef_poly[x_id], interest_area[y_id], inter_area]
                    end
                else
                    data = [reef_poly[x_id], missing, missing]
                end

                intersecting[y_i, :] = data
            end

            if all(ismissing, intersecting.area_ID)
                if nearest
                    distances = ArchGDAL.distance.([reef_poly.geometry], y[:, y_geom_col])

                    if units == "km"
                        distances .*= 111.1
                    end
                    if minimum(distances) < threshold
                        x_data = [intersecting[1, x_id], y[argmin(distances), y_id]]
                    else
                        x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
                    end
                else
                    x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
                end
            else
                dropmissing!(intersecting)
                max_inter_area = argmax(intersecting.inter_area)
                x_data = [intersecting[max_inter_area, x_id], intersecting[max_inter_area, :area_ID]]
            end
            rel_areas[x_i, :] = x_data
        end

    return rel_areas
end

function nonunique2(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end

# Insert directories for canonical-reefs geopackage and ltmp data directory.
rme_reefs = GDF.read(CANONICAL_GPKG_FN)
