"""
Matching RME reef locations to the manta_tow_data survey results.
"""

using NetCDF, YAXArrays, CSV, DataFrames, ArchGDAL, Statistics
import GeoDataFrames as GDF
import GeoFormatTypes as GFT

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
                if ArchGDAL.intersects(reef_poly.geometry, interest_area[y_geom_col])
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
rme_reefs = GDF.read("../canonical-reefs/output/rrap_canonical_2024-07-24-T12-38-38.gpkg")

DATA_DIR = "c:/Users/bgrier/Documents/Projects/ltmp_reefmod_comb/"
