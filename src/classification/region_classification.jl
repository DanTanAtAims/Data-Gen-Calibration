"""
For each ReefMod Location, idenfity the shelf position and ltmp region.

# Shelf Position
1 - Inner Shelf
2 - Middle Shelf
3 - Outer Shelf

# LTMP Region
1 - North
2 - Central
3 - South
"""

include("common.jl")

using ADRIA
using ADRIA: GDF

using CSV, DataFrames

import ArchGDAL as AG

# Data Loading
inner_shelf = GDF.read(INNER_SHELF_GPKG)
outer_shelf = GDF.read(OUTER_SHELF_GPKG)

# Domain loading can be slow. Only load domain if it is not already loaded.
if !isdefined(Main, :dom)
    dom = ADRIA.load_domain(RMEDomain, "C:\\Users\\dtan\\repos\\rme_ml_2024_01_08", "45")
end

features = CSV.read(OUTPUT_CSV, DataFrame)

# Data Transformations
@info "Classifying ltmp and shelf regions."

centroids = AG.centroid.(dom.site_data.geom)

function contained_in(polys,centroid)::Bool
    return any([AG.contains(p, centroid) for p in polys])
end
inner_loc = [contained_in(inner_shelf.geom, cent) for cent in centroids]
outer_loc = [contained_in(outer_shelf.geom, cent) for cent in centroids]

features[!,         :shelf_position] .= 2
features[outer_loc, :shelf_position] .= 3
features[inner_loc, :shelf_position] .= 1

shape_file = "C:/users/dtan/data/ltmp_shape/gbr_3Zone 2.shp"

region_shps = GDF.read(shape_file)
region_shps[!, :convex_hull] = AG.convexhull.(region_shps.geometry)
reload_shp = false

centroids = [AG.centroid(polygn) for polygn in dom.site_data.geom]

north_mask = BitVector([AG.contains(region_shps.convex_hull[1], cent) || AG.gety(cent, 0) > -11.5 for cent in centroids])
central_mask = BitVector([AG.contains(region_shps.convex_hull[2], cent) for cent in centroids])
south_mask = BitVector([AG.contains(region_shps.convex_hull[3], cent) || AG.gety(cent, 0) < -21.0 for cent in centroids])
not_contained = (!).(north_mask .|| central_mask .|| south_mask)

features[!,            :ltmp_region] .= 0
features[north_mask,   :ltmp_region] .= 1
features[central_mask, :ltmp_region] .= 2
features[south_mask,   :ltmp_region] .= 3

@info "Writing updated location data to $(OUTPUT_CSV)"
CSV.write(OUTPUT_CSV, features)
