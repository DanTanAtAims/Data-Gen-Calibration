# Output Data

All data is written to the `Outputs` directory.

The repository generates data in a usable format for the calibration of ADRIA/CoralBlox.

## Location Classification Data

Location from the canonical geopackage were grouped by LTMP region (north, central, south),
shelf position (inner, mid, outer) and split by 33%, 66% quantiles of depth, turbidity and
wave activity.

The file generated is written to `Outputs/location_classification.csv`

The classification is listed in two columns, `classification` and `consecutive_
classification`. `classification` number can be as large as 243 but not all classes have
locations. `consecutive_classification` removes classes with no location.

## Long Term Monitoring Program Data

The following files are written to the Outputs directory.

- `ltmp_reefmod_manta.gpkg`

Manta tow data aligned with reefmod polygons and `RME_UNIQUE_ID`. There may be multiple manta
tow sites to a single `RME_UNIQUE_ID`, and LTMP sites with no `RME_UNIQUE_ID`. The
geopackage also contains the `GBRMPA ID` and `REEF ID`.

- `ltmp_reefmod_photo.gpkg`

LTMP fixed site photogrammetry data with added `GBRMPA_ID` and `RME_UNIQUE_ID`.
Furthermore, data points occuring and the same canonical reefs location were aggregated
using the mean.

- `ltmp_reefmod_photo_hc.gpkg`

LTMP fixed site photogrammetry data with added `GBRMPA_ID` and `RME_UNIQUE_ID`.
Furthermore, data points occuring and the same canonical reefs location were aggregated
using the mean. This geopackage only contains hard coral data is has columns for each year
observations were collected.

The data frame column names are are

```julia-repl
geometry | REEF_ID | LATITUDE_mean | LONGITUDE_mean | GBRMPA_ID | RME_UNIQUE_ID | 1992 | 1994 | 1995 ...         1996
```

- `ltmp_reefmod_photo_hc.gpkg`

LTMP fixed site photogrammetry data with added `GBRMPA_ID` and `RME_UNIQUE_ID`.
Furthermore, data points occuring and the same canonical reefs location were aggregated
using the mean. This geopackage only contains the sum of hard and soft coral covere is has columns for each year
observations were collected.

## Classification Level Observation Data

*To do*
