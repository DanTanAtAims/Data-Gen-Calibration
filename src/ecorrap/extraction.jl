"""
Extract mean and standard deviation of growth and mortality values and save them to NetCDF.
"""

include("common.jl")
include("preprocessing.jl")
include("species.jl")

using CSV, DataFrames, YAXArrays, NetCDF

using CairoMakie

using Statistics

growth_data = get_growth_entries(CSV.read(CORAL_OBS_PATH, DataFrame))

# Functional group specific data
tabular_acropora_data       = growth_data[is_tabular_Acropora.(growth_data.Taxa), :]
corymbose_acropora_data     = growth_data[is_corymbose_Acropora.(growth_data.Taxa), :]
corymbose_non_acropora_data = growth_data[is_corymbose_non_Acropora.(growth_data.Taxa), :]
small_massive_data          = growth_data[is_small_massive.(growth_data.Taxa), :]
large_massive_data          = growth_data[is_large_massive.(growth_data.Taxa), :]

# Construct a list to make iteration easier
taxa_dfs = [
    tabular_acropora_data,
    corymbose_acropora_data,
    corymbose_non_acropora_data,
    small_massive_data,
    large_massive_data
]

"""
    bin_edges()

Different bin edges used in calibration processes to allow for different growth rates.
"""
function bin_edges()
    return Matrix(
        [
            2.5 7.5 12.5 25.0 50.0 80.0 120.0 160.0;
            2.5 7.5 12.5 20.0 30.0 60.0 100.0 150.0;
            2.5 7.5 12.5 20.0 30.0 40.0 50.0 60.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    )
end

dims = (
    Dim{:taxa}(1:5),
    Dim{:size_class}(1:7)
)

# Create arrays for growth statistics
mean_lin_ext = zeros(Float64, 5, 7)
std_lin_ext = zeros(Float64, 5, 7)
stdmn_lin_ext = zeros(Float64, 5, 7)
lin_ext_counts = zeros(Int64, 5, 7)

bins = bin_edges()

# Iterate over taxa and size classes
for (idx, df) in enumerate(taxa_dfs)
    for s_idx in 1:7
        sc_fg_mask = (df.diam .>= bins[idx, s_idx]) .&& (df.diam .< bins[idx, s_idx+1])
        tmp_df = df[sc_fg_mask, :]
        mean_lin_ext[idx, s_idx]   = mean(tmp_df.lin_ext)
        std_lin_ext[idx, s_idx]    = std(tmp_df.lin_ext)
        lin_ext_counts[idx, s_idx] = count(sc_fg_mask)
        stdmn_lin_ext[idx, s_idx]  = std(tmp_df.lin_ext) ./ (count(sc_fg_mask) - 1)
   end
end

lin_ext_μ = YAXArray(dims, mean_lin_ext)
lin_ext_σ = YAXArray(dims, std_lin_ext)
lin_ext_n = YAXArray(dims, lin_ext_counts)

survival_data = get_survival_entries(CSV.read(CORAL_OBS_PATH, DataFrame))

# Functional group specific data
tabular_acropora_data       = survival_data[is_tabular_Acropora.(survival_data.Taxa), :]
corymbose_acropora_data     = survival_data[is_corymbose_Acropora.(survival_data.Taxa), :]
corymbose_non_acropora_data = survival_data[is_corymbose_non_Acropora.(survival_data.Taxa), :]
small_massive_data          = survival_data[is_small_massive.(survival_data.Taxa), :]
large_massive_data          = survival_data[is_large_massive.(survival_data.Taxa), :]

# Construct a list to make iteration easier
taxa_dfs = [
    tabular_acropora_data,
    corymbose_acropora_data,
    corymbose_non_acropora_data,
    small_massive_data,
    large_massive_data
]

surv_prob = zeros(Float64, 5, 7)
surv_count = zeros(Int64, 5, 7)

for (idx, df) in enumerate(taxa_dfs)
    for s_idx in 1:7
        sc_fg_mask = (df.diam .>= bins[idx, s_idx]) .&& (df.diam .< bins[idx, s_idx+1])
        tmp_df = df[sc_fg_mask, :]
        surv_prob[idx, s_idx]  = mean(tmp_df.surv)
        surv_count[idx, s_idx] = count(sc_fg_mask)
    end
end

mb_rate_p = YAXArray(dims, 1 .- surv_prob)
mb_rate_n = YAXArray(dims, surv_count)

# Fill missing values

# Small massive second size class onwards to match large massives
mb_rate_p[4, 2:end] .= mb_rate_p[5, 2:end]

# Tabular Acropora largest mb rate 5% less then second largest
mb_rate_p[1, 7] = mb_rate_p[1, end-1]

# Corymbose Acropora last two mb rates equal (less than 10 samples)
mb_rate_p[2, 6:7] .= mb_rate_p[2, 5]

# Corymbose non-Acropora last three rates equal (all less then 11 samples)
mb_rate_p[3, 5:7] .= mb_rate_p[3, 4]

# LINEAR EXTENSION

# Fill last Corymbose Acropora but my match tabular % increase
lin_ext_μ[2, 6] = lin_ext_μ[2, 5] .* (lin_ext_μ[1, 6] / lin_ext_μ[1, 5])
lin_ext_σ[2, 6] = lin_ext_σ[2, 5] .* (lin_ext_σ[1, 6] / lin_ext_σ[1, 5])

# Fill Corymbose non-Acopora with previous size class growth values
lin_ext_μ[3, 5:6] .= lin_ext_μ[3, 4]
lin_ext_σ[3, 5:6] .= lin_ext_σ[3, 4]

# Fill Small massive with large massives
lin_ext_μ[4, 2:6] .= lin_ext_μ[5, 2:6]
lin_ext_σ[4, 2:6] .= lin_ext_σ[5, 2:6]

filled_ds = Dataset(
    ;
    :filled_lin_ext_mean=>lin_ext_μ,
    :filled_lin_ext_stdev=>lin_ext_σ,
    :filled_mb_rate_mean=>mb_rate_p
)

savedataset(
    filled_ds;
    path=ECORRAP_SAVE_PATH,
    driver=:netcdf, overwrite=true
)
