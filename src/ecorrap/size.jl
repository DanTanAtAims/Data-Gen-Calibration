"""
"""

include("common.jl")
include("species.jl")
include("preprocessing.jl")
include("utils.jl")

using CSV, DataFrames

using CairoMakie

using Distributions, Statistics
using ProgressMeter

all_fg_data = get_survival_entries(CSV.read(CORAL_OBS_PATH, DataFrame))
all_fg_data[!, :Cluster]  .= String.(all_fg_data.Cluster)
all_fg_data[!, :Reef]     .= String.(all_fg_data.Reef)
all_fg_data[!, :Site_UID] .= String.(all_fg_data.Site_UID)

tabular_acropora_data       = all_fg_data[is_tabular_Acropora.(all_fg_data.Taxa), :]
corymbose_acropora_data     = all_fg_data[is_corymbose_Acropora.(all_fg_data.Taxa), :]
corymbose_non_acropora_data = all_fg_data[is_corymbose_non_Acropora.(all_fg_data.Taxa), :]
small_massive_data          = all_fg_data[is_small_massive.(all_fg_data.Taxa), :]
large_massive_data          = all_fg_data[is_large_massive.(all_fg_data.Taxa), :]

taxa_nms = [
    "Tabular Acropora",
    "Coymbose Acropora",
    "Corymbose non-Acropora",
    "Small Massive",
    "Large Massive"
]

taxa_dfs = [
    tabular_acropora_data,
    corymbose_acropora_data,
    corymbose_non_acropora_data,
    small_massive_data,
    large_massive_data
]

unique_clusters = unique(all_fg_data.Cluster)
unique_reefs    = unique(all_fg_data.Reef)
unique_sites    = unique(all_fg_data.Site_UID)

OUT_SIZE_DIR = joinpath(@__DIR__, "..", "..", "Outputs", "Plots", "size")

function plot_size_distribution(
    df::DataFrame,
    location::String,
    taxa_name::String;
    unit::String="cm",
    fld::Symbol=:diam
)::Figure
    f = Figure(; size=(1200, 900))
    Axis(
        f[1, 1],
        xlabel=unit,
        ylabel="Colony Count",
        title="$(location): Distribution of $(taxa_name) colony size"
    )
    hist!(
        df[:, fld],
        bins=100,
        color=(:blue, 0.3),
        strokecolor=:black,
        strokewidth=1
    )
    return f
end
function plot_log_size_distribution(
    df::DataFrame,
    location::String,
    taxa_name::String
)::Figure
    return plot_size_distribution(df, location, taxa_name; unit="log cm", fld=:logdiam)
end

for tn in taxa_nms
    mkpath(joinpath(OUT_SIZE_DIR, "clusters", tn))
    mkpath(joinpath(OUT_SIZE_DIR, "reefs", tn))
    mkpath(joinpath(OUT_SIZE_DIR, "sites", tn))
end

@showprogress desc="Cluster Size Dists" for cluster in unique_clusters
    for (taxa_name, taxa_df) in zip(taxa_nms, taxa_dfs)
        msk = taxa_df.Cluster .== cluster
        count(msk) .== 0 && continue
        f = plot_size_distribution(taxa_df[msk, :], cluster, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "clusters", taxa_name, "size_dist_$(cluster).png"), f)
        f = plot_log_size_distribution(taxa_df[msk, :], cluster, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "clusters", taxa_name, "log_size_dist_$(cluster).png"), f)
    end
end
@showprogress desc="Reef Size Dists" for reef in unique_reefs
    for (taxa_name, taxa_df) in zip(taxa_nms, taxa_dfs)
        msk = taxa_df.Reef .== reef
        count(msk) .== 0 && continue
        f = plot_size_distribution(taxa_df[msk, :], reef, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "reefs", taxa_name, "size_dist_$(reef).png"), f)
        f = plot_log_size_distribution(taxa_df[msk, :], reef, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "reefs", taxa_name, "log_size_dist_$(reef).png"), f)
    end
end
@showprogress desc="Size Size Dists" for site in unique_sites
    for (taxa_name, taxa_df) in zip(taxa_nms, taxa_dfs)
        msk = taxa_df.Site_UID .== site
        count(msk) .== 0 && continue
        f = plot_size_distribution(taxa_df[msk, :], site, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "sites", taxa_name, "size_dist_$(site).png"), f)
        f = plot_log_size_distribution(taxa_df[msk, :], site, taxa_name)
        save(joinpath(OUT_SIZE_DIR, "sites", taxa_name, "log_size_dist_$(site).png"), f)
    end
end

reef_means = zeros(Float64, 5, length(unique_reefs))

f = Figure()
ax = Axis(f[1, 1], xlabel="Functional Group", ylabel="Mean Logdiam", title="Mean LogDiam")

@showprogress desc="Calculating size mean" for (r_idx, reef) in enumerate(unique_reefs)
    for (idx, taxa_df) in enumerate(taxa_dfs)
        reef_msk = taxa_df.Reef .== reef
        reef_means[idx, r_idx] = mean(taxa_df[reef_msk, :logdiam])
        scatter!(ax, 1:5, reef_means[:, r_idx])
    end
end

