"""
Calculate and plot statistics relating to the mortality of rates of corals observed in the
EcoRRAP dataset.
"""

include("common.jl")
include("species.jl")
include("preprocessing.jl")
include("utils.jl")

using CSV, DataFrames

using CairoMakie

using Distributions, Statistics

# Get data relating to growth statistics
mortality_data = get_survival_entries(CSV.read(CORAL_OBS_PATH, DataFrame))

# Functional group specific data
tabular_acropora_data       = mortality_data[is_tabular_Acropora.(mortality_data.Taxa), :]
corymbose_acropora_data     = mortality_data[is_corymbose_Acropora.(mortality_data.Taxa), :]
corymbose_non_acropora_data = mortality_data[is_corymbose_non_Acropora.(mortality_data.Taxa), :]
small_massive_data          = mortality_data[is_small_massive.(mortality_data.Taxa), :]
large_massive_data          = mortality_data[is_large_massive.(mortality_data.Taxa), :]

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

unique_sites = unique(mortality_data.Site_UID)

# Plot the distribution of all data points across all sites
site_counts = [
    count(mortality_data.Site_UID .== site_uid) for site_uid in unique_sites
]

mkpath(joinpath(OUT_PLOT_DIR, "mortality", "site_count_dists"))

f = plot_observation_over_sites(site_counts, "", "site", unique_sites)
save(joinpath(OUT_PLOT_DIR, "mortality", "site_count_dists", "ecorrap_growth_site_barplot.png"), f)

# Plot the distribution of data points across all sites limited to functional group
for (nm, df) in zip(taxa_nms, taxa_dfs)
    local site_counts = [
        count(df.Site_UID .== site_uid) for site_uid in unique_sites
    ]

    local f = plot_observation_over_sites(site_counts, nm, "site", unique_sites)
    save(joinpath(OUT_PLOT_DIR, "mortality", "site_count_dists", "ecorrap_growth_site_$(nm).png"), f)
end

# Plot the distribution of all data points across all clusters (offshore far northern, etc.)
unique_clusters = unique(mortality_data.Cluster)

cluster_counts = [
    count(mortality_data.Cluster .== cluster_name) for cluster_name in unique_clusters
]

mkpath(joinpath(OUT_PLOT_DIR, "mortality", "cluster_count_dists"))

# Plot the distribution of data points across all clusters limited to functional group
f = plot_observation_over_sites(cluster_counts, "", "Cluster", unique_clusters)
save(joinpath(OUT_PLOT_DIR, "mortality", "cluster_count_dists", "ecorrap_growth_cluster_barplot.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local cluster_counts = [
        count(df.Cluster .== cluster_name) for cluster_name in unique_clusters
    ]

    local f = plot_observation_over_sites(cluster_counts, nm, "Cluster", unique_clusters)
    save(joinpath(OUT_PLOT_DIR, "mortality", "cluster_count_dists", "ecorrap_growth_cluster_$(nm).png"), f)
end

# Plot the distribution of all data points across all reefs
unique_reefs = unique(mortality_data.Reef)

reef_counts = [
    count(mortality_data.Reef .== reef_name) for reef_name in unique_reefs
]

mkpath(joinpath(OUT_PLOT_DIR, "mortality", "reef_count_dists"))

# Plot the distribution of data points across all clusters limited to functional group
f = plot_observation_over_sites(reef_counts, "", "Reef", unique_reefs)
save(joinpath(OUT_PLOT_DIR, "mortality", "reef_count_dists", "ecorrap_growth_reef_barplot.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local reef_counts = [
        count(df.Reef .== reef_name) for reef_name in unique_reefs
    ]

    local f = plot_observation_over_sites(reef_counts, nm, "Reef", unique_reefs)
    save(joinpath(OUT_PLOT_DIR, "mortality", "reef_count_dists", "ecorrap_growth_reef_$(nm).png"), f)
end

# Plot log size distribution
mkpath(joinpath(OUT_PLOT_DIR, "mortality", "log_size_dists"))
mkpath(joinpath(OUT_PLOT_DIR, "mortality", "size_dists"))

f = plot_size_distribution(mortality_data.diam, "All Taxa")
save(joinpath(OUT_PLOT_DIR, "mortality", "size_dists", "all_taxa_size_dist.png"), f)

f = plot_size_distribution(mortality_data.logdiam, "All Taxa"; is_log=true)
save(joinpath(OUT_PLOT_DIR, "mortality", "log_size_dists", "all_taxa_size_dist_log.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local f = plot_size_distribution(df.diam, nm)
    save(joinpath(OUT_PLOT_DIR, "mortality", "size_dists", "$(nm)_size_dist.png"), f)
    local f = plot_size_distribution(df.logdiam, nm; is_log=true)
    save(joinpath(OUT_PLOT_DIR, "mortality", "log_size_dists", "$(nm)_size_dist_log.png"), f)
end

n_windows = 200
widths = 0.4

all_taxa_windows = construct_windows(mortality_data, :logdiam, widths, n_windows)

tabular_acropora_win       = construct_windows(tabular_acropora_data,       :logdiam, widths, n_windows)
corymbose_acropora_win     = construct_windows(corymbose_acropora_data,     :logdiam, widths, n_windows)
corymbose_non_acropora_win = construct_windows(corymbose_non_acropora_data, :logdiam, widths, n_windows)
small_massive_win          = construct_windows(small_massive_data,          :logdiam, widths, n_windows)
large_massive_win          = construct_windows(large_massive_data,          :logdiam, widths, n_windows)

taxa_wins = [
    tabular_acropora_win,
    corymbose_acropora_win,
    corymbose_non_acropora_win,
    small_massive_win,
    large_massive_win
]

mkpath(joinpath(OUT_PLOT_DIR, "mortality", "survival"))

f = plot_survival(mortality_data, all_taxa_windows, "All Taxa Survival Probability")
save(joinpath(OUT_PLOT_DIR, "mortality", "survival", "All_Taxa_survival.png"), f)


for (nm, df, t_win) in zip(taxa_nms, taxa_dfs, taxa_wins)
    local f = plot_survival(df, t_win, "$(nm) Survival Probability")
    save(joinpath(OUT_PLOT_DIR, "mortality", "survival", "$(nm)_survival.png"), f)
end

spatial_levels = ["clusters", "reefs", "sites"]
spatial_cols = [:Cluster, :Reef, :Site_UID]
spatial_unique = [unique_clusters, unique_reefs, unique_sites]

for (lvls, lvl_name, loc_nms) in zip(spatial_levels, spatial_cols, spatial_unique)

    growth_spat_masks         = [mortality_data[:, lvl_name] .== clst for clst in loc_nms]
    tab_acro_spat_masks       = [tabular_acropora_data[:, lvl_name] .== clst for clst in loc_nms]
    corym_acro_spat_masks     = [corymbose_acropora_data[:, lvl_name] .== clst for clst in loc_nms]
    corym_non_acro_spat_masks = [corymbose_non_acropora_data[:, lvl_name] .== clst for clst in loc_nms]
    small_massive_spat_masks  = [small_massive_data[:, lvl_name] .== clst for clst in loc_nms]
    large_massive_spat_masks  = [large_massive_data[:, lvl_name] .== clst for clst in loc_nms]

    local taxa_masks = [
        tab_acro_spat_masks,
        corym_acro_spat_masks,
        corym_non_acro_spat_masks,
        small_massive_spat_masks,
        large_massive_spat_masks
    ]

    mkpath(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)"))

    for (nm, df, t_win, msks) in zip(taxa_nms, taxa_dfs, taxa_wins, taxa_masks)
        local f = plot_locs_survival(df, t_win, msks, loc_nms, "$(nm) Survival"; count_threshold=20)
        save(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)", "$(nm)_survival_$(lvls).png"), f)
        local f = plot_implied_coefficient_survival(df, t_win, msks, loc_nms, "$(nm) survival coefficient"; count_threshold=20)
        save(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)", "$(nm)_survival_coef_$(lvls).png"), f)
        local f = plot_implied_quantile_survival(df, t_win, msks, loc_nms, "$(nm) survival quantile"; count_threshold=20)
        save(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)", "$(nm)_survival_quantile_$(lvls).png"), f)
    end
end

function apply_scale_coefficient(ys::Vector{Float64}, coef::Float64)::Vector{Float64}
    if coef >= 0.0
        return ys .+ (1 .- ys) .* coef
    else
        return ys .* (1 + coef)
    end
end

function plot_scale_impact()::Figure
    surv = x -> -1 / (x + 2) + 1
    xs = 0.0:0.1:30.0
    ys = surv.(xs)
    f = Figure()
    Axis(
        f[1, 1],
        xlabel="Size",
        ylabel="Probability of Survival",
        title="Application of Survival Coefficients"
    )
    lines!(xs, ys, color=:black)

    for coef in -0.3:0.01:0.0
        lines!(xs, apply_scale_coefficient(ys, coef), color=(:red, 0.3))
    end
    for coef in 0.0:0.05:0.5
        lines!(xs, apply_scale_coefficient(ys, coef), color=(:blue, 0.3))
    end
    return f
end

function plot_constant_scale_impact()::Figure
    surv = x -> -1 / (x + 2) + 1
    xs = 0.0:0.1:30.0
    ys = surv.(xs)
    f = Figure()
    Axis(
        f[1, 1],
        xlabel="Size",
        ylabel="Probability of Survival",
        title="Application of Survival Coefficients"
    )
    lines!(xs, ys, color=:black)

    for coef in 0.8:0.01:1.0
        lines!(xs, coef .* ys, color=(:red, 0.3))
    end
    for coef in 1.0:0.025:1.5
        lines!(xs, coef .* ys, color=(:blue, 0.3))
    end
    return f
end

function beta_quantile_prop(coef::Float64, prop::Float64)::Float64
    alpha::Float64 = 2 / (1 - prop)
    beta::Float64 = 2 / prop
    dist = Beta(alpha, beta)
    return quantile(dist, coef)
end

function plot_beta_scale_impact()::Figure
    surv = x -> -1 / (x + 2) + 1.0
    xs = 0.0:0.1:30.0
    ys = surv.(xs)
    f = Figure()
    Axis(
        f[1, 1],
        xlabel="Size",
        ylabel="Probability of Survival",
        title="Application of Survival Coefficients"
    )
    lines!(xs, ys, color=:black)

    for coef in 0.1:0.05:0.5
        lines!(xs, beta_quantile_prop.(coef, ys), color=(:red, 0.4))
    end
    for coef in 0.5:0.05:0.90
        lines!(xs, beta_quantile_prop.(coef, ys), color=(:blue, 0.4))
    end
    return f
end

function calculate_adjusted(prop::Float64, val::Float64, proportion_max::Float64)::Float64
    width::Float64 = (1 - prop) * proportion_max
    return prop + width * val
end

function plot_top_bounded_impact(proportion_max::Float64)::Figure
    surv = x -> -1 / (x + 2) + 1.0
    xs = 0.0:0.1:30.0
    ys = surv.(xs)
    f = Figure()
    Axis(
        f[1, 1],
        xlabel="Size",
        ylabel="Probability of Survival",
        title="Application of Survival Coefficients"
    )
    lines!(xs, ys, color=:black)
    for coef in 0.0:0.05:1.0
        lines!(xs, calculate_adjusted.(ys, coef, proportion_max), color=(:red, 0.4))
    end
    for coef in -1.0:0.05:0.0
        lines!(xs, calculate_adjusted.(ys, coef, proportion_max), color=(:blue, 0.4))
    end
    return f
end
