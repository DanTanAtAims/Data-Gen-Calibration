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

using Statistics

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

all_taxa_windows = construct_windows(growth_data, :logdiam, widths, n_windows)

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

    growth_spat_masks         = [growth_data[:, lvl_name] .== clst for clst in loc_nms]
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
        local f = plot_locs_survival(df, t_win, msks, loc_nms, "$(nm) Survival")
        save(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)", "$(nm)_survival_$(lvls).png"), f)
        local f = plot_implied_coefficient_survival(df, t_win, msks, loc_nms, "$(nm) survival coefficient")
        save(joinpath(OUT_PLOT_DIR, "mortality", "survival_$(lvls)", "$(nm)_survival_coef_$(lvls).png"), f)
    end
end
