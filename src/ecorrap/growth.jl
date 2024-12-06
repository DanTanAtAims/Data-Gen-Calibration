"""
Calculate the variation in linear extension between sites and/or regions.
"""

include("common.jl")
include("species.jl")
include("preprocessing.jl")
include("utils.jl")

using CSV, DataFrames

using CairoMakie

using Statistics

# Get data relating to growth statistics
growth_data = get_growth_entries(CSV.read(CORAL_OBS_PATH, DataFrame))

# Functional group specific data
tabular_acropora_data       = growth_data[is_tabular_Acropora.(growth_data.Taxa), :]
corymbose_acropora_data     = growth_data[is_corymbose_Acropora.(growth_data.Taxa), :]
corymbose_non_acropora_data = growth_data[is_corymbose_non_Acropora.(growth_data.Taxa), :]
small_massive_data          = growth_data[is_small_massive.(growth_data.Taxa), :]
large_massive_data          = growth_data[is_large_massive.(growth_data.Taxa), :]

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

unique_sites = unique(growth_data.Site_UID)

# Plot the distribution of all data points across all sites
site_counts = [
    count(growth_data.Site_UID .== site_uid) for site_uid in unique_sites
]

f = plot_observation_over_sites(site_counts, "", "site", unique_sites)
save(joinpath(OUT_PLOT_DIR, "growth", "site_count_dists", "ecorrap_growth_site_barplot.png"), f)

# Plot the distribution of data points across all sites limited to functional group
for (nm, df) in zip(taxa_nms, taxa_dfs)
    local site_counts = [
        count(df.Site_UID .== site_uid) for site_uid in unique_sites
    ]

    local f = plot_observation_over_sites(site_counts, nm, "site", unique_sites)
    save(joinpath(OUT_PLOT_DIR, "growth", "site_count_dists", "ecorrap_growth_site_$(nm).png"), f)
end

# Plot the distribution of all data points across all clusters (offshore far northern, etc.)
unique_clusters = unique(growth_data.Cluster)

cluster_counts = [
    count(growth_data.Cluster .== cluster_name) for cluster_name in unique_clusters
]

# Plot the distribution of data points across all clusters limited to functional group
f = plot_observation_over_sites(cluster_counts, "", "Cluster", unique_clusters)
save(joinpath(OUT_PLOT_DIR, "growth", "cluster_count_dists", "ecorrap_growth_cluster_barplot.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local cluster_counts = [
        count(df.Cluster .== cluster_name) for cluster_name in unique_clusters
    ]

    local f = plot_observation_over_sites(cluster_counts, nm, "Cluster", unique_clusters)
    save(joinpath(OUT_PLOT_DIR, "growth", "cluster_count_dists", "ecorrap_growth_cluster_$(nm).png"), f)
end

# Plot the distribution of all data points across all reefs
unique_reefs = unique(growth_data.Reef)

reef_counts = [
    count(growth_data.Reef .== reef_name) for reef_name in unique_reefs
]

# Plot the distribution of data points across all clusters limited to functional group
f = plot_observation_over_sites(reef_counts, "", "Reef", unique_reefs)
save(joinpath(OUT_PLOT_DIR, "growth", "reef_count_dists", "ecorrap_growth_reef_barplot.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local reef_counts = [
        count(df.Reef .== reef_name) for reef_name in unique_reefs
    ]

    local f = plot_observation_over_sites(reef_counts, nm, "Reef", unique_reefs)
    save(joinpath(OUT_PLOT_DIR, "growth", "reef_count_dists", "ecorrap_growth_reef_$(nm).png"), f)
end

# Plot log size distribution

f = plot_size_distribution(growth_data.diam, "All Taxa")
save(joinpath(OUT_PLOT_DIR, "growth", "size_dists", "all_taxa_size_dist.png"), f)

f = plot_size_distribution(growth_data.logdiam, "All Taxa"; is_log=true)
save(joinpath(OUT_PLOT_DIR, "growth", "log_size_dists", "all_taxa_size_dist_log.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local f = plot_size_distribution(df.diam, nm)
    save(joinpath(OUT_PLOT_DIR, "growth", "size_dists", "$(nm)_size_dist.png"), f)
    local f = plot_size_distribution(df.logdiam, nm; is_log=true)
    save(joinpath(OUT_PLOT_DIR, "growth", "log_size_dists", "$(nm)_size_dist_log.png"), f)
end

# construct moving mean linear extension values for all locations

# Take a moving average of linear extension rates
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


f = plot_linear_extension(growth_data, all_taxa_windows, "All Taxa Growth Rate")
save(joinpath(OUT_PLOT_DIR, "growth", "linear_extension", "All_Taxa_Llnear_Extension.png"), f)


for (nm, df, t_win) in zip(taxa_nms, taxa_dfs, taxa_wins)
    local f = plot_linear_extension(df, t_win, "$(nm) linear extension")
    save(joinpath(OUT_PLOT_DIR, "growth", "linear_extension", "$(nm)_lin_ext.png"), f)
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

    taxa_masks = [
        tab_acro_spat_masks,
        corym_acro_spat_masks,
        corym_non_acro_spat_masks,
        small_massive_spat_masks,
        large_massive_spat_masks
    ]

    for (nm, df, t_win, msks) in zip(taxa_nms, taxa_dfs, taxa_wins, taxa_masks)
        local f = plot_locs_linear_extension(df, t_win, msks, loc_nms, "$(nm) linear extension")
        save(joinpath(OUT_PLOT_DIR, "growth", "linear_extension_$(lvls)", "$(nm)_lin_ext_$(lvls).png"), f)
        local f = plot_implied_coefficient(df, t_win, msks, loc_nms, "$(nm) linear extension coefficient")
        save(joinpath(OUT_PLOT_DIR, "growth", "linear_extension_$(lvls)", "$(nm)_lin_ext_coef_$(lvls).png"), f)
    end

    f = plot_mean_coefficients(taxa_dfs, taxa_wins, taxa_masks, loc_nms, taxa_nms)
    save(joinpath(OUT_PLOT_DIR, "growth", "linear_extension_$(lvls)", "taxa_lin_ext_coef_$(lvls).png"), f)
end
