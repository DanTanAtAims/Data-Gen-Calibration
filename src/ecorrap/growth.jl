"""
Calculate the variation in linear extension between sites and/or regions.
"""

include("common.jl")
include("species.jl")
include("preprocessing.jl")

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

"""
    plot_observation_over_sites(counts::Vector{Int64}, taxa_name::String)::Figure
"""
function plot_observation_over_sites(counts::Vector{Int64}, taxa_name::String)::Figure

    f = Figure(; size=(1200, 900))
    Axis(
        f[1, 1],
        xlabel="Observation Site",
        ylabel="Observation Counts",
        title="EcoRRAP $(taxa_name) Observation site distributions"
    )
    barplot!(
        1:length(counts),
        counts,
        color=(:blue, 0.3),
        strokecolor=:black,
        strokewidth=1
    )

    return f
end

# Plot the distribution of all data points across all sites
site_counts = [
    count(growth_data.Site_UID .== site_uid) for site_uid in unique_sites
]

f = plot_observation_over_sites(site_counts, "")
save(joinpath(OUT_PLOT_DIR, "site_count_dists", "ecorrap_growth_site_barplot.png"), f)

# Plot the distribution of data points across all sites limited to functional group
for (nm, df) in zip(taxa_nms, taxa_dfs)
    local site_counts = [
        count(df.Site_UID .== site_uid) for site_uid in unique_sites
    ]

    local f = plot_observation_over_sites(site_counts, nm)
    save(joinpath(OUT_PLOT_DIR, "site_count_dists", "ecorrap_growth_site_$(nm).png"), f)
end

# Plot log size distribution

"""
    plot_size_distribution(sizes::Vector{Float64}, taxa_name::String)::Figure
"""
function plot_size_distribution(sizes::Vector{Float64}, taxa_name::String; is_log=false)::Figure
    f = Figure(; size=(1200, 900))
    Axis(
        f[1, 1],
        xlabel=is_log ? "log cm" : "cm",
        ylabel="Colony Count",
        title="Distribution of $(taxa_name) colony size"
    )
    hist!(
        sizes,
        bins=100,
        color=(:blue, 0.3),
        strokecolor=:black,
        strokewidth=1
    )
    return f
end

f = plot_size_distribution(growth_data.diam, "All Taxa")
save(joinpath(OUT_PLOT_DIR, "size_dists", "all_taxa_size_dist.png"), f)

f = plot_size_distribution(growth_data.logdiam, "All Taxa"; is_log=true)
save(joinpath(OUT_PLOT_DIR, "log_size_dists", "all_taxa_size_dist_log.png"), f)

for (nm, df) in zip(taxa_nms, taxa_dfs)
    local f = plot_size_distribution(df.diam, nm)
    save(joinpath(OUT_PLOT_DIR, "size_dists", "$(nm)_size_dist.png"), f)
    local f = plot_size_distribution(df.logdiam, nm; is_log=true)
    save(joinpath(OUT_PLOT_DIR, "log_size_dists", "$(nm)_size_dist_log.png"), f)
end

# construct moving mean linear extension values for all locations

# convenient type aliases
WIN_TYPE = Tuple{Vector{Float64}, Vector{Float64}}
MFLOAT = Union{Float64, Missing} # Maybe Float
"""
    construct_windows(min_val::Float64, max_val::Float64, width::Float64, n_windows::Int64)::Tuple{Vector{Float64}, Vector{Float64}}
    construct_windows(df::DataFrame, field::Symbol, width::Float64, n_windows::Int64)

Construct window lower and upper bounds for moving averages.
"""
function construct_windows(
    df::DataFrame, field::Symbol, width::Float64, n_windows::Int64
)
    return construct_windows(
        minimum(df[:, field]),
        maximum(df[:, field]),
        width,
        n_windows
    )
end
function construct_windows(
    min_val::Float64,
    max_val::Float64,
    width::Float64,
    n_windows::Int64
)::Tuple{Vector{Float64}, Vector{Float64}}
    if max_val - min_val <= width
        msg  = "Given width is grater than the difference "
        msg *= "between maximum and minimum"
        throw(ArgumentError(msg))
    end

    highest_lb = max_val - width
    lower_bounds::Vector{Float64} = collect(
        min_val:((highest_lb - min_val) / n_windows):highest_lb
    )[1:n_windows]

    lowest_ub = min_val + width
    upper_bounds::Vector{Float64} = collect(
        lowest_ub:((max_val - lowest_ub) / n_windows):max_val
    )[1:n_windows]

    return lower_bounds, upper_bounds
end

n_windows = 200
widths = 0.4

all_taxa_windows = construct_windows(growth_data, :logdiam, widths, n_windows)

tabular_acropora_win       = construct_windows(tabular_acropora_data, :logdiam, widths, n_windows)
corymbose_acropora_win     = construct_windows(corymbose_acropora_data, :logdiam, widths, n_windows)
corymbose_non_acropora_win = construct_windows(corymbose_non_acropora_data, :logdiam, widths, n_windows)
small_massive_win          = construct_windows(small_massive_data, :logdiam, widths, n_windows)
large_massive_win          = construct_windows(large_massive_data, :logdiam, widths, n_windows)

taxa_wins = [
    tabular_acropora_win,
    corymbose_acropora_win,
    corymbose_non_acropora_win,
    small_massive_win,
    large_massive_win
]

function plot_linear_extension(
    df::DataFrame,
    windows::WIN_TYPE,
    title::String;
    count_threshold::Int64=20,
)::Figure
    linexts_means::Vector{MFLOAT}  = Vector{MFLOAT}(missing, length(windows[1]))
    linexts_stdevs::Vector{MFLOAT} = Vector{MFLOAT}(missing, length(windows[1]))

    for (idx, win) in enumerate(zip(windows...))
        win_mask = win[1] .< df.logdiam .< win[2]
        if count(win_mask) .< count_threshold
            continue
        end
        linexts_means[idx]  = mean(df[win_mask, :lin_ext])
        linexts_stdevs[idx] = std(df[win_mask, :lin_ext])
    end

    xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    non_missing_mask = (!).(ismissing.(linexts_means))
    xs             = xs[non_missing_mask]
    linexts_means  = linexts_means[non_missing_mask]
    linexts_stdevs = linexts_stdevs[non_missing_mask]

    f = Figure(; size=(1200, 900))
    Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Linear Extension (cm)", title=title)
    lines!(xs, linexts_means; color=:black)
    band!(
        xs, linexts_means .- linexts_stdevs, linexts_means .+ linexts_stdevs;
        color=(:blue, 0.3)
    )

    return f
end

f = plot_linear_extension(growth_data, all_taxa_windows, "All Taxa Growth Rate")
save(joinpath(OUT_PLOT_DIR, "linear_extension", "All_Taxa_Llnear_Extension.png"), f)


for (nm, df, t_win) in zip(taxa_nms, taxa_dfs, taxa_wins)
    local f = plot_linear_extension(df, t_win, "$(nm) linear extension")
    save(joinpath(OUT_PLOT_DIR, "linear_extension", "$(nm)_lin_ext.png"), f)
end
