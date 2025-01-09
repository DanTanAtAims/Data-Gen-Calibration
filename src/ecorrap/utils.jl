using CairoMakie

"""
    plot_observation_over_sites(counts::Vector{Int64}, taxa_name::String)::Figure
"""
function plot_observation_over_sites(counts::Vector{Int64}, taxa_name::String, loc_level::String, x_nms)::Figure

    f = Figure(; size=(1200, 900))
    Axis(
        f[1, 1],
        xticks=(1:length(x_nms), x_nms),
        xticksvisible=true,
        xlabel= loc_level,
        ylabel="Observation Counts",
        title="EcoRRAP $(taxa_name) Observation " * loc_level * " distributions"
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

# convenient type aliases
WIN_TYPE = Tuple{Vector{Float64}, Vector{Float64}}
MFLOAT = Union{Float64, Missing} # Maybe Float
LINEXT_STATS = Tuple{Vector{MFLOAT}, Vector{MFLOAT}} # mean and standard deviation

"""
    construct_windows(min_val::Float64, max_val::Float64, width::Float64, n_windows::Int64)::Tuple{Vector{Float64}, Vector{Float64}}
    construct_windows(df::DataFrame, field::Symbol, width::Float64, n_windows::Int64)

Construct window lower and upper bounds for moving averages.
"""
function construct_windows(
    df::DataFrame, field::Symbol, width::Float64, n_windows::Int64
)::WIN_TYPE
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
)::WIN_TYPE
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

function calculate_linear_extension(
    df::DataFrame,
    windows::WIN_TYPE;
    count_threshold::Int64=20
)::LINEXT_STATS
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
    return linexts_means, linexts_stdevs
end

function plot_linear_extension(
    df::DataFrame,
    windows::WIN_TYPE,
    title::String;
    count_threshold::Int64=20,
)::Figure

    linexts_means, linexts_stdevs = calculate_linear_extension(
        df, windows; count_threshold=count_threshold
    )

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

"""
    plot_locs_linear_extension(df::DataFrame, windows::WIN_TYPE, row_masks, title::String; count_threshold::Int64=20)::Figure

Plot the linear extension mean over all locations and plot the linear extension means for
each location as well.
"""
function plot_locs_linear_extension(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks,
    row_names,
    title::String;
    count_threshold::Int64=20
)::Figure

    linexts_means, linexts_stdevs = calculate_linear_extension(
        df, windows; count_threshold=count_threshold
    )

    _xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    non_missing_mask = (!).(ismissing.(linexts_means))
    xs             = _xs[non_missing_mask]
    linexts_means  = linexts_means[non_missing_mask]
    linexts_stdevs = linexts_stdevs[non_missing_mask]

    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Linear Extension (cm)", title=title)
    lines!(xs, linexts_means; color=:black)
    band!(
        xs, linexts_means .- linexts_stdevs, linexts_means .+ linexts_stdevs;
        color=(:blue, 0.3)
    )
    # Only display the legend if there was data to plot to avoid error
    plotted_a_location = false

    for (nm, r_msk) in zip(row_names, row_masks)
        linexts_means, linexts_stdevs = calculate_linear_extension(
            df[r_msk, :], windows; count_threshold=count_threshold
        )

        non_missing_mask = (!).(ismissing.(linexts_means))
        xs             = _xs[non_missing_mask]
        linexts_means  = linexts_means[non_missing_mask]
        if count(non_missing_mask) == 0
            continue
        end
        plotted_a_location = true

        scatter!(ax, xs, linexts_means, label=nm)
    end
    if plotted_a_location
        axislegend(ax)
    end
    return f
end

"""
    plot_implied_coefficient(df::DataFrame, windows::WIN_TYPE, row_masks, row_names, title::String; count_threshold::Int64=20)

Plot the implied scale coefficient over the log diameter.
"""
function plot_implied_coefficient(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks,
    row_names,
    title::String;
    count_threshold::Int64=20
)
    gbr_linexts_means, gbr_linexts_stdevs = calculate_linear_extension(
        df, windows; count_threshold=count_threshold
    )

    _xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Scale Coefficient", title=title, limits=(0.0, nothing, 0.0, nothing))
    lines!(_xs[[1, end]], [1.0, 1.0])
    plotted_a_location = false

    for (nm, r_msk) in zip(row_names, row_masks)
        linexts_means, linexts_stdevs = calculate_linear_extension(
            df[r_msk, :], windows; count_threshold=count_threshold
        )
        coefs = linexts_means ./ gbr_linexts_means

        non_missing_mask = (!).(ismissing.(coefs))
        xs             = _xs[non_missing_mask]
        coefs  = coefs[non_missing_mask]

        if count(non_missing_mask) == 0
            continue
        end
        plotted_a_location = true

        scatter!(ax, xs, coefs, label=nm)
    end
    if plotted_a_location
        axislegend(ax)
    end
    return f
end

function mean_coefficients(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks;
    count_threshold::Int64=20
)::Vector{MFLOAT}
    gbr_linexts_means, gbr_linexts_stdevs = calculate_linear_extension(
        df, windows; count_threshold=count_threshold
    )

    mean_coefficients::Vector{MFLOAT} = Vector{MFLOAT}(missing, length(row_masks))

    counts = 0
    for r_msk in row_masks
        counts += 1
        linexts_means, linexts_stdevs = calculate_linear_extension(
            df[r_msk, :], windows; count_threshold=count_threshold
        )
        coefs = linexts_means ./ gbr_linexts_means

        non_missing_mask = (!).(ismissing.(coefs))
        coefs  = coefs[non_missing_mask]

        if count(non_missing_mask) == 0
            continue
        end
        mean_coefficients[counts] = mean(coefs)
    end
    return mean_coefficients
end

function plot_mean_coefficients(
    dfs::Vector{DataFrame},
    wins::Vector{WIN_TYPE},
    row_maskss,
    row_nms,
    taxa_nms;
    count_threshold::Int64=20
)::Figure
    n_taxa::Int64    = length(dfs)
    n_spatial::Int64 = length(row_maskss[1])

    mean_coefs::Matrix{MFLOAT} = Matrix{MFLOAT}(missing, n_taxa, n_spatial)

    for (i, (df, win, msks)) in enumerate(zip(dfs, wins, row_maskss))
        mean_coefs[i, :] .= mean_coefficients(df, win, msks)
    end

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1],
        xticks=(1:n_taxa, taxa_nms),
        xticksvisible=true,
        xlabel="Taxa",
        ylabel="Coefficient",
        title="Mean Scaling Coefficient",
        limits=(nothing, nothing, 0, nothing)
    )

    for (idx, nm) in enumerate(row_nms)
        scatter!(ax, 1:n_taxa, mean_coefs[:, idx], label=nm)
    end
    axislegend(ax)
    return f
end

function calculate_survival(
    df::DataFrame,
    windows::WIN_TYPE;
    count_threshold::Int64=20
)::LINEXT_STATS
    linexts_means::Vector{MFLOAT}  = Vector{MFLOAT}(missing, length(windows[1]))
    linexts_stdevs::Vector{MFLOAT} = Vector{MFLOAT}(missing, length(windows[1]))

    for (idx, win) in enumerate(zip(windows...))
        win_mask = win[1] .< df.logdiam .< win[2]
        if count(win_mask) .< count_threshold
            continue
        end
        linexts_means[idx]  = mean(df[win_mask, :surv])
        linexts_stdevs[idx] = std(df[win_mask, :surv])
    end
    return linexts_means, linexts_stdevs
end

function plot_survival(
    df::DataFrame,
    windows::WIN_TYPE,
    title::String;
    count_threshold::Int64=20,
)::Figure

    surv_means, surv_stdevs = calculate_survival(
        df, windows; count_threshold=count_threshold
    )

    xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    non_missing_mask = (!).(ismissing.(surv_means))
    xs             = xs[non_missing_mask]
    surv_means  = surv_means[non_missing_mask]
    surv_stdevs = surv_stdevs[non_missing_mask]

    f = Figure(; size=(1200, 900))
    Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Survival Probability", title=title, limits=(0, nothing, 0, 1))
    lines!(xs, surv_means; color=:black)
    band!(
        xs, surv_means .- surv_stdevs, surv_means .+ surv_stdevs;
        color=(:blue, 0.3)
    )

    return f
end

function plot_locs_survival(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks,
    row_names,
    title::String;
    count_threshold::Int64=20
)::Figure

    surv_means, surv_stdevs = calculate_survival(
        df, windows; count_threshold=count_threshold
    )

    _xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    non_missing_mask = (!).(ismissing.(surv_means))
    xs             = _xs[non_missing_mask]
    surv_means  = surv_means[non_missing_mask]
    surv_stdevs = surv_stdevs[non_missing_mask]

    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Survival Probability", title=title, limits=(0, nothing, 0, 1))
    lines!(xs, surv_means; color=:black)
    band!(
        xs, surv_means .- surv_stdevs, surv_means .+ surv_stdevs;
        color=(:blue, 0.3)
    )
    # Only display the legend if there was data to plot to avoid error
    plotted_a_location = false

    for (nm, r_msk) in zip(row_names, row_masks)
        surv_means, surv_stdevs = calculate_survival(
            df[r_msk, :], windows; count_threshold=count_threshold
        )

        non_missing_mask = (!).(ismissing.(surv_means))
        xs             = _xs[non_missing_mask]
        surv_means  = surv_means[non_missing_mask]
        if count(non_missing_mask) == 0
            continue
        end
        plotted_a_location = true

        scatter!(ax, xs, surv_means, label=nm)
    end
    if plotted_a_location
        axislegend(ax; position=:rb)
    end
    return f
end

function plot_implied_coefficient_survival(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks,
    row_names,
    title::String;
    count_threshold::Int64=20
)
    gbr_surv_means, gbr_surv_stdevs = calculate_survival(
        df, windows; count_threshold=count_threshold
    )

    _xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    f = Figure(; size=(1200, 900))
    ax = Axis(f[1, 1], xlabel="Log Diameter (cm)", ylabel="Scale Coefficient", title=title, limits=(0.0, nothing, -1.0, 1.0))
    lines!(_xs[[1, end]], [0.0, 0.0])
    plotted_a_location = false

    for (nm, r_msk) in zip(row_names, row_masks)
        surv_means, surv_stdevs = calculate_survival(
            df[r_msk, :], windows; count_threshold=count_threshold
        )
        coefs = (surv_means .- gbr_surv_means)

        non_missing_mask = (!).(ismissing.(coefs))
        xs             = _xs[non_missing_mask]
        coefs  = coefs[non_missing_mask]
        coefs = coefs ./ (
            (1 .- gbr_surv_means[non_missing_mask]) .* (coefs .>= 0.0) .+
            (gbr_surv_means[non_missing_mask] .* (coefs .< 0.0))
        )

        if count(non_missing_mask) == 0
            continue
        end
        plotted_a_location = true

        lines!(ax, xs, coefs, label=nm)
    end
    if plotted_a_location
        axislegend(ax)
    end
    return f
end

function calculate_quantile(gbr_mean::Float64, reg_mean::Float64)
    alpha::Float64 = 2 / (1 - gbr_mean)
    beta::Float64 = 2 / gbr_mean
    dist = Beta(alpha, beta)
    return cdf(dist, reg_mean)
end

function plot_implied_quantile_survival(
    df::DataFrame,
    windows::WIN_TYPE,
    row_masks,
    row_names,
    title::String;
    count_threshold::Int64=20
)::Figure
    gbr_surv_means, gbr_surv_stdevs = calculate_survival(
        df, windows; count_threshold=count_threshold
    )

    _xs = [
        ((lb + ub) / 2) for (lb, ub) in zip(windows[1], windows[2])
    ]

    f = Figure(; size=(1200, 900))
    ax = Axis(
        f[1, 1],
        xlabel="Log Diameter (cm)",
        ylabel="Quantile",
        title=title,
        limits=(0.0, nothing, 0.0, 1.0)
    )
    lines!(_xs[[1, end]], [0.0, 0.0])
    plotted_a_location = false

    for (nm, r_msk) in zip(row_names, row_masks)
        surv_means, surv_stdevs = calculate_survival(
            df[r_msk, :], windows; count_threshold=count_threshold
        )

        non_missing_mask = (!).(ismissing.(surv_means)) .&& (!).(ismissing.(gbr_surv_means))
        xs             = _xs[non_missing_mask]
        quantiles = calculate_quantile.(
            gbr_surv_means[non_missing_mask],
            surv_means[non_missing_mask]
        )

        if count(non_missing_mask) == 0
            continue
        end
        plotted_a_location = true

        lines!(ax, xs, quantiles, label=nm)
    end
    if plotted_a_location
        axislegend(ax)
    end
    return f
end
