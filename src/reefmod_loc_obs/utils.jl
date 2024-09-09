function lb_simple(encoding::AbstractString)::Union{Missing, Float64}
    tmp = encoding
    if length(encoding) == 1 && encoding[1] != '0'
        tmp *= "L"
    end
    encodings::Dict{String, Union{Missing, Float64}} = Dict(
        "" => missing,
        "0"  => 0.0,
        "1L" => 1.0,
        "1U" => 5.0,
        "2L" => 10.0,
        "2U" => 20.0,
        "3L" => 30.0,
        "3U" => 40.0,
        "4L" => 50.0,
        "4U" => 62.5,
        "5L" => 75.0,
        "5U" => 87.5
    )
    return missing
end
function ub_simple(encoding::AbstractString)::Union{Missing, Float64}
    tmp = encoding
    if length(encoding) == 1 && encoding[1] != '0'
        tmp *= "U"
    end
    encodings::Dict{String, Union{Missing, Float64}} = Dict(
        "" => missing,
        "0"  => 0.0,
        "1L" => 5.0,
        "1U" => 10.0,
        "2L" => 20.0,
        "2U" => 30.0,
        "3L" => 40.0,
        "3U" => 50.0,
        "4L" => 62.5,
        "4U" => 75.0,
        "5L" => 87.5,
        "5U" => 100.0
    )
    return encodings[tmp]
end

function lower_bound(encoding::Union{Missing, AbstractString})::Union{Missing, Float64}
    ismissing(encoding) && return missing
    tmp = encoding
    if contains(tmp, '/')
        tmp = split(tmp, '/')[1]
    end
    return lb_simple(tmp)
end

function upper_bound(encoding::Union{Missing, AbstractString})::Union{Missing, Float64}
    ismissing(encoding) && return missing
    tmp = encoding
    if contains(encoding, '/')
        tmp = split(encoding, '/')[2]
    end
    return ub_simple(tmp)
end
