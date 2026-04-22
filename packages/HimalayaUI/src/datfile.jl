using DelimitedFiles

"""
    load_dat(path) -> (q, I, σ)

Parse a three-column whitespace-separated SAXS integration file.
Returns vectors of scattering vector q, intensity I, and uncertainty σ.
No header rows expected.
"""
function load_dat(path::AbstractString)
    data = readdlm(path, Float64)
    size(data, 2) >= 3 || error("$path: expected ≥3 columns, got $(size(data,2))")
    data[:, 1], data[:, 2], data[:, 3]
end
