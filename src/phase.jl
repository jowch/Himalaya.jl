"""
    Phase

A Phase describes the liquid crystalline organization for some material. The
material's organization can result in a diffraction pattern when exposed to
light that is unique to it's phase. Each phase's diffraction pattern is
described by a set of ratios of the pattern's peaks divided by the first peak
(basis).

┌────────────────────────────────────────┐
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⢸⠀⠀⠀⠀⠀⢸⠀⠀⡆⢰⡇⠀⠀⠀⠀⠀⣤⠟⠦⡴⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠘⡄⠀⠀⠀⠀⢸⠀⢰⡇⢸⡇⠀⠀⣷⠀⢰⠃⠀⠀⠀⠀⠙⠲⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⡇⠀⠀⠀⠀⢸⠀⢸⡇⢸⡇⠀⠀⣿⡴⠃⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⢇⠀⠀⠀⠀⢸⠀⢸⡇⡜⢣⠀⣰⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠳⣄⠀⠀⠀⠀⠀⠀⠀│
│⠀⢸⠀⠀⠀⠀⣿⠀⡼⠑⠃⠘⠚⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀│
│⠀⠈⡇⠀⠀⠀⡏⠾⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⡀⠀⠀⠀│
│⠀⠀⢣⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀│
│⠀⠀⠸⡄⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⡀│
│⠀⠀⠀⢧⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹│
└────────────────────────────────────────┘

# Subtypes
    Lamellar
    Hexagonal
    Pn3m
    Im3m
    Ia3d
    Fd3m
"""
abstract type Phase end
abstract type Lamellar <: Phase end
abstract type Hexagonal <: Phase end
abstract type Pn3m <: Phase end
abstract type Im3m <: Phase end
abstract type Ia3d <: Phase end
abstract type Fd3m <: Phase end

"""
    phaseratios(::Phase)

Returns the expected peak ratios for the given `Phase`.
"""
function phaseratios(::Type{P}) where {P<:Phase}
    phaseratios(P)
end

phaseratios(::Type{Lamellar}) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
phaseratios(::Type{Hexagonal}) = [1, √3, √4, √7, √9, √11, √12]
phaseratios(::Type{Pn3m}) = [√2, √3, √4, √6, √8, √9, √10, √11] ./ √2
phaseratios(::Type{Im3m}) = [√2, √4, √6, √8, √10, √12, √14, √16, √18] ./ √2
phaseratios(::Type{Ia3d}) = [√6, √8, √14, √16, √20, √22, √24, √26] ./ √6
phaseratios(::Type{Fd3m}) = [√3, √8, √11, √12, √16, √19, √24, √27, √32, √35, √36] ./ √3

"""
    minpeaks(::Phase)

Returns the minimum number of peaks to assign a given `Phase`.

See also `indexpeaks`.
"""
function minpeaks(::Type{P}) where {P<:Phase}
    minpeaks(P)
end

minpeaks(::Type{Lamellar}) = 2
minpeaks(::Type{Hexagonal}) = 3
minpeaks(::Type{Pn3m}) = 4
minpeaks(::Type{Im3m}) = 4
minpeaks(::Type{Ia3d}) = 4
minpeaks(::Type{Fd3m}) = 4
