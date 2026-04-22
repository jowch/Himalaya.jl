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
abstract type Cubic <: Phase end
abstract type Lamellar <: Phase end
abstract type Hexagonal <: Phase end
abstract type Square <: Phase end
abstract type Pn3m <: Cubic end
abstract type Im3m <: Cubic end
abstract type Ia3d <: Cubic end
abstract type Fm3m <: Cubic end
abstract type Fd3m <: Cubic end

"""
    phaseratios(::Phase; normalize = false)

Returns the expected peak ratios for the given `Phase`. If `normalize`, then
return the ratios divided by the first ratio so that the first ratio is 1.
"""
function phaseratios(::Type{P}; normalize = false) where {P<:Phase}
    if normalize
        ratios = phaseratios(P)
        return ratios ./ first(ratios)
    end

    phaseratios(P)
end

phaseratios(::Type{Lamellar}) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
phaseratios(::Type{Hexagonal}) = [1, √3, √4, √7, √9, √11, √12, √13, √16, √19,
                                  √21, √25, √27, √28]
phaseratios(::Type{Square}) = [1, √2, √4, √5, √8, √9, √10, √13, √16, √17, √18, √20]

# Cubic phases
phaseratios(::Type{Pn3m}) = [√2, √3, √4, √6, √8, √9, √10, √11, √12, √14, √16,
                             √17, √18, √19, √20, √21] #, √22, √24, √26, √27]
phaseratios(::Type{Im3m}) = [√2, √4, √6, √8, √10, √12, √14, √16, √18, √20] #, √22, √24, √26]
phaseratios(::Type{Ia3d}) = [√6, √8, √14, √16, √20, √22, √24, √26]

# Micellar phases
phaseratios(::Type{Fm3m}) = [√3, √4, √8, √11, √12]
phaseratios(::Type{Fd3m}) = [√3, √8, √11, √12, √16, √19, √24, √27, √32, √35, √36]

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
minpeaks(::Type{Square}) = 3
minpeaks(::Type{Pn3m}) = 4
minpeaks(::Type{Im3m}) = 4
minpeaks(::Type{Ia3d}) = 4
minpeaks(::Type{Fm3m}) = 4
minpeaks(::Type{Fd3m}) = 6
