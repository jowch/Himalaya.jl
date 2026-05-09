# validate.jl — pure-function manifest validation. No DB, no IO.

const _VALID_NAME_REGEX = r"^[A-Za-z0-9._-]+$"

"""
    ManifestViolation(kind, sample_index, sample_name, detail)

One validation problem found in a parsed manifest. `kind` is one of
`:empty_name`, `:bad_name_chars`, `:duplicate_name`,
`:duplicate_filename_in_sample`, `:overlapping_filenames`.
"""
struct ManifestViolation
    kind         ::Symbol
    sample_index ::Int        # 1-based row in manifest
    sample_name  ::String     # may be "" for :empty_name
    detail       ::String     # human-readable specifics
end

"""
    validate_manifest(samples) -> Vector{ManifestViolation}

Apply all five rules, collect every violation, return them in stable order
(rule 1, rule 2, rule 3, rule 4, rule 5). No fail-fast — operator sees every
fix needed in one pass.
"""
function validate_manifest(samples::Vector{ManifestSample})::Vector{ManifestViolation}
    out = ManifestViolation[]

    # Rule 1: name non-empty.
    for (i, s) in enumerate(samples)
        if isempty(s.name)
            push!(out, ManifestViolation(:empty_name, i, "", "row $i has an empty name"))
        end
    end

    # Rule 2: name matches allowed character set.
    for (i, s) in enumerate(samples)
        isempty(s.name) && continue  # already caught by rule 1
        if !occursin(_VALID_NAME_REGEX, s.name)
            push!(out, ManifestViolation(
                :bad_name_chars, i, s.name,
                "sample name \"$(s.name)\" contains characters outside [A-Za-z0-9._-]"
            ))
        end
    end

    # Rule 3: name unique within manifest.
    seen_names = Dict{String, Int}()  # name => first-seen index
    for (i, s) in enumerate(samples)
        isempty(s.name) && continue  # empty names already caught
        if haskey(seen_names, s.name)
            push!(out, ManifestViolation(
                :duplicate_name, i, s.name,
                "sample name \"$(s.name)\" at row $i duplicates row $(seen_names[s.name])"
            ))
        else
            seen_names[s.name] = i
        end
    end

    # Rule 4: filename ranges within a sample have no internal duplicates after expansion.
    for (i, s) in enumerate(samples)
        seen = Set{String}()
        for fn in s.filenames
            if fn in seen
                push!(out, ManifestViolation(
                    :duplicate_filename_in_sample, i, s.name,
                    "filename \"$fn\" appears more than once in sample \"$(s.name)\" (row $i)"
                ))
                break  # one violation per sample is enough
            end
            push!(seen, fn)
        end
    end

    # Rule 5: filename ranges across samples don't overlap.
    # Build a map: filename => (sample_index, sample_name) for the first occurrence.
    claimed = Dict{String, Tuple{Int, String}}()  # filename => (index, name)
    for (i, s) in enumerate(samples)
        for fn in s.filenames
            if haskey(claimed, fn)
                (j, other_name) = claimed[fn]
                push!(out, ManifestViolation(
                    :overlapping_filenames, i, s.name,
                    "filename \"$fn\" in sample \"$(s.name)\" (row $i) already claimed by sample \"$other_name\" (row $j)"
                ))
            else
                claimed[fn] = (i, s.name)
            end
        end
    end

    out
end
