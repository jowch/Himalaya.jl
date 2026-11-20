is_kept_exposure(e) = !(e.status isa AbstractString && e.status == "rejected")

function representative_exposure_id(exposures)
    kept = [e for e in exposures if is_kept_exposure(e)]
    isempty(kept) && return nothing
    picked = [e for e in kept if coalesce(e.selected, 0) == 1]
    maximum(e -> Int(e.id), isempty(picked) ? kept : picked)
end
