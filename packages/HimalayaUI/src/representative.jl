is_kept_exposure(e) = !(e.status isa AbstractString && e.status == "rejected")

function representative_exposure_id(exposures)
    kept = [Int(e.id) for e in exposures if is_kept_exposure(e)]
    isempty(kept) ? nothing : maximum(kept)
end
