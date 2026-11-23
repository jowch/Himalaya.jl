is_kept_exposure(e) = !(e.status isa AbstractString && e.status == "rejected")

shows_combined(e) = hasproperty(e, :kind) && is_combined_kind(e.kind) &&
    hasproperty(e, :sources) && !ismissing(e.sources) && e.sources >= MIN_COMBINED_SOURCES

function representative_exposure_id(exposures)
    combined = [Int(e.id) for e in exposures if shows_combined(e)]
    isempty(combined) || return first(combined)
    kept = [Int(e.id) for e in exposures if is_kept_exposure(e) && !(hasproperty(e, :kind) && is_combined_kind(e.kind))]
    isempty(kept) ? nothing : maximum(kept)
end
