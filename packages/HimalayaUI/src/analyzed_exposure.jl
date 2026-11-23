is_kept_exposure(e) = !(e.status isa AbstractString && e.status == "rejected")

_frame_kind(e) = hasproperty(e, :kind) && e.kind isa AbstractString ? e.kind : "file"

_frame_selected(e) =
    !(hasproperty(e, :selected) && e.selected === false) &&
    !(hasproperty(e, :deselected) && !ismissing(e.deselected) && e.deselected != 0)

_frame_sources(e) = hasproperty(e, :sources) && !ismissing(e.sources) && !(e.sources isa Number) ?
    sort(Int[Int(s) for s in e.sources]) : Int[]

function analyzed_exposure_id(frames)
    sel = sort([Int(e.id) for e in frames
                if is_kept_exposure(e) && !is_combined_kind(_frame_kind(e)) && _frame_selected(e)])
    isempty(sel) && return nothing
    length(sel) == 1 && return sel[1]
    for e in frames
        is_combined_kind(_frame_kind(e)) && _frame_sources(e) == sel && return Int(e.id)
    end
    nothing
end

function analysis_gap(frames)
    images = [e for e in frames if !is_combined_kind(_frame_kind(e))]
    isempty(images) && return nothing
    kept = filter(is_kept_exposure, images)
    isempty(kept) && return "all_dropped"
    any(_frame_selected, kept) ? nothing : "none_selected"
end
