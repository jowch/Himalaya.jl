using Printf

const DETECTORS = (:hybrid, :persistence, :recipe, :oracle)

_oracle_picks(exps) = HimalayaUI.oracle_picks(exps)

function _detector_picks(detector::Symbol, exps)
    picks = Vector{Vector{Int}}(undef, length(exps))
    method = detector
    Threads.@threads for k in eachindex(exps)
        e = exps[k]
        picks[k] = collect(Int, Himalaya.findpeaks(e.q, e.I, e.σ; method).indices)
    end
    picks
end

_threshold_list(::Nothing, recipe) = [recipe.threshold]
_threshold_list(t::Symbol, _) = t === :sweep ? collect(HimalayaUI.THRESHOLD_GRID) :
    throw(ArgumentError("thresholds must be a number or :sweep"))
_threshold_list(t::Real, _) = [Float64(t)]

function _entry(name, detector, threshold, exps, picks, ctx)
    ev = HimalayaUI.peak_eval(exps, picks; ctx.targets, ctx.label_ids)
    (; name, detector, threshold, ev.summary, ev.phase_rows, ev.peak_count, ev.ff_peak_count)
end

"""
    evaluate_database(db_path; detectors = [:hybrid], recipe = default_recipe(),
                      thresholds = nothing, phase_auto = :user_events, phase_targets = :reachable,
                      group_by = nothing, io = stderr) -> NamedTuple

Scores each detector (`:hybrid`, `:persistence`, `:recipe`, `:oracle`, the user's own peaks as a reference) on every loaded
exposure of the database, opened read-only, with the training gate's own
scoring (`HimalayaUI.peak_eval`). `:recipe` is out of fold: one model per held
out group, as in training. `thresholds` (`nothing` = the recipe's, a number, or
`:sweep` over the gate's grid) gives one `:recipe` entry per threshold.
`phase_targets = :reachable` scores a confirmed call only when the oracle's peaks or
some detector in this run give a candidate of its phase, so the targets depend on
the detector set (the oracle is always included); `:all` scores every call.
`group_by` defaults to what training would pick for this database. Returns
`counts`, `phase_auto`, `phase_targets`, `entries` and `paired` (every entry after the first
against the first).
"""
function evaluate_database(db_path::AbstractString; detectors = [:hybrid],
                           recipe::Recipe = default_recipe(), thresholds = nothing,
                           phase_auto::Symbol = :user_events, phase_targets::Symbol = :reachable,
                           group_by::Union{Symbol, Nothing} = nothing,
                           io::IO = stderr)
    detectors = Symbol[Symbol(d) for d in detectors]
    isempty(detectors) && throw(ArgumentError("no detectors given"))
    bad = setdiff(detectors, DETECTORS)
    isempty(bad) || throw(ArgumentError("unknown detector $(first(bad)); use hybrid, persistence, recipe or oracle"))
    allunique(detectors) || throw(ArgumentError("detectors must be distinct"))
    phase_auto in (:user_events, :all) || throw(ArgumentError("phase_auto must be :user_events or :all"))
    phase_targets in (:reachable, :all) || throw(ArgumentError("phase_targets must be :reachable or :all"))
    thresholds === nothing || :recipe in detectors ||
        @warn "thresholds only applies to the recipe detector, which was not requested; ignoring it"
    thrs = _threshold_list(thresholds, recipe)
    db = open_readonly(db_path)
    try
        t0 = time()
        if group_by === nothing
            group_by = something(HimalayaUI.cold_start_check(db).group_by, :sample)
        end
        exps = HimalayaUI.labeled_exposures(db; group_by, include_unlabelled = true, io)
        println(io, "loaded $(length(exps)) exposures in $(round(time() - t0; digits = 1)) s")
        ds = nothing
        if :recipe in detectors
            ds = build_dataset(exps; features = recipe.features, float32 = true)
            println(io, "built dataset: $(length(ds.label)) candidates in $(round(time() - t0; digits = 1)) s")
        end
        runs = Any[]
        for d in detectors
            if d === :hybrid
                push!(runs, ("hybrid", d, nothing, ds === nothing ? _detector_picks(:hybrid, exps) : ds.hybrid))
            elseif d === :oracle
                push!(runs, ("oracle", d, nothing, _oracle_picks(exps)))
            elseif d === :persistence
                push!(runs, ("persistence", d, nothing, _detector_picks(:persistence, exps)))
            else
                oof = oof_proba(recipe, ds)
                println(io, "out-of-fold probabilities in $(round(time() - t0; digits = 1)) s")
                rows_of = [Int[] for _ in exps]
                for (r, k) in enumerate(ds.exposure); push!(rows_of[k], r); end
                for t in thrs
                    picks = [ds.cand[k][HPM.select_peaks(ds.cand[k], oof[rows_of[k]], t, recipe.nms)]
                             for k in eachindex(exps)]
                    push!(runs, ("recipe@$t", d, t, picks))
                end
            end
        end
        ctx = HimalayaUI.eval_context(db, exps; phase_auto, phase_targets, detector_picks = [r[4] for r in runs])
        entries = Any[_entry(n, d, t, exps, picks, ctx) for (n, d, t, picks) in runs]
        counts = (exposures = length(exps), labelled = count(e -> e.labelled, exps),
                  form_factor = count(e -> e.form_factor, exps), group_by,
                  groups = length(unique(e.group for e in exps)),
                  phase_targets = length(ctx.targets), phase_calls = ctx.counts, label_exposures = length(ctx.label_ids))
        println(io, "total $(round(time() - t0; digits = 1)) s")
        (; counts, phase_auto, phase_targets, entries, paired = [_pair(entries[1], e) for e in entries[2:end]],
           recipe_id = :recipe in detectors && length(thrs) == 1 ? recipe_id(Recipe(recipe.features, recipe.learner, recipe.weak_keep, thrs[1], recipe.nms)) : nothing)
    finally
        close(db)
    end
end

"""Per phase (the exposure's first confirmed call): `(top-1 correct, exposures)`."""
function phase_table(e)
    out = Dict{String, Tuple{Int, Int}}()
    for r in e.phase_rows
        c, n = get(out, r.phase, (0, 0))
        out[r.phase] = (c + r.ok, n + 1)
    end
    out
end

function _pair(base, e)
    ok_base = Dict(r.id => r.ok for r in base.phase_rows)
    cell = Dict{String, NTuple{4, Int}}()
    for r in e.phase_rows
        a = ok_base[r.id]
        w, l, both, neither = get(cell, r.phase, (0, 0, 0, 0))
        cell[r.phase] = (w + (r.ok && !a), l + (!r.ok && a), both + (r.ok && a), neither + (!r.ok && !a))
    end
    (; name = e.name, against = base.name, phases = cell)
end

_frac(a, b) = b == 0 ? "-" : @sprintf("%d/%d (%.2f)", a, b, a / b)

function _metric_cells(e)
    s = e.summary
    ("phase top-1" => _frac(s.phase...),
     "add recall" => _frac(s.adds...),
     "excludes detected" => _frac(s.excluded...),
     "peaks/exposure" => _frac(e.peak_count, s.exposures),
     "peaks/form-factor" => s.form_factor == 0 ? "-" : _frac(e.ff_peak_count, s.form_factor))
end

function _print_grid(io, corner, cols, rows)
    w0 = maximum(length.([corner; first.(rows)]))
    ws = [maximum(length.([cols[j]; [r[2][j] for r in rows]])) for j in eachindex(cols)]
    println(io, rpad(corner, w0), "  ", join((lpad(c, ws[j]) for (j, c) in enumerate(cols)), "  "))
    for (label, cells) in rows
        println(io, rpad(label, w0), "  ", join((lpad(c, ws[j]) for (j, c) in enumerate(cells)), "  "))
    end
end

function print_evaluation(io::IO, r)
    c = r.counts
    println(io, "exposures: $(c.exposures) loaded ($(c.labelled) labelled, $(c.form_factor) form factor); ",
            "folds by $(c.group_by) ($(c.groups) groups)")
    println(io, "phase calls ($(r.phase_auto), targets $(r.phase_targets)): $(c.phase_calls.reachable) of ",
            "$(c.phase_calls.all) exposures with a confirmed call are reachable from the user's peaks or a scored detector; ",
            "scoring $(c.phase_targets), on $(c.label_exposures) labelled exposures shown without a model")
    println(io, "reachable/all by confirmed phase: ",
            join(("$p $(v.reachable)/$(v.all)" for (p, v) in sort(collect(c.phase_calls.per_phase); by = first)), ", "))
    r.recipe_id === nothing || println(io, "recipe $(r.recipe_id), out of fold")
    println(io)
    names = [e.name for e in r.entries]
    labels = first.(_metric_cells(r.entries[1]))
    cells = [last.(_metric_cells(e)) for e in r.entries]
    _print_grid(io, "numerator/denominator (rate)", names,
                [(l, [cells[j][i] for j in eachindex(names)]) for (i, l) in enumerate(labels)])
    band = [e.summary.adds_in_band for e in r.entries]
    all(==(0), band) || println(io, "adds in the trimmed high-q band, not scored: ", join(unique(band), ", "))
    println(io)
    tables = phase_table.(r.entries)
    phases = sort(collect(union((keys(t) for t in tables)...)))
    println(io, "top-1 phase calls correct/total, by confirmed phase")
    _print_grid(io, "phase", names,
                [(p, [haskey(t, p) ? @sprintf("%d/%d", t[p]...) : "-" for t in tables]) for p in phases])
    isempty(r.paired) && return
    println(io)
    println(io, "paired against $(r.paired[1].against): wins/losses per phase (exposures only the row detector, or only the first, got right)")
    rows = [(p, [(w = get(x.phases, p, (0, 0, 0, 0)); @sprintf("%d/%d", w[1], w[2])) for x in r.paired]) for p in phases]
    push!(rows, ("all", [(w = reduce((a, b) -> a .+ b, values(x.phases); init = (0, 0, 0, 0));
                         @sprintf("%d/%d", w[1], w[2])) for x in r.paired]))
    _print_grid(io, "phase", [x.name for x in r.paired], rows)
end

_rate(a, b) = b == 0 ? nothing : a / b
_ratio_dict(c, n) = Dict("correct" => c, "total" => n, "rate" => _rate(c, n))
_count_dict(c, n) = Dict("count" => c, "total" => n, "rate" => _rate(c, n))

function evaluation_dict(r)
    c = r.counts
    Dict{String, Any}(
        "counts" => Dict{String, Any}("exposures" => c.exposures, "labelled" => c.labelled,
            "form_factor" => c.form_factor, "group_by" => String(c.group_by), "groups" => c.groups,
            "phase_targets" => c.phase_targets, "label_exposures" => c.label_exposures,
            "phase_calls" => Dict{String, Any}("all" => c.phase_calls.all, "reachable" => c.phase_calls.reachable,
                "per_phase" => Dict{String, Any}(p => Dict("reachable" => v.reachable, "all" => v.all)
                                                 for (p, v) in c.phase_calls.per_phase))),
        "phase_auto" => String(r.phase_auto),
        "phase_targets" => String(r.phase_targets),
        "recipe_id" => r.recipe_id,
        "detectors" => [Dict{String, Any}("name" => e.name, "detector" => String(e.detector), "threshold" => e.threshold,
            "metrics" => Dict{String, Any}("phase" => _ratio_dict(e.summary.phase...),
                "adds" => _ratio_dict(e.summary.adds...),
                "excluded_detected" => _ratio_dict(e.summary.excluded...),
                "adds_in_band" => e.summary.adds_in_band,
                "peaks_per_exposure" => _count_dict(e.peak_count, e.summary.exposures),
                "peaks_per_form_factor" => _count_dict(e.ff_peak_count, e.summary.form_factor)),
            "per_phase" => Dict{String, Any}(p => _ratio_dict(v...) for (p, v) in phase_table(e))) for e in r.entries],
        "paired" => [Dict{String, Any}("name" => x.name, "against" => x.against,
            "phases" => Dict{String, Any}(p => Dict("wins" => v[1], "losses" => v[2], "both" => v[3], "neither" => v[4])
                                          for (p, v) in x.phases)) for x in r.paired])
end

function write_evaluation(path::AbstractString, r)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON3.write(io, evaluation_dict(r))
    end
    path
end
