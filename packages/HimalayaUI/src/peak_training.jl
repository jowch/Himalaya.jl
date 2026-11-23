using SHA

const GATE_PEAKS_SLACK = 0.5
const GATE_FF_SLACK = 0.1
const PHASE_BASIS_TOL = 0.02
const THRESHOLD_GRID = 0.3:0.05:0.9

peak_feature_cache_dir(path::AbstractString) =
    joinpath(dirname(abspath(path)), "cache", "peak-features", splitext(basename(path))[1])
function peak_feature_cache_dir(db::SQLite.DB)
    f = db.file
    (f === nothing || isempty(f) || f == ":memory:" || startswith(f, "file:")) ? nothing : peak_feature_cache_dir(f)
end

function cold_start_check(db::SQLite.DB; min_curations::Integer = 100, min_exposures::Integer = 30,
                          min_experiments::Integer = 3, min_samples::Integer = 5)
    r = only(sqlrows(DBInterface.execute(db, """
        SELECT COUNT(*) AS curations, COUNT(DISTINCT c.exposure_id) AS curated_exposures,
               COUNT(DISTINCT e.experiment_id) AS experiments, COUNT(DISTINCT e.sample_id) AS samples
          FROM peak_curations c JOIN exposures e ON e.id = c.exposure_id
         WHERE c.kind IN ('add', 'exclude') AND e.status IS NOT 'rejected' AND e.kind = 'file'""")))
    counts = (curations = Int(r.curations), curated_exposures = Int(r.curated_exposures),
              experiments = Int(r.experiments), samples = Int(r.samples))
    reasons = String[]
    counts.curations >= min_curations ||
        push!(reasons, "needs at least $min_curations add and exclude curations, has $(counts.curations)")
    counts.curated_exposures >= min_exposures ||
        push!(reasons, "needs at least $min_exposures curated exposures, has $(counts.curated_exposures)")
    group_by = counts.experiments >= min_experiments ? :experiment :
               counts.samples >= min_samples ? :sample : nothing
    group_by === nothing &&
        push!(reasons, "needs at least $min_experiments experiments or $min_samples samples, has $(counts.experiments) experiments and $(counts.samples) samples")
    (ok = isempty(reasons), reasons, counts, group_by)
end

_check(name, passed, skipped, reason) = (; name, passed, skipped, reason = passed ? "" : reason)
_rate(t) = t[2] == 0 ? NaN : t[1] / t[2]
_f(x, d = 2) = string(round(x; digits = d))

function evaluate_gate(model, hybrid; peaks_slack::Real = GATE_PEAKS_SLACK, ff_slack::Real = GATE_FF_SLACK)
    ε = 1e-9
    checks = Any[]
    cap(name, label, m, h, slack) = push!(checks, isnan(m) || isnan(h) ? _check(name, true, true, "") :
        _check(name, m <= h + slack + ε, false,
               "$label $(_f(m)) is more than $slack above the hybrid's $(_f(h))"))
    cap("peaks_per_exposure", "peaks per exposure", model.peaks, hybrid.peaks, peaks_slack)
    cap("peaks_per_form_factor", "peaks per form-factor exposure", model.ff_peaks, hybrid.ff_peaks, ff_slack)
    push!(checks, model.phase[2] == 0 ? _check("phase_calls", true, true, "") :
        _check("phase_calls", model.phase[1] >= hybrid.phase[1], false,
               "top-1 phase calls $(model.phase[1])/$(model.phase[2]) are below the hybrid's $(hybrid.phase[1])/$(hybrid.phase[2])"))
    push!(checks, model.adds[2] == 0 ? _check("add_recall", true, true, "") :
        _check("add_recall", _rate(model.adds) >= _rate(hybrid.adds) - ε, false,
               "add recall $(_f(_rate(model.adds))) is below the hybrid's $(_f(_rate(hybrid.adds)))"))
    push!(checks, model.excluded[2] == 0 ? _check("excludes_detected", true, true, "") :
        _check("excludes_detected", _rate(model.excluded) <= _rate(hybrid.excluded) + ε, false,
               "excludes still detected $(_f(_rate(model.excluded))) exceed the hybrid's $(_f(_rate(hybrid.excluded)))"))
    (passed = all(c -> c.passed, checks), checks, reasons = String[c.reason for c in checks if !c.passed])
end

function top_call_matches(q::AbstractVector, sharpness::AbstractVector, confirmed)
    isempty(q) && return false
    candidates = Himalaya.indexpeaks(Float64.(q), Float64.(sharpness))
    isempty(candidates) && return false
    top = first(auto_group(candidates))
    name = string(nameof(Himalaya.phase(top)))
    any(p == name && abs(Himalaya.basis(top) / b - 1) < PHASE_BASIS_TOL for (p, b) in confirmed)
end

"""
Which exposures a comparison scores: `label_ids` are labelled exposures whose
latest direct review was shown without a model, `targets` their confirmed calls.
"""
function eval_context(db::SQLite.DB, exps; phase_auto::Symbol = :user_events)
    shown = no_model_shown_ids(db)
    label_ids = intersect(Set(e.id for e in exps if e.labelled), shown)
    targets = confirmed_calls(db, [e.id for e in exps if e.labelled && e.id in shown]; auto = phase_auto)
    (; label_ids, targets)
end

function confirmed_calls(db::SQLite.DB, exposure_ids; auto::Symbol = :user_events)
    auto in (:user_events, :all) || throw(ArgumentError("auto must be :user_events or :all"))
    out = Dict{Int, Vector{Tuple{String, Float64}}}()
    wanted = Set(Int.(exposure_ids))
    for r in sqlrows(DBInterface.execute(db, """
            SELECT m.exposure_id AS eid, i.phase AS phase, i.basis AS basis
              FROM assignment_members m
              JOIN assignments a ON a.exposure_id = m.exposure_id AND a.state = 'indexed'
              JOIN indices i ON i.id = m.index_id AND i.kind IN ('auto', 'speculative')
              JOIN exposures e ON e.id = m.exposure_id
              LEFT JOIN sample_calls sc ON sc.sample_id = e.sample_id
             WHERE COALESCE(sc.state, 'indexed') = 'indexed'
               AND (i.kind = 'speculative' OR $(auto === :all ? "1" : "EXISTS (SELECT 1 FROM user_actions ua WHERE ua.entity_type = 'exposure' AND ua.entity_id = m.exposure_id AND ua.action IN ('assignment_add', 'assignment_remove', 'assignment_set_state'))"))
             ORDER BY m.exposure_id, i.id"""))
        Int(r.eid) in wanted && push!(get!(out, Int(r.eid), Tuple{String, Float64}[]), (String(r.phase), Float64(r.basis)))
    end
    out
end

no_model_shown_ids(db::SQLite.DB) = Set(Int(r.exposure_id) for r in sqlrows(DBInterface.execute(db, """
    SELECT r.exposure_id FROM peak_reviews r
     WHERE r.direct = 1 AND r.shown_model_id IS NULL
       AND r.id = (SELECT MAX(id) FROM peak_reviews WHERE exposure_id = r.exposure_id AND direct = 1)""")))

function _median(v)
    s = sort(v); n = length(s)
    isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

function _sharpness_at(sfull, idx)
    s = sfull[idx]
    pos = filter(>(0), s)
    fill = isempty(pos) ? eps() : _median(pos)
    [x > 0 ? x : fill for x in s]
end

function _near(x, set, step)
    tol = HimalayaPeakModel.peak_match_tol(x, step)
    any(p -> abs(x - p) <= tol, set)
end

function peak_summary(exps, picks; kwargs...)
    peak_eval(exps, picks; kwargs...).summary
end

"""
Scores `picks` (grid indices per exposure). Returns the `summary` the training
gate compares, one `phase_rows` entry per exposure with a confirmed call (its
first confirmed phase and whether the top call matches any confirmed call), and
the total peak counts behind the per-exposure rates.
"""
function peak_eval(exps, picks; targets = Dict{Int, Vector{Tuple{String, Float64}}}(),
                   label_ids = Set(e.id for e in exps if e.labelled))
    phase_rows = NamedTuple{(:id, :phase, :ok), Tuple{Int, String, Bool}}[]
    adds_hit = adds_n = band = excl_hit = excl_n = 0
    phase_ok = phase_n = 0
    ff_n = 0; ff_sum = 0
    for (e, pk) in zip(exps, picks)
        step = HimalayaPeakModel.grid_step(e.q)
        pq = e.q[pk]
        if e.form_factor
            ff_n += 1; ff_sum += length(pk)
        end
        if e.id in label_ids
            qlo, qhi = extrema(e.q)
            qcut = qlo + (1 - HimalayaPeakModel.Q_TRIM_HIGH) * (qhi - qlo)
            for a in e.adds
                a > qcut && (band += 1; continue)
                adds_n += 1; adds_hit += _near(a, pq, step)
            end
            for x in e.excluded
                excl_n += 1; excl_hit += _near(x, pq, step)
            end
        end
        calls = get(targets, e.id, nothing)
        if calls !== nothing
            keep = [i for i in pk if !_near(e.q[i], e.excluded, step)]
            ok = top_call_matches(e.q[keep], _sharpness_at(Himalaya.sharpness(e.I), keep), calls)
            phase_n += 1; phase_ok += ok
            push!(phase_rows, (id = e.id, phase = first(calls)[1], ok = ok))
        end
    end
    summary = (peaks = isempty(exps) ? NaN : sum(length, picks) / length(exps),
               ff_peaks = ff_n == 0 ? NaN : ff_sum / ff_n,
               phase = (phase_ok, phase_n), adds = (adds_hit, adds_n), excluded = (excl_hit, excl_n),
               adds_in_band = band, exposures = length(exps), form_factor = ff_n)
    (; summary, phase_rows, peak_count = sum(length, picks; init = 0), ff_peak_count = ff_sum)
end

function _json_safe(x)
    x isa AbstractFloat ? (isfinite(x) ? x : nothing) :
    x isa NamedTuple ? Dict{String, Any}(String(k) => _json_safe(v) for (k, v) in pairs(x)) :
    x isa AbstractDict ? Dict{String, Any}(String(k) => _json_safe(v) for (k, v) in x) :
    x isa Union{Tuple, AbstractVector} ? Any[_json_safe(v) for v in x] :
    x isa Symbol ? String(x) : x
end

_content_key(e) = bytes2hex(sha256(vcat(reinterpret(UInt8, e.q), reinterpret(UInt8, e.I), reinterpret(UInt8, e.σ))))

function run_peak_training(db::SQLite.DB; recipe::HimalayaPeakModel.Recipe = HimalayaPeakModel.default_recipe(),
                           cache_dir::Union{AbstractString, Nothing} = nothing, io::IO = stderr,
                           cold_start = (;), phase_auto::Symbol = :user_events, peaks_slack = GATE_PEAKS_SLACK, ff_slack = GATE_FF_SLACK)
    t0 = time()
    watermark = label_watermark(db)
    base = Dict{String, Any}("watermark" => watermark, "core_version" => string(pkgversion(Himalaya)),
        "feature_version" => HimalayaPeakModel.feature_version(recipe.features))
    function finish(status, counts, recipe_used, metrics, fitted = nothing)
        merge(base, Dict{String, Any}("status" => status, "counts" => _json_safe(counts),
            "recipe" => HimalayaPeakModel.to_dict(recipe_used), "recipe_id" => HimalayaPeakModel.recipe_id(recipe_used),
            "threshold" => recipe_used.threshold, "nms" => recipe_used.nms,
            "metrics" => _json_safe(metrics), "fitted" => fitted))
    end
    cold = cold_start_check(db; cold_start...)
    if !cold.ok
        return finish("rejected", cold.counts, recipe,
            (gate = (passed = false, reasons = cold.reasons, checks = []), cold_start = cold.counts))
    end

    exps = labeled_exposures(db; group_by = cold.group_by, include_unlabelled = true, io)
    println(io, "loaded $(length(exps)) exposures in $(round(time() - t0; digits = 1)) s")
    features = recipe.features
    cache = nothing
    keys = Dict(e.id => _content_key(e) for e in exps)
    if cache_dir !== nothing
        cache = HimalayaPeakModel.DiskFeatureCache(String(cache_dir), HimalayaPeakModel.feature_version(features), e -> keys[e.id])
    end
    ds = HimalayaPeakModel.build_dataset(exps; features, cache, float32 = true)
    cache === nothing || HimalayaPeakModel.prune!(cache, values(keys))
    labels = Dict(n => count(==(l), ds.label) for (n, l) in
        ("positive" => HimalayaPeakModel.POS, "negative" => HimalayaPeakModel.NEG,
         "weak_negative" => HimalayaPeakModel.WEAK, "dropped" => HimalayaPeakModel.DROP))
    counts = (; cold.counts..., group_by = cold.group_by, exposures = length(exps),
        labelled_exposures = count(e -> e.labelled, exps), form_factor_exposures = count(e -> e.form_factor, exps),
        candidates = length(ds.label), labels)
    println(io, "built dataset: $(counts.candidates) candidates in $(round(time() - t0; digits = 1)) s")

    oof = HimalayaPeakModel.oof_proba(recipe, ds; verbose = false)
    println(io, "out-of-fold probabilities in $(round(time() - t0; digits = 1)) s")
    cv = HimalayaPeakModel.cross_validate(exps, recipe; ds, oof, verbose = false)

    rows_of = [Int[] for _ in exps]
    for (r, k) in enumerate(ds.exposure); push!(rows_of[k], r); end
    picks_at(thr) = [ds.cand[k][HimalayaPeakModel.select_peaks(ds.cand[k], oof[rows_of[k]], thr, recipe.nms)]
                     for k in eachindex(exps)]
    (; label_ids, targets) = eval_context(db, exps; phase_auto)
    hybrid = peak_summary(exps, ds.hybrid; targets, label_ids)

    sweep = [(threshold = t, summary = peak_summary(exps, picks_at(t); label_ids)) for t in THRESHOLD_GRID]
    within(s) = s.peaks <= hybrid.peaks + peaks_slack + 1e-9 &&
                (isnan(s.ff_peaks) || isnan(hybrid.ff_peaks) || s.ff_peaks <= hybrid.ff_peaks + ff_slack + 1e-9)
    ok = [x for x in sweep if within(x.summary)]
    threshold = isempty(ok) ? last(THRESHOLD_GRID) :
        sort(ok; by = x -> (isnan(_rate(x.summary.adds)) ? 0.0 : _rate(x.summary.adds), x.threshold))[end].threshold
    model = peak_summary(exps, picks_at(threshold); targets, label_ids)
    gate = evaluate_gate(model, hybrid; peaks_slack, ff_slack)
    isempty(ok) && pushfirst!(gate.reasons, "no threshold keeps peaks within the caps")
    gate = (; gate..., passed = gate.passed && !isempty(ok))

    final = HimalayaPeakModel.Recipe(recipe.features, recipe.learner, recipe.weak_keep, threshold, recipe.nms)
    metrics = (
        cv = (candidate_auc = cv.candidate_auc, groups = cv.groups,
              by_threshold = [(; threshold = t, s...) for (t, s) in cv.by_threshold],
              hybrid = cv.hybrid, candidates = cv.candidates),
        comparison = (threshold = threshold, phase_auto = phase_auto, model = model, hybrid = hybrid, phase_targets = length(targets),
                      label_exposures_no_model_shown = length(label_ids)),
        sweep = [(; threshold = x.threshold, peaks = x.summary.peaks, ff_peaks = x.summary.ff_peaks,
                  adds = x.summary.adds, within_caps = within(x.summary)) for x in sweep],
        gate = (; gate..., caps = (peaks_per_exposure = peaks_slack, peaks_per_form_factor = ff_slack)),
        seconds = round(time() - t0; digits = 1))
    gate.passed || return finish("rejected", counts, final, metrics)

    trained = HimalayaPeakModel.train(final, ds)
    println(io, "trained the final model in $(round(time() - t0; digits = 1)) s")
    finish("passed", merge(counts, (; positive = trained.trained_on.positive, negative = trained.trained_on.negative)),
           final, merge(metrics, (; seconds = round(time() - t0; digits = 1))), HimalayaPeakModel.to_dict(trained.fitted))
end
