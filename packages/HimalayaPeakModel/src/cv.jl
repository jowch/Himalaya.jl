mutable struct Tally
    adds::Int; adds_hit::Int
    kept::Int; kept_hit::Int
    excl::Int; excl_hit::Int
    curated::Int; peaks::Int; extra::Int
    ff::Int; ff_peaks::Int
end
Tally() = Tally(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

function tally!(t::Tally, e::LabeledExposure, peak_q::AbstractVector)
    q_step = grid_step(e.q)
    near(x, set) = any(p -> abs(x - p) <= match_tol(x, q_step), set)
    if e.curated
        t.curated += 1
        t.adds += length(e.adds);       t.adds_hit += count(a -> near(a, peak_q), e.adds)
        t.kept += length(e.kept);       t.kept_hit += count(a -> near(a, peak_q), e.kept)
        t.excl += length(e.excluded);   t.excl_hit += count(a -> near(a, peak_q), e.excluded)
        pos = vcat(e.adds, e.kept)
        t.peaks += length(peak_q)
        t.extra += count(x -> !near(x, pos), peak_q)
    end
    if e.form_factor
        t.ff += 1; t.ff_peaks += length(peak_q)
    end
    t
end

ratio(a, b) = b == 0 ? NaN : a / b

"""Rank-based area under the ROC curve; `NaN` if either class is empty."""
function auc(score, y)
    n1 = count(y); n0 = length(y) - n1
    (n1 == 0 || n0 == 0) && return NaN
    r = tiedrank_(score)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
end

function tiedrank_(x)
    o = sortperm(x); r = zeros(length(x)); i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[o[j + 1]] == x[o[i]]; j += 1; end
        r[o[i:j]] .= (i + j) / 2
        i = j + 1
    end
    r
end

"""
    cross_validate(exps, ds; thresholds = [0.3, 0.5, 0.7, 0.85, 0.9], kwargs...) -> NamedTuple

Leave-one-group-out cross-validation (groups are `LabeledExposure.group`,
an experiment by default). Each held-out group is scored by a model that
never saw it. Reports, per threshold, the classifier against the hybrid
detector and against the unfiltered candidate pool (the recall ceiling), on
the same exposures: recall of adds, kept peaks and excluded peaks (an
excluded peak that is still detected is a miss), and predicted peaks per
curated exposure that match neither an add nor a kept peak (`extra`;
unlabelled real peaks count here, so this overstates false positives).
`kwargs` go to `train`; its `nms` is applied.
"""
function cross_validate(exps::AbstractVector{LabeledExposure}, ds::Dataset;
                        thresholds = [0.3, 0.5, 0.7, 0.85, 0.9], verbose::Bool = true, kwargs...)
    nms = get(Dict(kwargs), :nms, 8)
    groups = sort(unique(ds.group))
    length(groups) >= 2 || error("cross-validation needs at least two groups, found $(length(groups))")
    oof = fill(NaN, length(ds.label))
    for g in groups
        tr = findall(!=(g), ds.group)
        m = train(ds; rows = tr, kwargs...)
        te = findall(==(g), ds.group)
        oof[te] = proba(m, ds.X[te, :])
        verbose && @printf("  held out group %d: %d exposures\n", g, length(unique(ds.exposure[te])))
    end
    rows_of = [Int[] for _ in exps]
    for (r, k) in enumerate(ds.exposure); push!(rows_of[k], r); end
    function run(pick)
        T = Tally()
        for k in eachindex(exps)
            tally!(T, exps[k], exps[k].q[pick(k)])
        end
        T
    end
    summary(T) = (add_recall = ratio(T.adds_hit, T.adds), adds = (T.adds_hit, T.adds),
        kept_recall = ratio(T.kept_hit, T.kept), kept = (T.kept_hit, T.kept),
        excluded_still_detected = ratio(T.excl_hit, T.excl), excluded = (T.excl_hit, T.excl),
        peaks_per_curated = ratio(T.peaks, T.curated), extra_per_curated = ratio(T.extra, T.curated),
        peaks_per_form_factor = ratio(T.ff_peaks, T.ff), curated = T.curated, form_factor = T.ff)
    by_threshold = map(thresholds) do thr
        thr => summary(run(k -> (idx = ds.cand[k];
            idx[select_peaks(idx, oof[rows_of[k]], thr, nms)])))
    end
    scored = findall(l -> l == POS || l == NEG || l == WEAK, ds.label)
    (by_threshold = by_threshold,
     hybrid = summary(run(k -> ds.hybrid[k])),
     candidates = summary(run(k -> ds.cand[k])),
     candidate_auc = auc(oof[scored], ds.label[scored] .== POS), groups = groups)
end

function print_cv(r)
    f(x) = (x[1] == 0 && x[2] == 0) ? "-" : @sprintf("%d/%d", x...)
    @printf("candidate AUC (positive vs any negative, out of fold): %.3f\n", r.candidate_auc)
    println(rpad("", 16), rpad("add recall", 18), rpad("kept recall", 18), rpad("excl detected", 18),
            rpad("extra/curated", 15), "peaks/FF")
    row(name, s) = @printf("%-16s%-18s%-18s%-18s%-15.2f%.2f\n", name,
        @sprintf("%.2f (%s)", s.add_recall, f(s.adds)),
        @sprintf("%.2f (%s)", s.kept_recall, f(s.kept)),
        @sprintf("%.2f (%s)", s.excluded_still_detected, f(s.excluded)),
        s.extra_per_curated, s.peaks_per_form_factor)
    for (thr, s) in r.by_threshold
        row("p >= $thr", s)
    end
    row("hybrid", r.hybrid)
    row("all candidates", r.candidates)
end
