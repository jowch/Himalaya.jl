using Himalaya
using SparseArrays
using SQLite, DBInterface, Tables

"""
Tolerance (relative to predicted q) used when snapping observed peaks to
predicted ratio positions during speculative-index construction. Matches
`Himalaya.indexpeaks`'s default `tol` so the snap UI agrees with the
phase-search engine on what counts as "close enough".
"""
const SNAP_TOL = 0.0025

"""
Tolerance used when a CLIENT-FITTED custom index claims its peaks at commit.

Deliberately looser than `SNAP_TOL`: the custom-index modal draws its comb by
eye and reports "N of M reflections land" using `landsOn`'s `relTol` in
`frontend/src/lib/customIndex.ts`. Committing at `SNAP_TOL` (9× tighter) would
claim fewer peaks than the modal just showed the user — a second tolerance
regime disagreeing with the surface they fitted on. MUST equal that `relTol`.
"""
const CUSTOM_SNAP_TOL = 0.022

"""
    resolve_phase(name) -> Union{Type{<:Himalaya.Phase}, Nothing}

Resolve a phase name string (e.g. `"Pn3m"`, `"Himalaya.Pn3m"`) to the phase
type. Returns `nothing` if the name isn't a known `Himalaya.Phase` subtype.
"""
function resolve_phase(name::AbstractString)
    bare = last(split(String(name), '.'))
    P = try
        getfield(Himalaya, Symbol(bare))
    catch
        return nothing
    end
    (P isa Type && P <: Himalaya.Phase) ? P : nothing
end

"""
    compute_snap(peak_rows, phase, anchor_q, anchor_ratio; tol = SNAP_TOL)
        -> Vector{NamedTuple}

For each ratio position of `phase`, return the predicted q (basis × ratio)
and the closest non-anchor peak whose relative deviation is within `tol`.
`peak_rows` is any iterable of NamedTuples with `id::Int` and `q::Real`
fields (e.g. the result of `Tables.rowtable` over a `peaks` query).

`max_order` caps how many ratio positions are scanned (default: all). Callers
that must not claim beyond a shorter series the user was actually shown pass
their own length — see `insert_custom_index!`.

Returns one row per scanned ratio position (1-indexed) with fields:
- `ratio_position::Int`
- `predicted_q::Float64`
- `suggested_peak_id::Union{Int, Nothing}`
- `suggested_peak_kind::Union{String, Nothing}`
- `suggested_q::Union{Float64, Nothing}`
- `suggested_residual::Union{Float64, Nothing}`
- `suggested_relresid::Union{Float64, Nothing}`

`suggested_peak_kind` is carried from the row's own `peak_kind` field when the
caller supplies one (`_effective_peak_rows` does). `auto_peaks` and
`peak_curations` are independent AUTOINCREMENT namespaces, so a bare
`suggested_peak_id` does NOT identify a peak on its own — consumers that
persist the result must key on the (id, kind) PAIR. Rows without the field
(the /speculative-snap route builds its own) yield `nothing` and behave as before.
"""
function compute_snap(peak_rows, phase::Type{P}, anchor_q::Real, anchor_ratio::Int;
                     tol::Real = SNAP_TOL,
                     max_order::Union{Int, Nothing} = nothing) where {P<:Himalaya.Phase}
    ratios = Himalaya.phaseratios(P; normalize = true)
    1 <= anchor_ratio <= length(ratios) || error("anchor_ratio $anchor_ratio out of range for $P (1..$(length(ratios)))")
    basis = Float64(anchor_q) / ratios[anchor_ratio]
    last_order = max_order === nothing ? length(ratios) : min(max_order, length(ratios))

    out = NamedTuple[]
    for rpos in 1:last_order
        ratio = ratios[rpos]
        predicted_q = basis * ratio
        best_id      = nothing
        best_kind    = nothing
        best_q       = nothing
        best_resid   = nothing
        best_relresid = Inf
        for pr in peak_rows
            q  = Float64(pr.q)
            relresid = abs(q - predicted_q) / predicted_q
            if relresid <= tol && relresid < best_relresid
                best_relresid = relresid
                best_id    = Int(pr.id)
                best_kind  = hasproperty(pr, :peak_kind) ? String(pr.peak_kind) : nothing
                best_q     = q
                best_resid = abs(q - predicted_q)
            end
        end
        push!(out, (
            ratio_position      = rpos,
            predicted_q         = predicted_q,
            suggested_peak_id   = best_id,
            suggested_peak_kind = best_kind,
            suggested_q         = best_q,
            suggested_residual  = best_resid,
            suggested_relresid  = best_id === nothing ? nothing : best_relresid,
        ))
    end
    out
end

"""
    build_speculative_index(peak_rows, phase, ratio_to_peak_id::Dict{Int,Int})
        -> NamedTuple

Construct the numeric fields a speculative `indices` row needs, given a
ratio-position → peak-id assignment. Reuses `Himalaya.fit` and
`Himalaya.score` so speculative indices stay numerically consistent with
auto indices.

`peak_rows` must contain rows with `id`, `q`, and `sharpness` fields.

Returns `(; basis, score, r_squared, lattice_d, residuals)` where
`residuals` is `Dict{Int, Float64}` mapping ratio_position → |q - ideal|.
"""
function build_speculative_index(peak_rows, phase::Type{P},
                                  ratio_to_peak_id::Dict{Int,Int}) where {P<:Himalaya.Phase}
    isempty(ratio_to_peak_id) && error("at least one ratio assignment required")

    ratios = Himalaya.phaseratios(P; normalize = true)

    by_id = Dict{Int, eltype(peak_rows)}()
    for pr in peak_rows
        by_id[Int(pr.id)] = pr
    end

    rpos_sorted = sort(collect(keys(ratio_to_peak_id)))
    qvals      = Float64[]
    sharpvals  = Float64[]
    for rpos in rpos_sorted
        pid = ratio_to_peak_id[rpos]
        haskey(by_id, pid) || error("peak id $pid not in supplied peak_rows")
        pr = by_id[pid]
        push!(qvals,     Float64(pr.q))
        push!(sharpvals, pr.sharpness === nothing || ismissing(pr.sharpness) ? 0.0 : Float64(pr.sharpness))
    end

    n = length(ratios)
    peaks_sv     = SparseVector{Float64, Int}(n, rpos_sorted, qvals)
    sharpness_sv = SparseVector{Float64, Int}(n, rpos_sorted, sharpvals)

    # Least-squares fit through assigned (ratio, q) pairs (intercept fixed at
    # 0), against NORMALIZED ratios so the stored basis means "q of the first
    # ratio position" — the same convention auto indices use (core:
    # predictpeaks = basis × phaseratios(normalize=true)). Every consumer
    # (predicted_q_for_phase, MillerPlot) assumes this scale; fitting against
    # un-normalized ratios shrank cubic predictions by the first raw ratio.
    observed_ratios_used = ratios[rpos_sorted]
    basis = observed_ratios_used \ qvals

    # Index basis is only carried, never read, by fit/score (fit refits d
    # internally from peaks; score reads peaks + sharpness) — but pass the
    # normalized value so the constructed Index is convention-correct.
    idx = Himalaya.Index{P}(basis, peaks_sv, sharpness_sv)
    fit_result = Himalaya.fit(idx)
    s          = Himalaya.score(idx)

    # Residuals are numerically identical under either convention
    # (ratios_unnorm·basis_unnorm ≡ ratios_normed·basis_normed).
    residuals = Dict{Int, Float64}()
    for (rpos, qv) in zip(rpos_sorted, qvals)
        ideal = ratios[rpos] * basis
        residuals[rpos] = abs(qv - ideal)
    end

    (; basis = basis,
       score = s,
       r_squared = fit_result.R²,
       lattice_d = fit_result.d,
       residuals = residuals)
end

"""
    _kind_for(db, exposure_id, peak_id) -> String

Determine whether `peak_id` references an `auto_peaks` row or a
`peak_curations(kind='add')` row for the given exposure.
Returns `"auto"` or `"curation"`. Errors if not found in either table.
"""
function _kind_for(db::SQLite.DB, exposure_id::Int, peak_id::Int)
    if !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM auto_peaks WHERE id = ? AND exposure_id = ?",
        [peak_id, exposure_id])))
        return "auto"
    elseif !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM peak_curations WHERE id = ? AND exposure_id = ? AND kind = 'add'",
        [peak_id, exposure_id])))
        return "curation"
    else
        error("peak_id $peak_id not found in auto_peaks or peak_curations for exposure $exposure_id")
    end
end

"""
    lattice_d_for(phase, basis) -> Float64

The lattice parameter a phase's `basis` (the normalized q₁ slope) implies,
derived the way Himalaya does: build a 1-peak `Index` at the first predicted
reflection (q₁ == basis, since the normalized first ratio is 1.0) and read
`fit().d`. Phase-correct for cubic/lamellar/hex.

`Himalaya.fit` recomputes d from the PEAK values, never from `Index.basis`, so
this is the only way to ask "what lattice does this basis mean?" without a peak
set. Shared by `insert_custom_index!` (at commit) and the pipeline's
basis-locked branch (at reanalysis) so both answer identically — a locked index
whose `lattice_d` was derived one way and refreshed the other would drift in
the rail while `basis` sat still.
"""
function lattice_d_for(phase::Type{P}, basis::Float64) where {P<:Himalaya.Phase}
    n = length(Himalaya.phaseratios(P))
    peaks_sv     = SparseVector{Float64, Int}(n, [1], [basis])
    sharpness_sv = SparseVector{Float64, Int}(n, [1], [1.0])
    Himalaya.fit(Himalaya.Index{P}(basis, peaks_sv, sharpness_sv)).d
end

"""
    _effective_peak_rows(db, exposure_id) -> Vector{NamedTuple}

The exposure's effective peak set as `(id, q, sharpness, peak_kind)` rows:
non-excluded auto peaks + curation adds. Mirrors the pipeline's
`effective_peaks` view, in the shape `compute_snap` /
`build_speculative_index` expect.

`peak_kind` is selected as a literal per UNION arm because `auto_peaks.id` and
`peak_curations.id` are independent AUTOINCREMENT sequences (db.jl:174, :184):
the same integer can name a row in BOTH tables for one exposure, so an id
without its kind is ambiguous. Resolving the kind separately (e.g. via
`_kind_for`, which checks `auto_peaks` first) can disagree with the row that
actually won a snap — the collision `insert_speculative_index!` guards against
by re-querying the resolved table.
"""
_effective_peak_rows(db::SQLite.DB, exposure_id::Int) =
    Tables.rowtable(DBInterface.execute(db, """
        SELECT a.id, a.q, a.sharpness, 'auto' AS peak_kind
        FROM auto_peaks a
        WHERE a.exposure_id = ?
          AND NOT EXISTS (
              SELECT 1 FROM peak_curations c
              WHERE c.exposure_id = a.exposure_id AND c.kind = 'exclude'
                AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
          )
        UNION ALL
        SELECT id, q, NULL AS sharpness, 'curation' AS peak_kind
        FROM peak_curations
        WHERE exposure_id = ? AND kind = 'add'
    """, [exposure_id, exposure_id]))

"""
    insert_speculative_index!(db, exposure_id, phase, ratio_to_peak_id)
        -> Int (new index id)

Builds the speculative index from supplied peak assignments and inserts
it into `indices` with `kind='speculative'`, plus the per-peak rows in
`index_peaks`. Returns the new id. Does **not** add the index to any
group — the caller is responsible for active-group membership.
"""
function insert_speculative_index!(db::SQLite.DB, exposure_id::Int,
                                    phase::Type{P},
                                    ratio_to_peak_id::Dict{Int,Int}) where {P<:Himalaya.Phase}
    peak_rows = _effective_peak_rows(db, exposure_id)

    built = build_speculative_index(peak_rows, P, ratio_to_peak_id)

    # Inherit the exposure's current analysis_inputs_hash so this new index
    # doesn't read as "stale" the moment it lands. The speculative is built
    # from the same effective peak set the exposure hash already covers, so
    # any inputs_hash other than `analysis_inputs_hash` would be misleading.
    # Without this, the stale-indices alert fires immediately after every
    # speculative create (NULL ≠ exposure hash → mismatch).
    current_hash = read_inputs_hash(db, exposure_id)

    res = DBInterface.execute(db,
        """INSERT INTO indices
             (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind, inputs_hash)
           VALUES (?, ?, ?, ?, ?, ?, 'candidate', 'speculative', ?)""",
        [exposure_id, string(nameof(P)), built.basis,
         built.score, built.r_squared, built.lattice_d, current_hash])
    new_id = Int(DBInterface.lastrowid(res))

    for (rpos, peak_id) in ratio_to_peak_id
        pk_kind = _kind_for(db, exposure_id, peak_id)
        # The intent q must come from the same table `_kind_for` resolved to:
        # auto_peaks and peak_curations ids live in independent AUTOINCREMENT
        # namespaces, so a bare-id dict over the UNION'd peak_rows would let a
        # colliding curation row's q shadow an auto peak's (curation rows come
        # last) — freezing the wrong q into durable intent state.
        q_row = pk_kind == "auto" ?
            Tables.rowtable(DBInterface.execute(db,
                "SELECT q FROM auto_peaks WHERE id = ?", [peak_id])) :
            Tables.rowtable(DBInterface.execute(db,
                "SELECT q FROM peak_curations WHERE id = ?", [peak_id]))
        DBInterface.execute(db,
            """INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
               VALUES (?, ?, ?, ?, ?)""",
            [new_id, peak_id, pk_kind, rpos, built.residuals[rpos]])
        # Durable intent: index_peaks is the per-analysis resolved view (wiped
        # and rebuilt by _persist_analysis_inner!); this row is the user's
        # assignment and survives every wipe. Frozen at creation.
        DBInterface.execute(db,
            """INSERT INTO speculative_peak_intents (index_id, ratio_position, q)
               VALUES (?, ?, ?)""",
            [new_id, rpos, Float64(q_row[1].q)])
    end
    new_id
end

"""
    insert_custom_index!(db, exposure_id, phase, basis) -> Int (new index id)

Plan D Task D-9 (B4): persist a CLIENT-FITTED custom index. The frontend's
custom-index modal computes `basis` directly from a symmetry + a user-chosen
lattice via real physics (`2π√N/a` cubic, `2πn/d` lamellar, `4π√M/(√3·a)` hex),
so the backend does NOT re-fit against observed peaks — it stores the supplied
basis verbatim and derives `lattice_d` via `Himalaya.fit` on a synthetic first
reflection, so `predicted_q_for_phase(phase, basis)` reproduces the modal's
client-side comb exactly.

CONVENTION (review finding #4): the supplied `basis` is the q₁ slope —
`2π/a × first(phaseratios(P))` — i.e. the value such that
`predicted_q_for_phase` (which multiplies by the NORMALIZED ratios, whose first
entry is 1.0) reproduces the physical first reflection. It is NOT `a` and NOT
`2π/a`. The √6 first reflection of Ia3d maximizes the convention-mismatch
signal, so the round-trip contract test pins it specifically.

The stored `basis`/`lattice_d` are the user's lattice VERBATIM — not refit
through the observed peaks (unlike `insert_speculative_index!`, which solves for
basis), and `basis_locked = 1` makes that DURABLE: `_persist_analysis_inner!`
re-resolves a locked index's peaks like any other speculative but skips the
least-squares refit, so no reanalysis can move the lattice the user chose. See
`docs/event-log.md` and `src/AGENTS.md` on why locking is what makes the
scan-derived intents below legitimate.

What IS written is the peak assignment the modal already showed the user: every
ratio position whose predicted q has an observed peak within `CUSTOM_SNAP_TOL`
claims that peak, in `index_peaks` (the per-analysis resolved view) AND
`speculative_peak_intents` (the durable record that survives the reanalysis
wipe). Without these rows the Focus comb, detector rings and assignment cart all
render a fully-fitted custom index as claiming nothing.

`orders` bounds the scan to the number of reflections the MODAL DREW. The
frontend's `SYMS[sym].Ms` (`lib/customIndex.ts`) is shorter than the core ratio
series for five of eight phases (Pn3m 6 vs 16, Im3m 6 vs 10, Ia3d 6 vs 8,
Lamellar 5 vs 11, Hexagonal 6 vs 14); scanning the full series would claim
reflections the user was never shown and let the rail's "explains N peaks"
exceed the modal's "N of M land" for the same fit. `nothing` scans the whole
series — correct only for a caller with no truncated display of its own.

`score`/`r_squared` stay NULL at commit. They are populated on the first
reanalysis that resolves any intent (pipeline.jl), computed against the LOCKED
basis. Returns the new index id.
"""
function insert_custom_index!(db::SQLite.DB, exposure_id::Int,
                              phase::Type{P}, basis::Float64;
                              orders::Union{Int, Nothing} = nothing) where {P<:Himalaya.Phase}
    basis > 0 || error("basis must be positive")
    orders === nothing || orders > 0 || error("orders must be positive")

    # Shared with the pipeline's basis-locked branch so a reanalysis refreshes
    # lattice_d to exactly the value committed here.
    lattice_d = lattice_d_for(P, basis)

    current_hash = read_inputs_hash(db, exposure_id)
    res = DBInterface.execute(db,
        """INSERT INTO indices
             (exposure_id, phase, basis, score, r_squared, lattice_d, status, kind,
              inputs_hash, basis_locked)
           VALUES (?, ?, ?, NULL, NULL, ?, 'candidate', 'speculative', ?, 1)""",
        [exposure_id, string(nameof(P)), basis, lattice_d, current_hash])
    new_id = Int(DBInterface.lastrowid(res))

    # Claim the peaks the modal showed landing. anchor_ratio = 1 with
    # anchor_q = basis reproduces the stored comb exactly (normalized ratios[1]
    # is 1.0), so compute_snap matches against predicted q's, not a refit.
    # `max_order` bounds the scan to the orders the modal actually DREW (see
    # the `orders` note above) — without it the backend would claim reflections
    # the user was never shown for every phase whose SYMS.Ms is truncated.
    snaps = compute_snap(_effective_peak_rows(db, exposure_id), P, basis, 1;
                         tol = CUSTOM_SNAP_TOL, max_order = orders)

    # Best fit wins a contested peak. One peak may sit inside CUSTOM_SNAP_TOL of
    # two orders less than 2·tol apart (Pn3m √16/√17 are 3.1% apart, √19/√20
    # 2.6%, Fd3m √35/√36 1.4%); claiming it for whichever order came first would
    # label it with the wrong Miller index. Sorting by relative residual makes
    # the closer order win, and the dedup below drops the looser one.
    hits = sort(filter(s -> s.suggested_peak_id !== nothing, snaps);
                by = s -> s.suggested_relresid)

    # One claim per PEAK, keyed on the (id, kind) PAIR: `suggested_peak_id`
    # alone is ambiguous across the two AUTOINCREMENT namespaces (see
    # `_effective_peak_rows`), so keying on the bare id would treat auto 7 and
    # curation 7 as one peak and silently drop a genuinely distinct claim.
    # Load-bearing beyond correctness: index_peaks' PK is
    # (index_id, peak_id, peak_kind) (db.jl:170) and this is a plain INSERT, so
    # a duplicate would RAISE, not be ignored — unlike the pipeline's call
    # sites, which use INSERT OR IGNORE because they tolerate re-resolution
    # collisions. Here a collision is a bug, so it should never be reached.
    claimed = Set{Tuple{Int, String}}()
    for s in hits
        # peak_kind comes from the row that actually won the snap, never from a
        # second lookup — and the intent q comes from that same row, so
        # index_peaks and speculative_peak_intents can never describe different
        # peaks (the durable-state hazard insert_speculative_index! documents).
        key = (s.suggested_peak_id, s.suggested_peak_kind)
        key in claimed && continue
        push!(claimed, key)
        DBInterface.execute(db,
            """INSERT INTO index_peaks
                 (index_id, peak_id, peak_kind, ratio_position, residual)
               VALUES (?, ?, ?, ?, ?)""",
            [new_id, s.suggested_peak_id, s.suggested_peak_kind,
             s.ratio_position, s.suggested_residual])
        DBInterface.execute(db,
            """INSERT INTO speculative_peak_intents (index_id, ratio_position, q)
               VALUES (?, ?, ?)""",
            [new_id, s.ratio_position, s.suggested_q])
    end

    new_id
end
