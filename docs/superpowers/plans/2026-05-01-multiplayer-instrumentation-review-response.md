# Multiplayer + Instrumentation — Response to Review

**Date:** 2026-05-01
**Reviewing:** [2026-05-01-multiplayer-instrumentation-design.md](../specs/2026-05-01-multiplayer-instrumentation-design.md), [2026-05-01-multiplayer-instrumentation.md](2026-05-01-multiplayer-instrumentation.md)

This is the structured rebuttal to the review of the spec and plan. Disposition tags: **Agreed — fixed** (substantive issue, change applied), **Agreed — noted** (small issue, change applied without further discussion), **Partial** (split agreement; clarification below), **Acknowledged** (reviewer was reading something I'd already meant to scrub or the issue is real but out of scope here).

Where I pushed back, the reasoning is in-line — not a hedge.

---

## Load-bearing issues

**1. `apply_event!` contract contradicts itself — Agreed — fixed.**

Real bug. The plan as written had routes INSERT into `peak_curations` *and* call `apply_event!`, with the dispatcher returning early. That made the property test "rebuild folds to incremental state" a tautology — the rebuild can't recreate a row that was never written by the dispatcher in the first place. Fixed in two places: the spec's §"Materialized view contract" now states that the dispatcher is the *sole* writer to view tables, and the plan's R4.2 example route is rewritten so it only validates input and calls `apply_event!`; the dispatcher branch performs the `peak_curations` INSERT. The property test is now a contract enforcer — bypassing `apply_event!` to write a view directly will fail the rebuild test, which is the whole point of having it.

**2. Negative-id sentinel + `peak_kind` column is double-encoding — Agreed — fixed.**

Right. The negative-id trick was a hangover from when I was thinking about the layers separately and then merged them without revisiting. Dropped from both spec and plan. `effective_peaks` now returns `peak_id::Vector{Int}` and `peak_kind::Vector{Symbol}` as parallel arrays; `index_peaks.peak_kind` is the only on-disk discriminator. Spec includes a one-line note explaining why an earlier draft used the sentinel, in case future readers wonder.

**3. R2 silently breaks index_peaks for speculatives anchored on manual peaks — Agreed — fixed (highest-impact catch in the review).**

This was a genuine correctness bug — silent data loss for users with hand-built speculatives. The original migration deleted `index_peaks` rows whose `peak_id` came from manual peaks and counted on reanalysis recovery, but recovery requires the speculative to score well under the new effective peaks; if it didn't, the user's hand-built indexing vanished without warning. The fix in R2.1 is structural: instead of `INSERT INTO peak_curations ... SELECT ...`, the migration iterates row-by-row to capture each `(old_peak_id, new_curation_id)` pair, then issues `UPDATE index_peaks SET peak_id = ?, peak_kind = 'curation'` for each. Followed by an orphan-detection check that fails loudly if any `index_peaks` row references a now-defunct id. New regression test added: "migrate_r2_split_peaks! preserves index_peaks for speculatives anchored on manual peaks."

**4. Sharpness on `peak_curations` is a stale-cache bug waiting to happen — Agreed — fixed.**

The reviewer's framing is exactly right: sharpness depends on the trace, so persisting it on the curation row decouples it from the trace and makes the R3 hash lie. I copied the column shape from the existing `peaks(source='manual').sharpness` without thinking about what it meant under the memoization model. Dropped the column from `peak_curations`. `effective_peaks` now takes `(q_grid, I)` as args and samples `Himalaya.sharpness(I)` for each `add` curation on every call. The `backfill_curation_sharpness!` helper in the plan was deleted entirely (it was solving a problem that no longer exists). Inline note in the spec schema explains the choice for future readers.

---

## SSE-specific issues (R5)

**5. Busy-poll subscriber loop with 1s `sleep` — Agreed — fixed.**

Replaced with a `Threads.Condition` that `broadcast_event!` notifies, plus a per-subscriber `Channel` for pending frames. Subscriber loop now blocks on the condition + a 15s timer (for heartbeats) instead of polling. Event latency is now bounded by network RTT instead of 1s polling interval; idle subscribers don't pin worker threads.

**6. SSE behind nginx will buffer to death — Agreed — fixed.**

`X-Accel-Buffering: no` header added; 15s heartbeat comment frames (`:heartbeat\n\n`) emit when no events arrive. Both the spec and plan now explicitly call out the reverse-proxy compatibility because the lab's `/opt/himalaya` deployment model implies nginx termination.

**7. Broadcast happens after commit, outside the transaction — Agreed — noted.**

Real concern but the right resolution is "best-effort delivery + reconcile on reconnect," not "make broadcast transactional." Spec now spells this out in §R5a Server side: SSE is a notification channel; `user_actions` is the durable log; if the process dies between commit and broadcast, the event survives but the frame is lost — clients refetch on reconnect via TanStack Query. Connected this explicitly to Open Question 4 (`Last-Event-ID` is the principled v2 fix and would also cover crash-during-broadcast). Did not promote `Last-Event-ID` to v1: it requires server-side replay logic that's not justified before R4 data shows reconnect-gap losses are happening at meaningful rates.

**8. `Oxygen.jl @stream` may not exist — Agreed — fixed.**

You're right that I was hedging. Looked at it: Oxygen 1.10.x's macro surface for streaming responses isn't `@stream` — the actual pattern is closer to `@get` with manual stream access via the underlying HTTP.jl primitives. Replaced the plan's R5.1 step 1 with a new R5.0 task: a 30-line scratch spike before any R5a code commits. The spike's deliverable is a documented decision: either Oxygen composes cleanly with SSE, or the fallback is a separate `HTTP.serve` for the event channel on a different port (sharing the same DB Ref). Downstream R5a tasks are gated on the spike resolving.

---

## Smaller things

**9. R2 sets `status='stale'`, R3 drops the value — Agreed — fixed.**

Added to R3.2 step 3: `UPDATE indices SET status = 'candidate' WHERE status = 'stale'` after the hash columns are added. The `status` column itself stays (`'candidate'` is still meaningful and may grow other values later); just the orphaned `'stale'` value is normalized away.

**10. Self-echo for shared username — Agreed — fixed.**

Added the comment in both the spec example and the plan's R5.3 step 1. Rationale documented (lab edge case, not worth per-tab client ids unless it actually bites someone).

**11. Time estimates are optimistic — Acknowledged.**

Earlier in the design conversation the user explicitly asked me to drop time estimates. I removed them from the recommendation block but missed the "Days" / "Days–weeks" tags in the migration-plan list at the bottom of the spec. Now scrubbed. To the reviewer's substantive point: yes, the estimates were optimistic, especially R2 and R4. Without estimates the plan is honest that scope-not-time is what's bounded.

**12. Scope question: SSE-first, defer If-Match until data justifies it — Partial agreement.**

**Strongly agreed on the SSE/If-Match split.** Restructured R5 into R5a (SSE broadcast + LWW) and R5b (If-Match + 409 retry). R5b is now explicitly gated on R4 instrumentation showing concurrent edits to the same exposure within a hash window happening at meaningful rates. If R4 data shows contention is rare, R5b can be skipped indefinitely — LWW under R5a is the resolution model. This is exactly the data-driven phasing the spec's own framing argues for ("the instrumentation we just built tells us whether the next layer is needed"), and getting it consistent with the spec's principles is a real improvement.

**Partial on the R2 deferral.** The reviewer also suggested "ship R0+R1+R3 first, decide on R2 vs. just-add-`created_by` once R1 is live." Two reasons I'm not committing to that order in the plan, though I added the fallback to Open Questions:

- **R1 alone doesn't fix the snapshot/restore dance.** With auto peak ids preserved (R1) but manual peaks still in the same table, `_persist_analysis_inner!` continues to snapshot excluded q-values, fold manual peaks into the indexpeaks input, and re-resolve speculative `index_peaks` references. The pipeline-clarity argument the spec makes — "180 lines of workaround go away" — only fires after R2. Deferring R2 means deferring the simplification that the spec is partly built around.
- **R4's payloads reference curation rows as first-class entities.** With `peaks` unified, R4 has to encode curation events as `(peak_id, excluded=1)` deltas — workable but harder to fold and harder to query historically. The argument "ship R4 instrumentation first" is real, but R4-after-R2 produces a more useful event log.

That said, the reviewer's escape hatch is genuinely real: if R2 implementation reveals unexpected complexity, dropping back to `peaks.created_by` + R1 + R5a is a viable shipping plan that delivers ~80% of the multiplayer story without R2's cost. Surfaced as Open Question 7 in the spec with concrete trade-offs spelled out, rather than buried as an implementation hedge.

---

## Items called out as good — kept

The four items the reviewer flagged as well-done (idempotency contract, R0.1 pre-flight check, R2.3 split-out, R4 payload-discipline table) are unchanged. Specifically the idempotency contract / DB state matrix and R2.3 separation reflect deliberate design choices the reviewer correctly identified — those stay because they are load-bearing for safe migration.

---

## Closing — what's different after this review

Net changes:

- One real correctness bug fixed (R2 silently dropping speculatives).
- One real contract bug fixed (`apply_event!` not actually enforcing what the spec implied).
- Two real architectural improvements (sharpness re-derivation, R5a/R5b split).
- Five real operational improvements (SSE concurrency primitive, nginx headers, Oxygen spike, orphan stale cleanup, shared-username comment).
- Two minor doc cleanups (time estimates scrubbed, negative-id sentinel removed).
- One open question opened (R2 vs. `created_by` fallback path).

The review caught the kind of issues that don't show up until implementation hits them — the speculative-breakage in particular would have shipped silent data loss. Worth the cost of writing the spec and plan deliberately enough to be reviewable.

Open work the review surfaced but didn't resolve:
- Oxygen 1.10.x streaming API spike (R5.0 task — must complete before R5a code commits).
- R5b decision is now data-driven on R4 — explicit gating is in the plan, but the actual go/no-go happens after R0–R5a have been live in production for a sustained-use period.

---

## Second-round review — response

The reviewer's second pass caught real cross-section drift: most fixes from round 1 landed locally, but their callers in adjacent sections weren't updated. Disposition tags as before.

**1. R3.2 calls `effective_peaks` with wrong arity — Agreed — fixed.** Real bug. Round 1 added `(q, I)` to `effective_peaks`'s signature in R2.2 but missed the call site in R3.2's hash-guarded `analyze_exposure!`. Fixed.

**2. R3.2 calls deleted `backfill_curation_sharpness!` — Agreed — fixed.** The rebuttal claimed deletion; the call survived in R3.2 step 2 because the snippet pre-dated the sharpness decoupling and wasn't re-read. Fixed: call deleted; sharpness derivation lives entirely in `effective_peaks` per the round 1 design.

**3. R2.2 step 4 GET handler SELECTs phantom `sharpness` column — Agreed — fixed.** Same drift class. Schema correctly omits the column, but the UNION query was still selecting it. Replaced with `NULL AS sharpness`.

**4. `synthesize_peaks_result` referenced never defined; argmin reconstruction is lossy — Agreed — fixed (highest-impact second-round catch).** The reviewer correctly noted argmin-nearest doesn't equal local-maximum-sample, which subtly wrong-answers downstream `intensity` lookups via `I[indices]`. Fix is structural rather than rhetorical: added `findpeaks_index INTEGER` to the `auto_peaks` schema, so `diff_update_auto_peaks!` can persist the exact grid index from `peaks_result.indices[i]` at INSERT/UPDATE time. `synthesize_peaks_result` now reads it directly; argmin is the fallback only for legacy rows from R2 migration where the column is NULL (those self-heal on next analyze). Helper is now spelled out, not hand-waved. R2.2 step 2 now explicitly notes that `diff_update_auto_peaks!` (introduced in R1 against legacy `peaks`) needs an R2 update to write the new column and target the new tables.

**5. R2 → R3 deployment gap leaves indices stale forever — Agreed — fixed.** Real concern: shipping R2 alone deletes the `mark_all_indices_stale!` call sites *and* sets `status='stale'` on every index, but R3's hash-mismatch banner isn't yet wired up to clear them — UI stuck in stale state. Added an explicit DEPLOYMENT NOTE at the migration's `UPDATE indices SET status = 'stale'` line: ship R2 and R3 together (single release / one `migrate_schema!` invocation). The reviewer's framing — "either bundle R2+R3 or add a transitional clear step" — pointed to the cleaner choice; bundling avoids a transitional-clear branch that would itself need testing and removal.

**6. R4.2 dispatcher migration non-atomic per kind — Agreed — fixed.** The plan said "ship `apply_event!` alongside `log_action!`, migrate call sites in batches" without making the per-kind atomicity explicit. A route flipped before its dispatcher branch lands silently drops the view-side write. Added an explicit "Atomicity requirement (per kind)" paragraph to R4.2's migration impact: the route refactor and matching dispatcher branch land in the same commit, never across two passes. Each route's commit now bundles route + dispatcher branch + property test together — this is what the rebuild-views-from-log property test is *supposed to* enforce, but only after the kind is covered, so per-kind landing is the load-bearing rule.

**7. R5.1 SSE loop has TOCTOU between `isready` and `wait` — Agreed — fixed (cleaner pattern adopted).** The reviewer's preferred design (per-subscriber `Channel` + `Timer` that `put!`s heartbeat frames; loop is just `take!(pending)`) is genuinely better — race-free *and* simpler. Replaced the hybrid Condition+isready+wait sketch with the Channel-only pattern. The shared `SSE_WAKEUP` Condition is now noted as superseded; R5.0 spike's job is to confirm the Channel pattern composes with whichever Oxygen streaming primitive it picks (and fall back to Condition if not, but with eyes-open about the race). Broadcast helper updated to drop the now-redundant `notify(SSE_WAKEUP)` call.

**8. R2.1 step numbering bug — Agreed — fixed.** Two "Step 3"s, then 4, 5. Renumbered to 3-4-5-6. An executing agent following checkboxes won't lose its place.

**9. `lookup_username` undefined — Agreed — fixed.** Trivial helper added inline next to `broadcast_event!` so the snippet is self-contained.

**10. `current_db()` "is not the project's pattern" — Disagreed.** Reviewer was reading CLAUDE.md too literally. The Ref-based DB state and the `current_db()` accessor coexist: the Ref is `_DB_REF` at [server.jl:7](packages/HimalayaUI/src/server.jl:7), and `current_db()` at [server.jl:9](packages/HimalayaUI/src/server.jl:9) is the project-wide accessor — every existing route file uses it (`grep current_db packages/HimalayaUI/src/routes_*.jl` shows ~30 hits). The plan's usage matches the convention; no change. CLAUDE.md says "a Ref holds the live DB" because that's the storage; `current_db()` is the read interface.

**11. R3.2 step 1 doesn't say where to add columns to fresh `SCHEMA` — Agreed — fixed.** Clarified that the new columns are *appended to the existing CREATE TABLE statements* in `SCHEMA`, not added as a separate ALTER block. The `migrate_schema!` ALTER calls handle existing DBs; `SCHEMA` describes the post-migration state for fresh DBs (where the ALTERs are no-ops because the columns already exist via CREATE).

**12. `persist_analysis!` rewrite is hand-waved — Agreed — fixed.** Replaced the "reference the current implementation closely" punt with a structural sketch: GONE block (snapshot/restore/merge — ~180 lines), KEEPS block (groups, speculatives, status), and a small `eff_lookup` Dict that replaces the `q_to_peak_id` lookup. Reader can see what shrinks vs. survives without reverse-engineering from the surrounding tasks.

---

## What's different after the second round

Net additional changes:
- One real correctness bug fixed (R3.2 calling helpers with wrong signatures / deleted helpers / phantom columns).
- One real improvement to a fragile reconstruction (`findpeaks_index` persisted, argmin demoted to legacy fallback).
- One deployment correctness fix (R2+R3 bundling note).
- One contract clarification (R4.2 per-kind atomicity).
- One operational improvement (Channel-only SSE pattern, race-free).
- Three doc/structural cleanups (step numbering, lookup_username inline, persist_analysis! sketch).
- One reviewer disagreement (current_db is the actual project pattern).

These are the bugs that escape because each section was edited in isolation. The takeaway is that round 1's verification ("did the fix land?") needs to be paired with a cross-section reading ("does the fix still match what other sections now expect?") before each handoff. None of the bugs are deep, but every one would have stopped an executing agent at first run. Worth the cost of a second pass.
