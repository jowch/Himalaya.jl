/**
 * Regression tests: mutator `onSuccess` paths must remain correct when SSE
 * beats HTTP and the response handed in is `synthesizeResponseFromSse`'s
 * minimal/synthetic shape rather than the route's full HTTP body.
 *
 * Fixed bugs covered here:
 *   - createSpeculative: SSE-wins gave `{event_id, client_op_id, index_id}`
 *     with no `id`/`phase`. Old code spliced this directly, landing a phantom
 *     IndexEntry with undefined fields. Fixed by guarding splice on `phase`
 *     presence and falling back to invalidate.
 *   - updateSample: SSE payload is the diff (only patched fields). Old code
 *     destructively spread `name/notes/label` from response, clobbering
 *     unpatched fields to undefined. Fixed by skipping undefined fields.
 *   - addSampleTag / addExposureTag: SSE uses `tag_id`; HTTP uses `id`.
 *     Old synth produced a tag with `id: undefined, source: undefined`.
 *     Fixed by adding kind-aware synthesis for `add_tag`.
 *   - peakSetExcluded (peak_excluded / peak_unexcluded): SSE payload is
 *     `{q, auto_peak_id}` — no `id`. Old fallback produced a response with
 *     `id: undefined`, so onSuccess's `pk.id === peakOnly.id` matched no row
 *     and the canonical state never replaced the optimistic one. Fixed by
 *     adding kind-aware synthesis mapping `auto_peak_id → id` and switching
 *     onSuccess to merge fields onto the existing row (preserves
 *     intensity/prominence/sharpness which the SSE payload omits).
 *   - addIndexToGroup / removeIndexFromGroup: SSE payload is `{group_id,
 *     index_id}` — no `id`. Old code spliced via `g.id === row.id` which
 *     reduced to `g.id === undefined` and matched nothing. Worse, the
 *     first confirmation of an exposure creates a NEW custom group via
 *     `ensure_custom_group!` and demotes the auto group, so even HTTP-wins
 *     missed the cache (cache had auto id=1; response had custom id=5).
 *     Fixed by detecting id-mismatch in onSuccess and invalidating the
 *     groups query so the next read fetches the canonical structure.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { updateSampleMutator } from "../../src/lib/queue/mutators/trivial";
import {
  peakExcludeMutator,
  peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator,
  removeIndexFromGroupMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";
import { queryKeys } from "../../src/queries";

describe("mutator onSuccess on SSE-wins synthetic responses", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient();
  });

  it("createSpeculative.onSuccess invalidates when response lacks IndexEntry shape (SSE-wins)", () => {
    qc.setQueryData(queryKeys.indices(42), [{ id: 1, phase: "Pn3m" } as any]);
    let invalidated = false;
    const origInvalidate = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 42 && k[2] === "indices") {
        invalidated = true;
      }
      return origInvalidate(filters);
    };
    // Synthesized SSE-wins shape: only event metadata + payload
    // (`{index_id}`). NO `id`/`phase`/`score`/etc.
    const sseSynth = {
      event_id: 7,
      client_op_id: "op-1",
      analysis_inputs_hash: undefined,
      index_id: 99,
    } as any;
    createSpeculativeMutator.onSuccess(
      { exposureId: 42, phase: "Pn3m", anchor_peak_id: 1, anchor_ratio: 1.0,
        additional: [], username: "u", clientId: "c", clientOpId: "op-1" } as any,
      sseSynth,
      qc,
    );
    expect(invalidated).toBe(true);
    // Cache was NOT polluted with a phantom row
    const indices = qc.getQueryData<any[]>(queryKeys.indices(42));
    expect(indices).toEqual([{ id: 1, phase: "Pn3m" }]);
  });

  it("createSpeculative.onSuccess splices when response has full IndexEntry shape (HTTP-wins)", () => {
    qc.setQueryData(queryKeys.indices(42), [{ id: 1, phase: "Pn3m" } as any]);
    const httpResponse = {
      id: 99, exposure_id: 42, phase: "Im3m", basis: 1.0, score: 0.9,
      r_squared: 0.99, lattice_d: 100, status: "speculative", kind: "speculative",
      inputs_hash: "h", peaks: [], predicted_q: [], ngc: -2.0,
    } as any;
    createSpeculativeMutator.onSuccess(
      { exposureId: 42, phase: "Im3m", anchor_peak_id: 1, anchor_ratio: 1.0,
        additional: [], username: "u", clientId: "c", clientOpId: "op-1" } as any,
      httpResponse,
      qc,
    );
    const indices = qc.getQueryData<any[]>(queryKeys.indices(42));
    expect(indices).toHaveLength(2);
    expect(indices![1]).toMatchObject({ id: 99, phase: "Im3m" });
  });

  it("peakExclude.onSuccess writes canonical state from synth (SSE-wins, id mapped from auto_peak_id)", () => {
    // Cache row reflects a stale state (excluded=false) — exercises the
    // canonical-replaces-optimistic contract. Pre-fix, synth lacks `id`,
    // onSuccess's match silently no-ops, and the stale state persists.
    qc.setQueryData(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
        sharpness: 30.0, source: "auto", excluded: false },
    ]);
    qc.setQueryData(queryKeys.exposure(5), {
      id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
      status: null, image_path: null, image_version: "",
      trace_hash: null, analysis_inputs_hash: "h0",
      tags: [], sources: [],
    });
    // Post-fix synthesizeResponseFromSse for peak_excluded:
    // {event_id, client_op_id, analysis_inputs_hash, id, q, source, excluded}.
    // intensity/prominence/sharpness intentionally absent — the SSE frame
    // doesn't carry them; onSuccess must merge, not replace.
    const sseSynth = {
      event_id: 7,
      client_op_id: "op-pe",
      analysis_inputs_hash: "h1",
      id: 7,
      q: 0.5,
      source: "auto",
      excluded: true,
    } as any;
    peakExcludeMutator.onSuccess(
      { exposureId: 5, peakId: 7, q: 0.5,
        username: "u", clientId: "c", clientOpId: "op-pe" } as any,
      sseSynth,
      qc,
    );
    const peaks = qc.getQueryData<any[]>(queryKeys.peaks(5))!;
    expect(peaks).toHaveLength(1);
    expect(peaks[0]).toEqual({
      id: 7, exposure_id: 5, q: 0.5,
      // Preserved from optimistic state — synth doesn't carry these:
      intensity: 1.2, prominence: 0.8, sharpness: 30.0,
      source: "auto",
      // Canonical from synth:
      excluded: true,
    });
    const exp = qc.getQueryData<any>(queryKeys.exposure(5))!;
    expect(exp.analysis_inputs_hash).toBe("h1");
  });

  it("addIndexToGroup.onSuccess invalidates when response id is not in cache (issue #37 Bug 1a, HTTP-wins, new custom group)", () => {
    // Cache only has the auto group (id=1). The first confirmation triggers
    // `ensure_custom_group!` which creates a fresh custom group (id=5) and
    // demotes the auto group. The HTTP response is the new custom group.
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    let invalidated = false;
    const origInvalidate = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return origInvalidate(filters);
    };
    // HTTP response shape (routes_analysis.jl POST /groups/:id/members).
    // The response id (5) does NOT exist in the cache.
    const httpResponse = {
      id: 5, exposure_id: 5, kind: "custom", active: true, members: [42],
      event_id: 11, view_row_id: 5,
    } as any;
    addIndexToGroupMutator.onSuccess(
      { exposureId: 5, groupId: 1, indexId: 42,
        username: "u", clientId: "c", clientOpId: "op-c1" } as any,
      httpResponse,
      qc,
    );
    expect(invalidated).toBe(true);
  });

  it("addIndexToGroup.onSuccess invalidates when synth lacks id (issue #37 Bug 1b, SSE-wins)", () => {
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    let invalidated = false;
    const origInvalidate = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return origInvalidate(filters);
    };
    // Synth output for index_confirmed when SSE wins: only event metadata
    // and the {group_id, index_id} payload. No `id`.
    const sseSynth = {
      event_id: 7,
      client_op_id: "op-c2",
      analysis_inputs_hash: undefined,
      group_id: 5,
      index_id: 42,
    } as any;
    addIndexToGroupMutator.onSuccess(
      { exposureId: 5, groupId: 1, indexId: 42,
        username: "u", clientId: "c", clientOpId: "op-c2" } as any,
      sseSynth,
      qc,
    );
    expect(invalidated).toBe(true);
  });

  it("addIndexToGroup.onSuccess splices when response id matches existing cache row (HTTP-wins, existing custom group)", () => {
    // Subsequent confirmations: custom group already exists in cache.
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: false, members: [] },
      { id: 5, exposure_id: 5, kind: "custom", active: true, members: [42] },
    ]);
    let invalidated = false;
    const origInvalidate = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return origInvalidate(filters);
    };
    const httpResponse = {
      id: 5, exposure_id: 5, kind: "custom", active: true, members: [42, 99],
      event_id: 12, view_row_id: 5,
    } as any;
    addIndexToGroupMutator.onSuccess(
      { exposureId: 5, groupId: 5, indexId: 99,
        username: "u", clientId: "c", clientOpId: "op-c3" } as any,
      httpResponse,
      qc,
    );
    expect(invalidated).toBe(false);
    const groups = qc.getQueryData<any[]>(queryKeys.groups(5))!;
    const custom = groups.find((g) => g.id === 5)!;
    expect(custom.members).toEqual([42, 99]);
  });

  it("removeIndexFromGroup.onSuccess invalidates when synth lacks id (SSE-wins)", () => {
    qc.setQueryData(queryKeys.groups(5), [
      { id: 5, exposure_id: 5, kind: "custom", active: true, members: [42] },
    ]);
    let invalidated = false;
    const origInvalidate = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return origInvalidate(filters);
    };
    const sseSynth = {
      event_id: 8,
      client_op_id: "op-u1",
      analysis_inputs_hash: undefined,
      group_id: 5,
      index_id: 42,
    } as any;
    removeIndexFromGroupMutator.onSuccess(
      { exposureId: 5, groupId: 5, indexId: 42,
        username: "u", clientId: "c", clientOpId: "op-u1" } as any,
      sseSynth,
      qc,
    );
    expect(invalidated).toBe(true);
  });

  it("peakUnexclude.onSuccess writes canonical state from synth (SSE-wins)", () => {
    qc.setQueryData(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
        sharpness: 30.0, source: "auto", excluded: true },
    ]);
    qc.setQueryData(queryKeys.exposure(5), {
      id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
      status: null, image_path: null, image_version: "",
      trace_hash: null, analysis_inputs_hash: "h0",
      tags: [], sources: [],
    });
    const sseSynth = {
      event_id: 8,
      client_op_id: "op-pue",
      analysis_inputs_hash: "h2",
      id: 7,
      q: 0.5,
      source: "auto",
      excluded: false,
    } as any;
    peakUnexcludeMutator.onSuccess(
      { exposureId: 5, peakId: 7, q: 0.5,
        username: "u", clientId: "c", clientOpId: "op-pue" } as any,
      sseSynth,
      qc,
    );
    const peaks = qc.getQueryData<any[]>(queryKeys.peaks(5))!;
    expect(peaks[0]).toEqual({
      id: 7, exposure_id: 5, q: 0.5,
      intensity: 1.2, prominence: 0.8, sharpness: 30.0,
      source: "auto",
      excluded: false,
    });
  });

  // -------------------------------------------------------------------------
  // Compare page mutators (Phase 3). Three event-shape rows.
  // -------------------------------------------------------------------------

  it("saveComparison.onSuccess invalidates when response lacks full Comparison shape (SSE-wins, comparison_created)", () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42, title: "stale" });
    let invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidatedKeys.push(arg.queryKey); return orig(arg);
    }) as typeof qc.invalidateQueries;
    // Synth from `synthesizeResponseFromSse`: title/description/members ride
    // as INPUT shape (no `id`s on members, no `is_stale`, no `content_hash`).
    // The mutator's onSuccess must NOT splice this into cache — it would
    // pollute with a malformed Comparison shape.
    const sseSynth = {
      event_id: 7, client_op_id: "op-cmp-create",
      analysis_inputs_hash: undefined,
      id: 42, title: "X", description: null,
      members: [{ exposure_id: 100, display_order: 0,
                  snapshot: { effective_peaks: [], confirmed_index: null,
                              analysis_inputs_hash: "h" } }],
    } as any;
    saveComparisonMutator.onSuccess(
      { id: 42, title: "X", members: [],
        username: "u", clientId: "c", clientOpId: "op-cmp-create" } as any,
      sseSynth, qc,
    );
    // The half-baked synth was NOT spliced. Cache is invalidated for refetch.
    expect(qc.getQueryData(queryKeys.comparison(42))).toEqual({ id: 42, title: "stale" });
    expect(invalidatedKeys).toContainEqual(queryKeys.comparison(42));
    expect(invalidatedKeys).toContainEqual(queryKeys.comparisonMembers(42));
    expect(invalidatedKeys).toContainEqual(["comparisons"]);
  });

  it("saveComparison.onSuccess splices when response has full Comparison shape (HTTP-wins, comparison_submitted)", () => {
    const fullResponse = {
      id: 42, title: "X edited", description: null,
      content_hash: "sha256:new",
      created_by: 1, created_at: "2026-05-06", updated_at: "2026-05-06",
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      members: [
        { id: 999, comparison_id: 42, exposure_id: 100, display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null,
          peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null,
                      analysis_inputs_hash: "sha256:zero" },
          is_stale: false, created_by: 1, created_at: "2026-05-06" },
      ],
    } as any;
    saveComparisonMutator.onSuccess(
      { id: 42, title: "X edited", members: [],
        expected_content_hash: "sha256:base",
        username: "u", clientId: "c", clientOpId: "op-cmp-submit" } as any,
      fullResponse, qc,
    );
    expect(qc.getQueryData(queryKeys.comparison(42))).toEqual(fullResponse);
    expect(qc.getQueryData(queryKeys.comparisonMembers(42))).toEqual(fullResponse.members);
  });

  it("deleteComparison.onSuccess removes entity caches and prunes listings (uniform across HTTP- and SSE-wins)", () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42 });
    qc.setQueryData(queryKeys.comparisons("all"), [
      { id: 42, title: "doomed", description: null, content_hash: "h",
        created_by: 1, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null },
    ]);
    // SSE-wins synth for comparison_deleted is just `{event_id, client_op_id,
    // analysis_inputs_hash, id}` — and the mutator only reads `p.id` from the
    // flat input, never from the response. Pin: same observable behavior.
    const sseSynth = {
      event_id: 7, client_op_id: "op-cmp-del",
      analysis_inputs_hash: undefined, id: 42,
    } as any;
    deleteComparisonMutator.onSuccess(
      { id: 42, username: "u", clientId: "c", clientOpId: "op-cmp-del" } as any,
      sseSynth, qc,
    );
    expect(qc.getQueryState(queryKeys.comparison(42))).toBeUndefined();
    const listing = qc.getQueryData<{ id: number }[]>(queryKeys.comparisons("all"))!;
    expect(listing).toEqual([]);
  });

  it("updateSample.onSuccess skips undefined fields (SSE-wins diff payload)", () => {
    const original = {
      id: 5, experiment_id: 1, label: "S5", name: "old name",
      notes: "old notes", tags: [{ id: 100, key: "k", value: "v", source: "manual" }],
    };
    qc.setQueryData(queryKeys.sample(5), original);
    qc.setQueryData(queryKeys.samples(1), [original]);
    // SSE payload for update_sample is the diff; if only `name` was patched,
    // `notes` and `label` are undefined in the synthesized response.
    const sseSynth = {
      event_id: 7,
      client_op_id: "op-update",
      analysis_inputs_hash: undefined,
      name: "new name",
      // notes: undefined, label: undefined  ← intentionally absent
    } as any;
    updateSampleMutator.onSuccess(
      { sampleId: 5, experimentId: 1, name: "new name",
        username: "u", clientId: "c", clientOpId: "op-update" } as any,
      sseSynth,
      qc,
    );
    const single = qc.getQueryData<any>(queryKeys.sample(5));
    // name updated; notes/label/tags preserved (NOT clobbered to undefined)
    expect(single).toEqual({
      id: 5, experiment_id: 1, label: "S5", name: "new name",
      notes: "old notes",
      tags: [{ id: 100, key: "k", value: "v", source: "manual" }],
    });
    const list = qc.getQueryData<any[]>(queryKeys.samples(1));
    expect(list![0]).toEqual({
      id: 5, experiment_id: 1, label: "S5", name: "new name",
      notes: "old notes",
      tags: [{ id: 100, key: "k", value: "v", source: "manual" }],
    });
  });
});
