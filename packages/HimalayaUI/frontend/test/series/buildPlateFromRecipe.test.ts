import { describe, it, expect } from "vitest";
import { buildPlateFromRecipe } from "../../src/lib/series/buildPlateFromRecipe";
import type { SeriesMember, SeriesSample } from "../../src/api";

// ── fixtures ─────────────────────────────────────────────────────────────────
function sample(id: number, sampleId: number, position: number, over: Partial<SeriesSample> = {}): SeriesSample {
  return { id, series_id: 5, sample_id: sampleId, position, pinned: false, excluded: false, ...over };
}

function member(exposureId: number | null, displayOrder: number, over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1000 + displayOrder,
    series_id: 5,
    exposure_id: exposureId,
    display_order: displayOrder,
    band_height: 1,
    y_offset: 0,
    normalization: "none",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: null,
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

// Picker map: sample n → exposure 100 + n (mirrors the corpus picker's
// indexing_exposure_id projection the page feeds in).
const resolve = (sampleId: number): number | null =>
  sampleId === 99 ? null : 100 + sampleId;

describe("buildPlateFromRecipe", () => {
  it("REPRO BU-RECIPENOOP (reorder): the plate follows RECIPE position order, not the stale members' display_order", () => {
    // Recipe says [2, 1] (a reorder just saved by the PATCH); the cached
    // members still hold the OLD order [exp 101, exp 102].
    const samples = [sample(11, 2, 0), sample(10, 1, 1)];
    const old = [member(101, 0), member(102, 1)];
    const { members, unresolvedSampleIds } = buildPlateFromRecipe(samples, old, resolve);
    expect(members.map((m) => m.exposure_id)).toEqual([102, 101]);
    expect(members.map((m) => m.display_order)).toEqual([0, 1]);
    expect(unresolvedSampleIds).toEqual([]);
  });

  it("REPRO BU-RECIPENOOP (add): an added recipe sample becomes a member with server-default props (identity fields only)", () => {
    // 3 recipe samples, only 2 old members — the dev-DB series-1 shape.
    const samples = [sample(10, 1, 0), sample(11, 2, 1), sample(12, 3, 2)];
    const old = [member(101, 0), member(102, 1)];
    const { members } = buildPlateFromRecipe(samples, old, resolve);
    expect(members).toHaveLength(3);
    const added = members[2]!;
    expect(added.exposure_id).toBe(103);
    expect(added.display_order).toBe(2);
    // New member sends ONLY identity — the commit route's
    // _series_member_payload fills the create-path defaults + snapshot.
    expect(Object.keys(added).sort()).toEqual(["display_order", "exposure_id"]);
  });

  it("a removed recipe sample's member drops from the plate", () => {
    const samples = [sample(10, 1, 0)];
    const old = [member(101, 0), member(102, 1)];
    const { members } = buildPlateFromRecipe(samples, old, resolve);
    expect(members.map((m) => m.exposure_id)).toEqual([101]);
  });

  it("carry-over: an old member's display props survive a reorder (matched by exposure_id)", () => {
    const samples = [sample(11, 2, 0), sample(10, 1, 1)];
    const old = [
      member(101, 0, { label_override: "ratio 1:50", band_height: 2.5, color_override: "#aa0000", q_window_min: 0.04 }),
      member(102, 1, { normalization: "max", y_offset: 0.3 }),
    ];
    const { members } = buildPlateFromRecipe(samples, old, resolve);
    // exp 101's props followed it to its NEW slot (display_order 1).
    const m101 = members.find((m) => m.exposure_id === 101)!;
    expect(m101.display_order).toBe(1);
    expect(m101.label_override).toBe("ratio 1:50");
    expect(m101.band_height).toBe(2.5);
    expect(m101.color_override).toBe("#aa0000");
    expect(m101.q_window_min).toBe(0.04);
    const m102 = members.find((m) => m.exposure_id === 102)!;
    expect(m102.normalization).toBe("max");
    expect(m102.y_offset).toBe(0.3);
  });

  it("carry-over keeps a non-null snapshot and omits a null one (server refills)", () => {
    const snap = { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" };
    const old = [member(101, 0, { snapshot: snap }), member(102, 1, { snapshot: null })];
    const { members } = buildPlateFromRecipe([sample(10, 1, 0), sample(11, 2, 1)], old, resolve);
    expect(members[0]!.snapshot).toEqual(snap);
    expect("snapshot" in members[1]!).toBe(false);
  });

  it("members are id-less on the wire (the dispatcher mints fresh PKs)", () => {
    const { members } = buildPlateFromRecipe([sample(10, 1, 0)], [member(101, 0)], resolve);
    for (const m of members) expect("id" in m).toBe(false);
  });

  it("an unresolvable sample is left out AND reported; display_order stays gapless", () => {
    // sample 99 resolves to null (no usable exposure).
    const samples = [sample(10, 1, 0), sample(11, 99, 1), sample(12, 2, 2)];
    const { members, unresolvedSampleIds } = buildPlateFromRecipe(samples, [], resolve);
    expect(members.map((m) => m.exposure_id)).toEqual([101, 102]);
    expect(members.map((m) => m.display_order)).toEqual([0, 1]);
    expect(unresolvedSampleIds).toEqual([99]);
  });

  it("a sample MISSING from the picker map (resolver returns undefined) counts as unresolvable", () => {
    const { members, unresolvedSampleIds } = buildPlateFromRecipe(
      [sample(10, 1, 0)], [], () => undefined,
    );
    expect(members).toEqual([]);
    expect(unresolvedSampleIds).toEqual([1]);
  });

  it("an excluded recipe row is skipped WITHOUT being reported (deliberate exclusion, not a failure)", () => {
    const samples = [sample(10, 1, 0), sample(11, 2, 1, { excluded: true })];
    const { members, unresolvedSampleIds } = buildPlateFromRecipe(samples, [], resolve);
    expect(members.map((m) => m.exposure_id)).toEqual([101]);
    expect(unresolvedSampleIds).toEqual([]);
  });

  it("recipe order ties on position break by recipe-row id (backend ORDER BY position, id)", () => {
    const samples = [sample(20, 2, 0), sample(10, 1, 0)];
    const { members } = buildPlateFromRecipe(samples, [], resolve);
    expect(members.map((m) => m.exposure_id)).toEqual([101, 102]);
  });

  it("consume-once: two recipe samples resolving to the SAME exposure cannot both inherit one old member", () => {
    const samples = [sample(10, 1, 0), sample(11, 7, 1)];
    const old = [member(101, 0, { label_override: "kept once" })];
    const sameExposure = (): number => 101;
    const { members } = buildPlateFromRecipe(samples, old, sameExposure);
    expect(members).toHaveLength(2);
    expect(members[0]!.label_override).toBe("kept once");
    expect("label_override" in members[1]!).toBe(false);
  });
});
