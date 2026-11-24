import { describe, expect, it } from "vitest";
import type { SampleBackground, Trace } from "../../src/api";
import { backgroundsToFetch, planMinusCompares, SKIP_NOTES } from "../../src/workspace/minusCompares";

const q = [0.05, 0.1, 0.15, 0.2, 0.25];
const trace = (level: number, monitors?: Trace["monitors"]): Trace => ({
  q: [...q], I: q.map(() => level), sigma: q.map(() => 0), ...(monitors ? { monitors } : {}),
});
const normalised = { transmitted_counts: 1, exposure_time: 1, counting_time: 1, normalised_by: "i2" as const, mode: "autosampler" as const };

const BUFFER_A = 90;
const BUFFER_B = 91;
const entry = (id: number, over: Partial<SampleBackground> = {}): SampleBackground => ({
  sample_id: id, ref_sample_id: BUFFER_A, source: "default", default_ref_sample_id: BUFFER_A, is_buffer: false, run_buffer_ids: [BUFFER_A], ...over,
});

const plan = (over: Partial<Parameters<typeof planMinusCompares>[0]> = {}) =>
  planMinusCompares({
    chips: [11],
    backgrounds: [entry(11)],
    traces: new Map([[11, trace(5, normalised)], [BUFFER_A, trace(2, normalised)]]),
    statuses: new Map(),
    hidden: new Set(),
    setupOf: () => 1,
    peakQs: [],
    ...over,
  });

describe("planMinusCompares", () => {
  it("subtracts a compared sample's own background, with the scale of that pair", () => {
    const [p] = plan();
    expect(p).toMatchObject({ kind: "drawn", id: 11, background: BUFFER_A });
    if (p?.kind !== "drawn") throw new Error("not drawn");
    expect(p.d).toEqual(q.map(() => 3));
    expect(p.q).toEqual(q);
    expect(p.scale).toEqual({ scale: 1, source: "normalised" });
  });

  it("gives two comparisons from two runs each their own buffer", () => {
    const out = plan({
      chips: [11, 12],
      backgrounds: [entry(11), entry(12, { ref_sample_id: BUFFER_B, default_ref_sample_id: BUFFER_B })],
      traces: new Map([[11, trace(5, normalised)], [12, trace(9, normalised)], [BUFFER_A, trace(2, normalised)], [BUFFER_B, trace(4, normalised)]]),
    });
    expect(out.map((p) => (p.kind === "drawn" ? p.d[0] : null))).toEqual([3, 5]);
  });

  it("uses the buffer the user chose instead of the default", () => {
    const out = plan({
      backgrounds: [entry(11, { ref_sample_id: BUFFER_B, source: "user" })],
      traces: new Map([[11, trace(5, normalised)], [BUFFER_A, trace(2, normalised)], [BUFFER_B, trace(4, normalised)]]),
    });
    expect(out[0]).toMatchObject({ kind: "drawn", background: BUFFER_B });
  });

  it("flags a pair that had to be matched at high q, for that pair alone", () => {
    const out = plan({
      chips: [11, 12],
      backgrounds: [entry(11), entry(12)],
      traces: new Map([[11, trace(5, normalised)], [12, trace(6, { ...normalised, beam_off: true })], [BUFFER_A, trace(2, normalised)]]),
    });
    expect(out.map((p) => (p.kind === "drawn" ? p.scale.source : null))).toEqual(["normalised", "beam_off"]);
  });

  it("matches at high q when the compared sample and its background are in different setups", () => {
    const out = plan({ setupOf: (id) => (id === BUFFER_A ? 2 : 1) });
    expect(out[0]).toMatchObject({ kind: "drawn", scale: { source: "other_setup" } });
  });

  it.each([
    ["is a buffer", { is_buffer: true, ref_sample_id: null, default_ref_sample_id: null }, "buffer"],
    ["has its background set to none", { ref_sample_id: null, source: "user" as const }, "none"],
    ["has no default and no choice", { ref_sample_id: null, default_ref_sample_id: null, run_buffer_ids: [] }, "no-default"],
  ])("does not draw a sample that %s", (_name, over, reason) => {
    expect(plan({ backgrounds: [entry(11, over)] })).toEqual([{ kind: "skipped", id: 11, reason }]);
  });

  it("does not draw a sample without a background module, such as a capillary", () => {
    expect(plan({ backgrounds: [] })).toEqual([{ kind: "skipped", id: 11, reason: "no-module" }]);
  });

  it("does not draw when the background's setup is hidden", () => {
    expect(plan({ hidden: new Set([BUFFER_A]) })).toEqual([{ kind: "skipped", id: 11, reason: "hidden" }]);
  });

  it("reads a background trace still loading as loading, and an absent one as missing", () => {
    const traces = new Map([[11, trace(5, normalised)]]);
    expect(plan({ traces, statuses: new Map([[BUFFER_A, "loading" as const]]) })).toEqual([{ kind: "skipped", id: 11, reason: "loading" }]);
    expect(plan({ traces, statuses: new Map([[BUFFER_A, "missing" as const]]) })).toEqual([{ kind: "skipped", id: 11, reason: "missing" }]);
    expect(plan({ traces })).toEqual([{ kind: "skipped", id: 11, reason: "missing" }]);
  });

  it("reads the compared sample's own trace the same way", () => {
    const traces = new Map([[BUFFER_A, trace(2, normalised)]]);
    expect(plan({ traces, statuses: new Map([[11, "loading" as const]]) })).toEqual([{ kind: "skipped", id: 11, reason: "loading-own" }]);
    expect(plan({ traces })).toEqual([{ kind: "skipped", id: 11, reason: "no-trace" }]);
  });

  it("has a plain note for every reason", () => {
    expect(Object.values(SKIP_NOTES).every((n) => n.length > 0)).toBe(true);
    expect(SKIP_NOTES.buffer).toBe("Not subtracted: this is a buffer");
    expect(SKIP_NOTES.loading).toBe("Loading the background trace");
    expect(SKIP_NOTES["loading-own"]).toBe("Loading this sample's trace");
  });
});

describe("backgroundsToFetch", () => {
  it("lists each background once, even when several samples of a run share it", () => {
    const backgrounds = [entry(11), entry(12), entry(13, { ref_sample_id: BUFFER_B })];
    expect(backgroundsToFetch([11, 12, 13], backgrounds, new Set())).toEqual([BUFFER_A, BUFFER_B]);
  });

  it("leaves out buffers, unset backgrounds, hidden setups and samples without a module", () => {
    const backgrounds = [
      entry(11, { is_buffer: true, ref_sample_id: null }),
      entry(12, { ref_sample_id: null }),
      entry(13, { ref_sample_id: 40 }),
    ];
    expect(backgroundsToFetch([11, 12, 13, 14], backgrounds, new Set([40]))).toEqual([]);
  });

  it("fetches nothing for chips that are not asked about", () => {
    expect(backgroundsToFetch([], [entry(11)], new Set())).toEqual([]);
  });
});
