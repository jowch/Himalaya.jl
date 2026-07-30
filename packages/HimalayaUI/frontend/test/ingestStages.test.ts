import { describe, it, expect } from "vitest";
import {
  deriveSegments,
  isIngestStage,
  stageLabel,
  INGEST_STAGES,
} from "../src/lib/ingestStages";

describe("deriveSegments", () => {
  it("returns null for an absent or unknown stage so callers fall back", () => {
    // An `ingest_started` frame carries no stage, and a backend predating stage
    // reporting never sends one. Rendering an all-empty segment strip would be
    // worse than the plain single-track bar.
    expect(deriveSegments(undefined, 0, 0)).toBeNull();
    expect(deriveSegments("", 0, 0)).toBeNull();
    expect(deriveSegments("nonsense", 5, 10)).toBeNull();
  });

  it("fills the active segment against its own denominator", () => {
    const segs = deriveSegments("discovery", 780, 1100);
    expect(segs).not.toBeNull();
    const discovery = segs!.find((s) => s.key === "discovery")!;
    expect(discovery.active).toBe(true);
    expect(discovery.fraction).toBeCloseTo(780 / 1100);
    expect(discovery.processed).toBe(780);
    expect(discovery.total).toBe(1100);
  });

  it("treats earlier stages as complete and later ones as untouched", () => {
    const segs = deriveSegments("analyzing", 92, 604)!;
    const byKey = Object.fromEntries(segs.map((s) => [s.key, s]));
    // discovery already ran to completion
    expect(byKey.discovery.fraction).toBe(1);
    expect(byKey.discovery.active).toBe(false);
    // thumbnails has not started
    expect(byKey.thumbnails.fraction).toBe(0);
    expect(byKey.thumbnails.active).toBe(false);
    // exactly one active segment at a time
    expect(segs.filter((s) => s.active)).toHaveLength(1);
  });

  it("renders a zero-total stage as COMPLETE, not empty", () => {
    // A clean rescan finds no new exposures, so scan_and_group! closes the
    // analyzing segment as 0-of-0. At 0% the bar would look stalled for the rest
    // of the scan; 0-of-0 means "nothing to do here", i.e. done.
    const segs = deriveSegments("analyzing", 0, 0)!;
    const analyzing = segs.find((s) => s.key === "analyzing")!;
    expect(analyzing.fraction).toBe(1);
    expect(analyzing.active).toBe(true);
  });

  it("clamps a fraction that would overflow or go negative", () => {
    expect(deriveSegments("discovery", 999, 10)![0].fraction).toBe(1);
    expect(deriveSegments("discovery", -5, 10)![0].fraction).toBe(0);
  });

  it("always returns one segment per stage, in execution order", () => {
    const segs = deriveSegments("thumbnails", 1, 2)!;
    expect(segs.map((s) => s.key)).toEqual([...INGEST_STAGES]);
    // Every segment carries a human label for the caption + a11y text.
    expect(segs.every((s) => s.label.length > 0)).toBe(true);
  });

  it("advances monotonically across a whole scan (never resets)", () => {
    // The point of segmenting: total fill only ever increases, even though each
    // stage restarts its own count.
    const totalFill = (stage: string, p: number, t: number) =>
      deriveSegments(stage, p, t)!.reduce((a, s) => a + s.fraction, 0);

    const timeline = [
      totalFill("discovery", 0, 1100),
      totalFill("discovery", 550, 1100),
      totalFill("discovery", 1100, 1100),
      totalFill("analyzing", 0, 604),
      totalFill("analyzing", 604, 604),
      totalFill("thumbnails", 0, 604),
      totalFill("thumbnails", 604, 604),
    ];
    for (let i = 1; i < timeline.length; i++) {
      expect(timeline[i]).toBeGreaterThanOrEqual(timeline[i - 1]);
    }
    expect(timeline.at(-1)).toBe(3); // all three segments full
  });
});

describe("isIngestStage / stageLabel", () => {
  it("accepts only the known wire values", () => {
    expect(isIngestStage("discovery")).toBe(true);
    expect(isIngestStage("analyzing")).toBe(true);
    expect(isIngestStage("thumbnails")).toBe(true);
    expect(isIngestStage("rescan")).toBe(false); // that's `phase`, not `stage`
    expect(isIngestStage(undefined)).toBe(false);
  });

  it("labels every stage", () => {
    for (const s of INGEST_STAGES) {
      expect(stageLabel(s)).toBeTruthy();
    }
  });
});
