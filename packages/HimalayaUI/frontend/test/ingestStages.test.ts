import { describe, it, expect } from "vitest";
import {
  isIngestStage,
  stageLabel,
  INGEST_STAGES,
} from "../src/lib/ingestStages";

describe("isIngestStage / stageLabel", () => {
  it("accepts only the known wire values", () => {
    expect(isIngestStage("discovery")).toBe(true);
    expect(isIngestStage("analyzing")).toBe(true);
    expect(isIngestStage("thumbnails")).toBe(false);
    expect(isIngestStage("rescan")).toBe(false); // that's `phase`, not `stage`
    expect(isIngestStage(undefined)).toBe(false);
  });

  it("labels every stage", () => {
    for (const s of INGEST_STAGES) {
      expect(stageLabel(s)).toBeTruthy();
    }
  });
});
