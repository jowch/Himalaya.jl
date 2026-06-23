import { describe, it, expect } from "vitest";
import { effectiveIngestStatus } from "../src/lib/ingestStatus";

describe("effectiveIngestStatus", () => {
  it("terminal persisted state wins over a stale scanning overlay (8c self-heal)", () => {
    // The missed-ingest_complete-frame case: inFlight stuck "scanning", DB says complete.
    expect(effectiveIngestStatus("scanning", "complete")).toBe("complete");
    expect(effectiveIngestStatus("scanning", "failed")).toBe("failed");
  });

  it("a rescan overlay (analyzing) leads even when persisted is terminal", () => {
    // A rescan keeps the persisted row "complete"; the overlay is the truth.
    expect(effectiveIngestStatus("analyzing", "complete")).toBe("analyzing");
    expect(effectiveIngestStatus("analyzing", "idle")).toBe("analyzing");
  });

  it("initial-scan overlay leads while persisted is non-terminal", () => {
    expect(effectiveIngestStatus("scanning", "scanning")).toBe("scanning");
    expect(effectiveIngestStatus("scanning", "idle")).toBe("scanning");
  });

  it("falls back to persisted then idle when no overlay", () => {
    expect(effectiveIngestStatus(undefined, "scanning")).toBe("scanning");
    expect(effectiveIngestStatus(undefined, undefined)).toBe("idle");
  });
});
