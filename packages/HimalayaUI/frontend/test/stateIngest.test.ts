// test/stateIngest.test.ts
import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";

describe("ingestInFlight store slice (Phase E1)", () => {
  beforeEach(() => {
    useAppState.setState({ ingestInFlight: null });
  });

  it("starts null", () => {
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("setIngestProgress upserts one experiment's progress", () => {
    useAppState.getState().setIngestProgress(7, { processed: 100, total: 680, status: "scanning" });
    expect(useAppState.getState().ingestInFlight).toEqual({
      7: { processed: 100, total: 680, status: "scanning" },
    });
    useAppState.getState().setIngestProgress(7, { processed: 300, total: 680, status: "scanning" });
    expect(useAppState.getState().ingestInFlight![7]!.processed).toBe(300);
  });

  it("clearIngestProgress removes one experiment, nulling the map when empty", () => {
    useAppState.getState().setIngestProgress(7, { processed: 1, total: 2, status: "scanning" });
    useAppState.getState().clearIngestProgress(7);
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("is NOT included in the persisted partition", () => {
    // partialize (state.ts:500-508) must not surface ingestInFlight (ephemeral).
    // "himalaya-ui:state" is the real persist key — it is `LS_KEY` (state.ts:24),
    // the `name` passed to zustand persist (state.ts:489). Assert the serialized
    // blob never carries the key.
    useAppState.getState().setIngestProgress(7, { processed: 1, total: 2, status: "scanning" });
    const raw = localStorage.getItem("himalaya-ui:state");
    expect(raw === null || !raw.includes("ingestInFlight")).toBe(true);
  });
});
