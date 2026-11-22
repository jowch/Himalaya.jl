import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";

// Spec §4.4 + §6 — state.ts additions for permalink URL handling.

describe("Zustand state — permalink slots", () => {
  beforeEach(() => {
    // Reset relevant slots without nuking the whole store (preserves persisted
    // username so other tests don't fight unexpectedly).
    useAppState.setState({
      staleUrlContext: null,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      paletteOpen: false,
    });
  });

  it("setStaleUrlContext stores and clears", () => {
    const ctx = {
      kind: "not_found" as const,
      missing: "experiment" as const,
      missing_value: "lipid-typo",
      experiment_resolved: undefined,
      sample_resolved: undefined,
    };
    useAppState.getState().setStaleUrlContext(ctx);
    expect(useAppState.getState().staleUrlContext).toEqual(ctx);
    useAppState.getState().setStaleUrlContext(null);
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  // I5.1 (#182): setActivePage is deleted with the dual-nav model; the
  // setActive{Experiment,Sample,Exposure} setters still clear staleUrlContext.
  it.each([
    ["setActiveExperiment", (id: number) => useAppState.getState().setActiveExperiment(id)],
    ["setActiveSample",     (id: number) => useAppState.getState().setActiveSample(id)],
    ["setActiveExposure",   (id: number) => useAppState.getState().setActiveExposure(id)],
  ])("%s clears staleUrlContext", (_label, fn) => {
    useAppState.getState().setStaleUrlContext({
      kind: "unknown_path", raw: "/foo/bar",
    });
    fn(42);
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("recoverFromStaleUrl clears stale and opens the palette", () => {
    useAppState.getState().setStaleUrlContext({
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    });
    useAppState.getState().recoverFromStaleUrl();
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.paletteOpen).toBe(true);
  });

  it("staleUrlContext is NOT in the persisted slice", () => {
    useAppState.getState().setStaleUrlContext({ kind: "unknown_path", raw: "/x" });
    const persisted = JSON.parse(localStorage.getItem("himalaya-ui:state") ?? "{}");
    expect(persisted.state?.staleUrlContext).toBeUndefined();
  });

  it("setStaleUnknownPath stores raw", () => {
    useAppState.getState().setStaleUnknownPath("/foo/bar/baz");
    const s = useAppState.getState();
    expect(s.staleUrlContext).toEqual({ kind: "unknown_path", raw: "/foo/bar/baz" });
  });

  it("setStaleNotFound commits 404 context atomically", () => {
    const ctx = {
      kind: "not_found" as const,
      missing: "sample" as const,
      missing_value: "lipid-typo",
      experiment_resolved: { id: 3, name: "ExpA" },
      sample_resolved: undefined,
    };
    useAppState.getState().setStaleNotFound(ctx);
    const s = useAppState.getState();
    expect(s.staleUrlContext).toEqual(ctx);
  });
});
