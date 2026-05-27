import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";

// Spec §4.4 + §6 — state.ts additions for permalink URL handling.

describe("Zustand state — permalink slots", () => {
  beforeEach(() => {
    // Reset relevant slots without nuking the whole store (preserves persisted
    // username so other tests don't fight unexpectedly).
    useAppState.setState({
      staleUrlContext: null,
      resolving: false,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",
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

  it("setResolving toggles", () => {
    expect(useAppState.getState().resolving).toBe(false);
    useAppState.getState().setResolving(true);
    expect(useAppState.getState().resolving).toBe(true);
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

  it("setActive* does NOT clear resolving", () => {
    useAppState.getState().setResolving(true);
    useAppState.getState().setActiveExperiment(5);
    expect(useAppState.getState().resolving).toBe(true);
  });

  it("recoverFromStaleUrl row 1 (experiment): clears stale, opens modal at experiment step", () => {
    useAppState.getState().setStaleUrlContext({
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    });
    useAppState.getState().recoverFromStaleUrl({
      step: "experiment",
      experimentId: undefined,
      sampleId: undefined,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeSampleId).toBeUndefined();
    expect(s.activeExposureId).toBeUndefined();
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("experiment");
  });

  it("recoverFromStaleUrl row 2 (sample): sets experimentId, opens modal at sample step", () => {
    useAppState.getState().recoverFromStaleUrl({
      step: "sample",
      experimentId: 17,
      sampleId: undefined,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBeUndefined();
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("sample");
  });

  it("recoverFromStaleUrl row 3 (exposure): preserves sample, openModal=false", () => {
    useAppState.getState().recoverFromStaleUrl({
      step: "sample",
      experimentId: 17,
      sampleId: 42,
      openModal: false,
    });
    const s = useAppState.getState();
    expect(s.staleUrlContext).toBeNull();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBe(42);
    expect(s.activeExposureId).toBeUndefined();
    expect(s.navModalOpen).toBe(false);
    expect(s.navModalStep).toBe("sample");
  });

  it("staleUrlContext and resolving are NOT in the persisted slice", () => {
    // Touch them, then read back the Zustand persist `partialize` output.
    useAppState.getState().setResolving(true);
    useAppState.getState().setStaleUrlContext({ kind: "unknown_path", raw: "/x" });
    const persisted = JSON.parse(localStorage.getItem("himalaya-ui:state") ?? "{}");
    expect(persisted.state?.resolving).toBeUndefined();
    expect(persisted.state?.staleUrlContext).toBeUndefined();
  });

  it("setStaleUnknownPath stores raw + clears resolving", () => {
    useAppState.getState().setResolving(true);
    useAppState.getState().setStaleUnknownPath("/foo/bar/baz");
    const s = useAppState.getState();
    expect(s.staleUrlContext).toEqual({ kind: "unknown_path", raw: "/foo/bar/baz" });
    expect(s.resolving).toBe(false);
  });

  it("setStaleNotFound commits 404 context atomically + clears resolving", () => {
    useAppState.getState().setResolving(true);
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
    expect(s.resolving).toBe(false);
  });
});
