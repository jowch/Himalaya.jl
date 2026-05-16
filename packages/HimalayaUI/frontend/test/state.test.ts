import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

describe("useAppState", () => {
  beforeEach(() => {
    localStorage.clear();
    useAppState.setState({
      username: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
    });
  });

  it("starts with undefined fields", () => {
    const s = useAppState.getState();
    expect(s.username).toBeUndefined();
    expect(s.activeSampleId).toBeUndefined();
    expect(s.activeExposureId).toBeUndefined();
  });

  it("setUsername updates state", () => {
    useAppState.getState().setUsername("alice");
    expect(useAppState.getState().username).toBe("alice");
  });

  it("setActiveSample clears activeExposureId", () => {
    useAppState.setState({ activeExposureId: 7 });
    useAppState.getState().setActiveSample(3);
    expect(useAppState.getState().activeSampleId).toBe(3);
    expect(useAppState.getState().activeExposureId).toBeUndefined();
  });

  it("persists to localStorage under the stable key", () => {
    useAppState.getState().setUsername("bob");
    const raw = localStorage.getItem(LS_KEY);
    expect(raw).not.toBeNull();
    const parsed = JSON.parse(raw!);
    expect(parsed.state.username).toBe("bob");
  });

  it("hoveredIndexId starts undefined and can be set/cleared", () => {
    useAppState.setState({ hoveredIndexId: undefined });
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
    useAppState.getState().setHoveredIndex(7);
    expect(useAppState.getState().hoveredIndexId).toBe(7);
    useAppState.getState().setHoveredIndex(undefined);
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  });

  it("hoveredIndexId is NOT in the persisted partition", () => {
    useAppState.setState({ hoveredIndexId: 42 });
    const raw = localStorage.getItem(LS_KEY);
    expect(raw ?? "").not.toContain("hoveredIndexId");
  });

  // ── new fields added by the three-card redesign ────────────────────────

  it("activeExperimentId starts undefined and can be set", () => {
    expect(useAppState.getState().activeExperimentId).toBeUndefined();
    useAppState.getState().setActiveExperiment(5);
    expect(useAppState.getState().activeExperimentId).toBe(5);
  });

  it("setActiveExperiment clears activeSampleId and activeExposureId", () => {
    useAppState.setState({ activeSampleId: 9, activeExposureId: 2 });
    useAppState.getState().setActiveExperiment(3);
    const s = useAppState.getState();
    expect(s.activeExperimentId).toBe(3);
    expect(s.activeSampleId).toBeUndefined();
    expect(s.activeExposureId).toBeUndefined();
  });

  it("activePage defaults to 'index' and can switch to 'compare'", () => {
    expect(useAppState.getState().activePage).toBe("index");
    useAppState.getState().setActivePage("compare");
    expect(useAppState.getState().activePage).toBe("compare");
  });

  it("tutorialSeen defaults to false and can be set true", () => {
    expect(useAppState.getState().tutorialSeen).toBe(false);
    useAppState.getState().setTutorialSeen(true);
    expect(useAppState.getState().tutorialSeen).toBe(true);
  });

  it("theme defaults to 'dark' and can toggle", () => {
    expect(useAppState.getState().theme).toBe("dark");
    useAppState.getState().setTheme("light");
    expect(useAppState.getState().theme).toBe("light");
  });

  it("navModal state is ephemeral — open/close + step transitions", () => {
    expect(useAppState.getState().navModalOpen).toBe(false);
    useAppState.getState().openNavModal();
    expect(useAppState.getState().navModalOpen).toBe(true);
    useAppState.getState().setNavModalStep("sample");
    expect(useAppState.getState().navModalStep).toBe("sample");
    useAppState.getState().closeNavModal();
    expect(useAppState.getState().navModalOpen).toBe(false);
  });

  it("persists the allow-listed fields", () => {
    useAppState.getState().setUsername("alice");
    useAppState.getState().setActiveExperiment(4);
    useAppState.setState({ activeSampleId: 12 });
    useAppState.getState().setActivePage("compare");
    useAppState.getState().setTutorialSeen(true);
    useAppState.getState().setTheme("light");
    const raw = localStorage.getItem(LS_KEY);
    const parsed = JSON.parse(raw!);
    expect(parsed.state.username).toBe("alice");
    expect(parsed.state.activeExperimentId).toBe(4);
    expect(parsed.state.activeSampleId).toBe(12);
    expect(parsed.state.activePage).toBe("compare");
    expect(parsed.state.tutorialSeen).toBe(true);
    expect(parsed.state.theme).toBe("light");
  });

  it("does NOT persist ephemeral UI fields (navModal*, hoveredIndexId)", () => {
    useAppState.setState({
      hoveredIndexId: 3,
      navModalOpen: true,
      navModalStep: "sample",
    });
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("hoveredIndexId");
    expect(raw).not.toContain("navModalOpen");
    expect(raw).not.toContain("navModalStep");
  });

  it("hoveredPeakId starts undefined and can be set/cleared", () => {
    useAppState.setState({ hoveredPeakId: undefined });
    expect(useAppState.getState().hoveredPeakId).toBeUndefined();
    useAppState.getState().setHoveredPeak(5);
    expect(useAppState.getState().hoveredPeakId).toBe(5);
    useAppState.getState().setHoveredPeak(undefined);
    expect(useAppState.getState().hoveredPeakId).toBeUndefined();
  });

  it("hoveredPeakId is NOT in the persisted partition", () => {
    useAppState.setState({ hoveredPeakId: 99 });
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("hoveredPeakId");
  });

  // ── compare-page review-mode UI state (Phase 9) ────────────────────────

  // C-4: groupingMode was removed from Zustand; it now lives on
  // ActiveDraft.viewGroupingMode and is resolved via effectiveGroupingMode.
  it("setDraftViewGroupingMode creates an empty draft and sets viewGroupingMode", () => {
    // Start with no draft.
    useAppState.setState({ activeDraft: null });
    useAppState.getState().setDraftViewGroupingMode("byPhase");
    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    expect(draft?.viewGroupingMode).toBe("byPhase");
  });

  it("setDraftViewGroupingMode updates existing draft's viewGroupingMode", () => {
    useAppState.getState().setDraftViewGroupingMode("bySample");
    useAppState.getState().setDraftViewGroupingMode("distinct");
    expect(useAppState.getState().activeDraft?.viewGroupingMode).toBe("distinct");
  });

  it("viewGroupingMode is NOT in the persisted partition (carried on draft, not LS_KEY state)", () => {
    // The draft is sessionStorage-persisted, not localStorage; the main store
    // key should not contain the raw string "viewGroupingMode".
    useAppState.getState().setDraftViewGroupingMode("byPhase");
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("viewGroupingMode");
  });

  it("showPeakTicks defaults to true and can be set", () => {
    expect(useAppState.getState().showPeakTicks).toBe(true);
    useAppState.getState().setShowPeakTicks(false);
    expect(useAppState.getState().showPeakTicks).toBe(false);
    useAppState.getState().setShowPeakTicks(true);
    expect(useAppState.getState().showPeakTicks).toBe(true);
  });

  it("showPeakLabels defaults to true and can be set", () => {
    expect(useAppState.getState().showPeakLabels).toBe(true);
    useAppState.getState().setShowPeakLabels(false);
    expect(useAppState.getState().showPeakLabels).toBe(false);
  });

  it("annotation toggles are NOT persisted", () => {
    useAppState.getState().setShowPeakTicks(false);
    useAppState.getState().setShowPeakLabels(false);
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("showPeakTicks");
    expect(raw).not.toContain("showPeakLabels");
  });

  it("highlightedCompareMemberId starts undefined and can set/clear via single setter", () => {
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
    useAppState.getState().setHighlightedCompareMemberId(42);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    useAppState.getState().setHighlightedCompareMemberId(undefined);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });

  it("highlightedCompareMemberId is NOT persisted", () => {
    useAppState.getState().setHighlightedCompareMemberId(7);
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("highlightedCompareMemberId");
  });

  // ── compareXDomains is keyed per-comparison (Phase 6 reviewer fix) ──────
  // A single shared zoom domain bled across comparisons: zooming A to a
  // narrow q-range and then navigating to B (with a different q-range) made
  // B look broken. Per-comparison keying preserves zoom across review/edit
  // toggles for the SAME comparison while isolating different ones.

  it("compareXDomains starts as an empty record and is per-comparison-id", () => {
    useAppState.setState({ compareXDomains: {} });
    expect(useAppState.getState().compareXDomains).toEqual({});
  });

  it("setCompareXDomain(id, d) scopes the domain to that id only", () => {
    useAppState.setState({ compareXDomains: {} });
    useAppState.getState().setCompareXDomain(1, [0.1, 0.3]);
    useAppState.getState().setCompareXDomain(2, [0.5, 1.0]);
    const s = useAppState.getState();
    expect(s.compareXDomains[1]).toEqual([0.1, 0.3]);
    expect(s.compareXDomains[2]).toEqual([0.5, 1.0]);
  });

  it("setting the domain for one comparison does not affect another", () => {
    useAppState.setState({ compareXDomains: {} });
    useAppState.getState().setCompareXDomain(1, [0.1, 0.3]);
    useAppState.getState().setCompareXDomain(2, [0.5, 1.0]);
    // Mutate id=1, ensure id=2 untouched.
    useAppState.getState().setCompareXDomain(1, [0.2, 0.4]);
    const s = useAppState.getState();
    expect(s.compareXDomains[1]).toEqual([0.2, 0.4]);
    expect(s.compareXDomains[2]).toEqual([0.5, 1.0]);
  });

  it("setCompareXDomain(id, null) clears the entry", () => {
    useAppState.setState({ compareXDomains: {} });
    useAppState.getState().setCompareXDomain(1, [0.1, 0.3]);
    useAppState.getState().setCompareXDomain(1, null);
    expect(useAppState.getState().compareXDomains[1] ?? null).toBeNull();
  });

  it("compareXDomains is NOT in the persisted partition (ephemeral UI state)", () => {
    useAppState.getState().setCompareXDomain(1, [0.1, 0.3]);
    const raw = localStorage.getItem(LS_KEY) ?? "";
    expect(raw).not.toContain("compareXDomains");
  });
});
