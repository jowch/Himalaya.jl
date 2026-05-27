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

  // I5.1 (#182): the `activePage` field + `setActivePage` are deleted with the
  // dual-nav model. The dedicated activePage test is removed (no field to test).

  it("tutorialSeen defaults to false and can be set true", () => {
    expect(useAppState.getState().tutorialSeen).toBe(false);
    useAppState.getState().setTutorialSeen(true);
    expect(useAppState.getState().tutorialSeen).toBe(true);
  });

  // R0a (#221): "The Print" is the single identity. The `theme` field,
  // `setTheme`, and `ThemeId` were removed; the app renders Print by default
  // with no toggle. Guard against a regression that re-adds the field.
  it("has no `theme` field or `setTheme` action (Print is the single identity)", () => {
    const s = useAppState.getState() as Record<string, unknown>;
    expect("theme" in s).toBe(false);
    expect("setTheme" in s).toBe(false);
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
    useAppState.getState().setTutorialSeen(true);
    const raw = localStorage.getItem(LS_KEY);
    const parsed = JSON.parse(raw!);
    expect(parsed.state.username).toBe("alice");
    expect(parsed.state.activeExperimentId).toBe(4);
    expect(parsed.state.activeSampleId).toBe(12);
    // I5.1 (#182): `activePage` is no longer in the persisted partition.
    expect(parsed.state.activePage).toBeUndefined();
    expect(parsed.state.tutorialSeen).toBe(true);
    // R0a (#221): `theme` is no longer persisted (field removed).
    expect(parsed.state.theme).toBeUndefined();
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

  // I5.3 (#184): the compareXDomains field + setCompareXDomain action were
  // removed (Compare-only; the series builder holds xDomain in local useState).
});
