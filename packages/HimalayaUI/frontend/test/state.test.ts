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
});
