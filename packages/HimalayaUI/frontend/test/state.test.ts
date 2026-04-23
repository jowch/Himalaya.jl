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
});
