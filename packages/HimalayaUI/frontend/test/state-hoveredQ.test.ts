import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ hoveredQ: undefined });
});

describe("hoveredQ ephemeral field", () => {
  it("defaults to undefined and round-trips via setHoveredQ", () => {
    expect(useAppState.getState().hoveredQ).toBeUndefined();
    useAppState.getState().setHoveredQ(0.123);
    expect(useAppState.getState().hoveredQ).toBe(0.123);
    useAppState.getState().setHoveredQ(undefined);
    expect(useAppState.getState().hoveredQ).toBeUndefined();
  });

  it("is NOT written to the persisted localStorage payload", () => {
    // Guarantee a fresh persist envelope independent of write-on-every-set
    // timing: drive a PARTITIONED field first, so the payload definitely
    // exists with a real field, then set the ephemeral one. (R0a #221: the
    // former `theme` driver is gone; `tutorialSeen` is another partitioned field.)
    useAppState.getState().setTutorialSeen(true);
    useAppState.getState().setHoveredQ(0.42);
    const raw = localStorage.getItem(LS_KEY);
    expect(raw).not.toBeNull();
    const persisted = JSON.parse(raw as string);
    // persist wraps state under `.state`
    expect("tutorialSeen" in (persisted.state ?? {})).toBe(true); // partitioned field present
    expect("hoveredQ" in (persisted.state ?? {})).toBe(false);    // ephemeral field omitted
  });
});
