import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

describe("useAppState", () => {
  beforeEach(() => {
    localStorage.clear();
    useAppState.setState({ username: undefined, tutorialSeen: false, paletteOpen: false });
  });

  it("setUsername updates state", () => {
    useAppState.getState().setUsername("alice");
    expect(useAppState.getState().username).toBe("alice");
  });

  it("tutorialSeen defaults to false and can be set true", () => {
    expect(useAppState.getState().tutorialSeen).toBe(false);
    useAppState.getState().setTutorialSeen(true);
    expect(useAppState.getState().tutorialSeen).toBe(true);
  });

  it("opens and closes the palette", () => {
    useAppState.getState().openPalette();
    expect(useAppState.getState().paletteOpen).toBe(true);
    useAppState.getState().closePalette();
    expect(useAppState.getState().paletteOpen).toBe(false);
  });

  it("persists only the user and tutorial fields", () => {
    useAppState.getState().setUsername("alice");
    useAppState.getState().setTutorialSeen(true);
    useAppState.getState().openPalette();
    const parsed = JSON.parse(localStorage.getItem(LS_KEY)!);
    expect(parsed.state).toEqual({ username: "alice", tutorialSeen: true });
  });
});
