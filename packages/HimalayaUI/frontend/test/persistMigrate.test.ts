import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

describe("persist migrate (v5 → v6)", () => {
  beforeEach(() => {
    localStorage.clear();
    useAppState.setState({ username: undefined, tutorialSeen: false });
  });

  it("keeps the user's prefs and strips the retired keys", async () => {
    localStorage.setItem(
      LS_KEY,
      JSON.stringify({
        state: {
          username: "alice",
          theme: "light",
          tutorialSeen: true,
          activeExperimentId: 7,
          activeSampleId: 12,
          activeExposureId: 3,
          activePage: "compare",
        },
        version: 5,
      }),
    );
    await useAppState.persist.rehydrate();
    const s: Record<string, unknown> = useAppState.getState() as unknown as Record<string, unknown>;
    expect(s["username"]).toBe("alice");
    expect(s["tutorialSeen"]).toBe(true);
    for (const key of ["activePage", "theme", "activeExperimentId", "activeSampleId", "activeExposureId"]) {
      expect(key in s).toBe(false);
    }
    const reparsed = JSON.parse(localStorage.getItem(LS_KEY)!);
    expect(reparsed.version).toBe(6);
    expect(reparsed.state.username).toBe("alice");
    expect(reparsed.state.activeExperimentId).toBeUndefined();
  });

  it("passes a non-object blob through untouched", () => {
    const migrate = useAppState.persist.getOptions().migrate!;
    expect(migrate(undefined, 5)).toBeUndefined();
    expect(migrate("corrupt", 5)).toBe("corrupt");
    expect(migrate(42, 5)).toBe(42);
  });
});
