import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

// I5.2 (#183): the Zustand `persist` `version` is bumped 3 → 4 WITH a real
// `migrate` that strips the dead `activePage` key (left inert by I5.1 #182's
// field deletion) while preserving every surviving pref. A version bump WITHOUT
// a migrate would make Zustand discard the whole persisted blob — the wipe risk
// these tests guard against.
//
// One-responsibility-per-file: persist-version/migrate behavior lives here,
// mirroring the mergePersistedState.test.ts split. `state.test.ts` stays focused
// on the persist partition + setters.

describe("persist migrate (v3 → v4, #183)", () => {
  beforeEach(() => {
    localStorage.clear();
    // The store is a module singleton shared across tests; reset the prefs we
    // assert on so a prior test's rehydrate can't leak in.
    useAppState.setState({
      username: undefined,
      theme: "dark",
      tutorialSeen: false,
      activeExperimentId: undefined,
    });
  });

  it("preserves prefs across a pre-cutover (v3) blob and strips activePage", async () => {
    // Seed a pre-cutover v3 envelope that still carries the dead `activePage`
    // key alongside real prefs.
    localStorage.setItem(
      LS_KEY,
      JSON.stringify({
        state: {
          username: "alice",
          theme: "light",
          tutorialSeen: true,
          activeExperimentId: 7,
          activePage: "compare",
        },
        version: 3,
      }),
    );

    // Re-read storage through the configured migrate/merge.
    await useAppState.persist.rehydrate();

    const s = useAppState.getState();
    // Prefs preserved — NOT wiped.
    expect(s.username).toBe("alice");
    expect(s.theme).toBe("light");
    expect(s.tutorialSeen).toBe(true);
    expect(s.activeExperimentId).toBe(7);
    // Dead key stripped from live state.
    expect("activePage" in s).toBe(false);

    // …and stripped from the re-persisted envelope (version now 4).
    const reparsed = JSON.parse(localStorage.getItem(LS_KEY)!);
    expect(reparsed.version).toBe(4);
    expect(reparsed.state.activePage).toBeUndefined();
    expect(reparsed.state.username).toBe("alice");
  });

  it("migrate is a pure data transform that drops activePage and keeps prefs", () => {
    const migrate = useAppState.persist.getOptions().migrate;
    expect(typeof migrate).toBe("function");
    const out = migrate!(
      { theme: "dark", username: "bob", activePage: "index" },
      3,
    ) as Record<string, unknown>;
    expect(out.activePage).toBeUndefined();
    expect(out.theme).toBe("dark");
    expect(out.username).toBe("bob");
    expect(out).not.toBeUndefined();
  });

  it("malformed-blob wipe-guard: a non-object blob passes through untouched", () => {
    // The else-branch must return the input verbatim — never {}/undefined,
    // which would be the only way the migrate could itself partial-wipe prefs.
    const migrate = useAppState.persist.getOptions().migrate!;
    expect(migrate(undefined, 3)).toBeUndefined();
    expect(migrate("corrupt", 3)).toBe("corrupt");
    expect(migrate(42, 3)).toBe(42);
  });

  it("no stored blob → defaults intact, no wipe", async () => {
    localStorage.clear();
    await useAppState.persist.rehydrate();
    const s = useAppState.getState();
    expect(s.theme).toBe("dark"); // default preserved
    expect("activePage" in s).toBe(false);
  });
});
