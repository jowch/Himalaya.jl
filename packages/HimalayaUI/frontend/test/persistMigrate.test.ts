import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

// I5.2 (#183): the Zustand `persist` `version` is bumped 3 → 4 WITH a real
// `migrate` that strips the dead `activePage` key (left inert by I5.1 #182's
// field deletion) while preserving every surviving pref. A version bump WITHOUT
// a migrate would make Zustand discard the whole persisted blob — the wipe risk
// these tests guard against.
//
// R0a (#221): "The Print" is the single identity. The dark `theme` field is
// removed, so the migrate now ALSO strips a lingering `theme` key from any
// pre-cutover (v4 and earlier) blob. `version` bumps 4 → 5. The same wipe-guard
// invariant holds: a malformed blob passes through untouched.
//
// One-responsibility-per-file: persist-version/migrate behavior lives here,
// mirroring the mergePersistedState.test.ts split. `state.test.ts` stays focused
// on the persist partition + setters.

describe("persist migrate (v4 → v5, #221)", () => {
  beforeEach(() => {
    localStorage.clear();
    // The store is a module singleton shared across tests; reset the prefs we
    // assert on so a prior test's rehydrate can't leak in.
    useAppState.setState({
      username: undefined,
      tutorialSeen: false,
      activeExperimentId: undefined,
    });
  });

  it("preserves prefs across a pre-cutover blob and strips activePage + theme", async () => {
    // Seed a pre-cutover envelope that still carries the dead `activePage` and
    // the retired `theme` key alongside real prefs.
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
        version: 4,
      }),
    );

    // Re-read storage through the configured migrate/merge.
    await useAppState.persist.rehydrate();

    const s = useAppState.getState();
    // Prefs preserved — NOT wiped.
    expect(s.username).toBe("alice");
    expect(s.tutorialSeen).toBe(true);
    expect(s.activeExperimentId).toBe(7);
    // Dead keys stripped from live state.
    expect("activePage" in s).toBe(false);
    expect("theme" in s).toBe(false);

    // …and stripped from the re-persisted envelope (version now 5).
    const reparsed = JSON.parse(localStorage.getItem(LS_KEY)!);
    expect(reparsed.version).toBe(5);
    expect(reparsed.state.activePage).toBeUndefined();
    expect(reparsed.state.theme).toBeUndefined();
    expect(reparsed.state.username).toBe("alice");
  });

  it("an old persisted `theme` value does not break hydration", async () => {
    // A blob carrying a stale `theme` must hydrate cleanly (no throw) and
    // surviving prefs must round-trip.
    localStorage.setItem(
      LS_KEY,
      JSON.stringify({
        state: { username: "carol", theme: "dark", activeSampleId: 12 },
        version: 4,
      }),
    );
    await expect(useAppState.persist.rehydrate()).resolves.not.toThrow();
    const s = useAppState.getState();
    expect(s.username).toBe("carol");
    expect(s.activeSampleId).toBe(12);
    expect("theme" in s).toBe(false);
  });

  it("migrate is a pure data transform that drops activePage + theme, keeps prefs", () => {
    const migrate = useAppState.persist.getOptions().migrate;
    expect(typeof migrate).toBe("function");
    const out = migrate!(
      { theme: "dark", username: "bob", activePage: "index" },
      4,
    ) as Record<string, unknown>;
    expect(out.activePage).toBeUndefined();
    expect(out.theme).toBeUndefined();
    expect(out.username).toBe("bob");
    expect(out).not.toBeUndefined();
  });

  it("malformed-blob wipe-guard: a non-object blob passes through untouched", () => {
    // The else-branch must return the input verbatim — never {}/undefined,
    // which would be the only way the migrate could itself partial-wipe prefs.
    const migrate = useAppState.persist.getOptions().migrate!;
    expect(migrate(undefined, 4)).toBeUndefined();
    expect(migrate("corrupt", 4)).toBe("corrupt");
    expect(migrate(42, 4)).toBe(42);
  });

  it("no stored blob → defaults intact, no wipe", async () => {
    localStorage.clear();
    await useAppState.persist.rehydrate();
    const s = useAppState.getState();
    expect("activePage" in s).toBe(false);
    expect("theme" in s).toBe(false);
  });
});
