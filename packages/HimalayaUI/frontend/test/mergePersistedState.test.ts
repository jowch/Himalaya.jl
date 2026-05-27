import { describe, it, expect } from "vitest";
import { mergePersistedState, useAppState } from "../src/state";

// I5.1 (#182): the dual-nav `activePage` model + `coerceActivePage` are deleted.
// `mergePersistedState` survives as the explicit persist `merge` strategy
// (shallow: persisted keys override the live defaults, store actions are
// preserved). These assertions are salvaged from the old coerceActivePage.test.ts
// minus the (now-removed) `activePage` coercion cases.
//
// R0a (#221): the dark `theme` field + `setTheme` are retired. These tests now
// assert the merge strategy against surviving prefs/actions instead.

describe("mergePersistedState (persist merge)", () => {
  it("lets persisted keys override the current defaults", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { tutorialSeen: true, username: "alice" },
      current,
    );
    expect(merged.tutorialSeen).toBe(true);
    expect(merged.username).toBe("alice");
  });

  it("preserves store actions through the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState({ tutorialSeen: true }, current);
    expect(typeof merged.setTutorialSeen).toBe("function");
    expect(typeof merged.recoverFromStaleUrl).toBe("function");
  });

  it("handles undefined persisted state (first run) without dropping defaults", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(undefined, current);
    expect(merged.tutorialSeen).toBe(current.tutorialSeen);
    expect(typeof merged.setTutorialSeen).toBe("function");
  });

  it("tolerates a stale extra key from a pre-cutover blob (activePage, theme)", () => {
    // A pre-cutover persisted blob may still carry dead `activePage` / `theme`
    // keys. The shallow merge leaves them as inert extra properties — it must
    // not throw and must not disturb the real persisted keys. #221's migrate
    // formally strips them via a version-bump + migrate.
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "compare", theme: "light", username: "dave" } as Record<
        string,
        unknown
      >,
      current,
    );
    expect(merged.username).toBe("dave");
  });
});
