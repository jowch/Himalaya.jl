import { describe, it, expect } from "vitest";
import { mergePersistedState, useAppState } from "../src/state";

// I5.1 (#182): the dual-nav `activePage` model + `coerceActivePage` are deleted.
// `mergePersistedState` survives as the explicit persist `merge` strategy
// (shallow: persisted keys override the live defaults, store actions are
// preserved). These assertions are salvaged from the old coerceActivePage.test.ts
// minus the (now-removed) `activePage` coercion cases.

describe("mergePersistedState (persist merge)", () => {
  it("lets persisted keys override the current defaults", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState({ theme: "light", username: "alice" }, current);
    expect(merged.theme).toBe("light");
    expect(merged.username).toBe("alice");
  });

  it("preserves store actions through the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState({ theme: "light" }, current);
    expect(typeof merged.setTheme).toBe("function");
    expect(typeof merged.setResolveSuccess).toBe("function");
  });

  it("handles undefined persisted state (first run) without dropping defaults", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(undefined, current);
    expect(merged.theme).toBe(current.theme);
    expect(typeof merged.setTheme).toBe("function");
  });

  it("tolerates a stale extra key from a pre-cutover blob (e.g. activePage)", () => {
    // A pre-I5.1 persisted blob still carries an `activePage` key. The shallow
    // merge leaves it as an inert extra property — it must not throw and must
    // not disturb the real persisted keys. #183 (I5.2) formally strips it via
    // a version-bump + migrate.
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "compare", theme: "light" } as Record<string, unknown>,
      current,
    );
    expect(merged.theme).toBe("light");
  });
});
