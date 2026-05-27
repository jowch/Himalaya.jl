import { describe, it, expect } from "vitest";
import {
  coerceActivePage,
  mergePersistedState,
  useAppState,
} from "../src/state";

// I3.6 (#177): every legacy AppShell surface is retired — Inspect (#163),
// Index (#181), Compare (#177) — so `PageId` is narrowed to the inert "none"
// sentinel. `coerceActivePage` always resolves to "none"; it is a pure type
// artifact (bare `/` redirects to `/samples`, no user lands on a page driven
// by `activePage`). I5.1 deletes the whole model.

describe("coerceActivePage", () => {
  it("passes the sole valid PageId 'none' through unchanged", () => {
    expect(coerceActivePage("none")).toBe("none");
  });

  it("coerces a stale 'compare' to 'none' after Compare retirement (#177)", () => {
    expect(coerceActivePage("compare")).toBe("none");
  });

  it("coerces a stale 'index' to 'none' after Index retirement (#181)", () => {
    expect(coerceActivePage("index")).toBe("none");
  });

  it("coerces a stale 'inspect' to 'none' after Inspect retirement (#163)", () => {
    expect(coerceActivePage("inspect")).toBe("none");
  });

  it("coerces an unknown surface name to 'none'", () => {
    expect(coerceActivePage("loupe")).toBe("none");
    expect(coerceActivePage("samples")).toBe("none");
    expect(coerceActivePage("")).toBe("none");
  });

  it("coerces non-string / missing values to 'none'", () => {
    expect(coerceActivePage(undefined)).toBe("none");
    expect(coerceActivePage(null)).toBe("none");
    expect(coerceActivePage(42)).toBe("none");
  });
});

describe("mergePersistedState (persist merge)", () => {
  it("coerces a stale persisted activePage during the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "ghost", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("none");
  });

  it("keeps the valid persisted activePage and other persisted keys", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "none", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("none");
    expect(merged.theme).toBe("light");
  });

  it("preserves store actions through the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState({ activePage: "inspect" }, current);
    expect(typeof merged.setActivePage).toBe("function");
  });

  it("handles undefined persisted state (first run)", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(undefined, current);
    expect(merged.activePage).toBe("none");
  });
});
