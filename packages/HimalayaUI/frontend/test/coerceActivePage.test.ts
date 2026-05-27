import { describe, it, expect } from "vitest";
import {
  coerceActivePage,
  mergePersistedState,
  useAppState,
} from "../src/state";

describe("coerceActivePage", () => {
  it("passes the sole valid PageId 'compare' through unchanged", () => {
    expect(coerceActivePage("compare")).toBe("compare");
  });

  it("coerces a stale 'index' to 'compare' after Index retirement (#181)", () => {
    // Index is retired; a persisted activePage:"index" must not strand the
    // user on an empty PageBody — coerceActivePage rewrites it to a live page.
    expect(coerceActivePage("index")).toBe("compare");
  });

  it("coerces a stale 'inspect' to 'compare' after Inspect retirement (#163)", () => {
    expect(coerceActivePage("inspect")).toBe("compare");
  });

  it("coerces an unknown surface name to 'compare'", () => {
    expect(coerceActivePage("loupe")).toBe("compare");
    expect(coerceActivePage("samples")).toBe("compare");
    expect(coerceActivePage("")).toBe("compare");
  });

  it("coerces non-string / missing values to 'compare'", () => {
    expect(coerceActivePage(undefined)).toBe("compare");
    expect(coerceActivePage(null)).toBe("compare");
    expect(coerceActivePage(42)).toBe("compare");
  });
});

describe("mergePersistedState (persist merge)", () => {
  it("coerces a stale persisted activePage during the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "ghost", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("compare");
  });

  it("keeps a valid persisted activePage and other persisted keys", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "compare", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("compare");
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
    expect(merged.activePage).toBe("compare");
  });
});
