import { describe, it, expect } from "vitest";
import {
  coerceActivePage,
  mergePersistedState,
  useAppState,
} from "../src/state";

describe("coerceActivePage", () => {
  it("passes valid PageIds through unchanged", () => {
    expect(coerceActivePage("index")).toBe("index");
    expect(coerceActivePage("compare")).toBe("compare");
  });

  it("coerces a stale 'inspect' to 'index' after Inspect retirement (#163)", () => {
    // Inspect is retired; a persisted activePage:"inspect" must not strand the
    // user on an empty PageBody — coerceActivePage rewrites it to a live page.
    expect(coerceActivePage("inspect")).toBe("index");
  });

  it("coerces an unknown surface name to 'index'", () => {
    expect(coerceActivePage("loupe")).toBe("index");
    expect(coerceActivePage("samples")).toBe("index");
    expect(coerceActivePage("")).toBe("index");
  });

  it("coerces non-string / missing values to 'index'", () => {
    expect(coerceActivePage(undefined)).toBe("index");
    expect(coerceActivePage(null)).toBe("index");
    expect(coerceActivePage(42)).toBe("index");
  });
});

describe("mergePersistedState (persist merge)", () => {
  it("coerces a stale persisted activePage during the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "ghost", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("index");
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
    expect(merged.activePage).toBe("index");
  });
});
