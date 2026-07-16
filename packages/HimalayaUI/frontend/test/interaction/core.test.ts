import { describe, it, expect } from "vitest";
import { CORE, core, page, assertNoCoreCollision } from "../../src/print/interaction/core";

describe("CORE vocabulary", () => {
  it("openFocus is bound to Enter app-wide", () => {
    expect(CORE.openFocus.keys).toEqual(["Enter"]);
  });
  it("undo is Mod+z, redo Mod+Shift+z", () => {
    expect(CORE.undo.keys).toEqual(["Mod+z"]);
    expect(CORE.redo.keys).toEqual(["Mod+Shift+z"]);
  });
});

describe("core()", () => {
  it("fills label/keys/group from CORE, page supplies the handler", () => {
    const run = (): void => {};
    const a = core("openFocus", { run, dock: "primary" });
    expect(a).toMatchObject({ id: "openFocus", label: "Focus", keys: ["Enter"], group: "Navigate", dock: "primary" });
    expect(a.run).toBe(run);
  });
  it("lets the page override the label only", () => {
    expect(core("openFocus", { run: () => {}, label: "Apply" }).label).toBe("Apply");
  });
});

describe("page()", () => {
  it("builds a page-local verb verbatim", () => {
    const a = page("cull", { label: "Cull", keys: ["x"], group: "Act", run: () => {} });
    expect(a).toMatchObject({ id: "cull", label: "Cull", keys: ["x"], group: "Act" });
  });
});

describe("assertNoCoreCollision", () => {
  it("throws when a page verb steals a core key", () => {
    const bad = page("cull", { label: "Cull", keys: ["Enter"], group: "Act", run: () => {} });
    expect(() => assertNoCoreCollision([bad])).toThrow(/Enter/);
  });
  it("passes when keys are disjoint from core", () => {
    const ok = page("cull", { label: "Cull", keys: ["x"], group: "Act", run: () => {} });
    expect(() => assertNoCoreCollision([ok])).not.toThrow();
  });
});
