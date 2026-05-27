import { describe, it, expect } from "vitest";
import { parseLocation } from "../src/lib/url/parseLocation";

// Spec §4.1 — discriminated-union output keyed by `kind`.

describe("parseLocation", () => {
  it("/ → root", () => {
    expect(parseLocation("/", "")).toEqual({ kind: "root" });
    expect(parseLocation("", "")).toEqual({ kind: "root" });
  });

  it("/index → empty index", () => {
    expect(parseLocation("/index", "")).toEqual({
      kind: "index", experiment: undefined, sample: undefined,
    });
  });

  it("/index/<exp> → experiment chosen", () => {
    expect(parseLocation("/index/lipid-screen", "")).toEqual({
      kind: "index", experiment: "lipid-screen", sample: undefined,
    });
  });

  it("/index/<exp>/<sample> → full", () => {
    expect(parseLocation("/index/lipid-screen/JC001", "")).toEqual({
      kind: "index", experiment: "lipid-screen", sample: "JC001",
    });
  });

  // I1.7 (#163): Inspect is retired. /inspect* is redirected to /samples at
  // the router layer, so parseLocation never legitimately sees it; the parse
  // arm is gone and any /inspect path now falls through to `stale`.
  it("/inspect is no longer parsed (falls through to stale)", () => {
    expect(parseLocation("/inspect/lipid/JC001", "?exposure=JC001-007")).toEqual({
      kind: "stale", raw: "/inspect/lipid/JC001?exposure=JC001-007",
    });
  });

  it("/compare (legacy list) → compare:list", () => {
    expect(parseLocation("/compare", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all → compare:list", () => {
    expect(parseLocation("/compare/all", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all/42 → compare:list (ComparePage handles numeric id)", () => {
    expect(parseLocation("/compare/all/42", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/compare/all/42/edit → compare:list", () => {
    expect(parseLocation("/compare/all/42/edit", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/experiments/17/compare → compare:list", () => {
    expect(parseLocation("/experiments/17/compare", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/experiments/17/compare/42 → compare:list", () => {
    expect(parseLocation("/experiments/17/compare/42", "")).toEqual({ kind: "compare", view: "list" });
  });

  it("/foo/bar → stale", () => {
    const r = parseLocation("/foo/bar", "?x=1");
    expect(r.kind).toBe("stale");
    if (r.kind === "stale") expect(r.raw).toBe("/foo/bar?x=1");
  });

  it("decodes percent-encoded slugs", () => {
    expect(parseLocation("/index/lipid%20screen/JC%20001", "")).toEqual({
      kind: "index", experiment: "lipid screen", sample: "JC 001",
    });
  });

  it("trailing slash tolerant", () => {
    expect(parseLocation("/index/", "")).toEqual({
      kind: "index", experiment: undefined, sample: undefined,
    });
  });
});
