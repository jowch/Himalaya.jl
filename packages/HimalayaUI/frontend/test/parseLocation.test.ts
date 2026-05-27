import { describe, it, expect } from "vitest";
import { parseLocation } from "../src/lib/url/parseLocation";

// Spec §4.1 — discriminated-union output keyed by `kind`.

describe("parseLocation", () => {
  it("/ → root", () => {
    expect(parseLocation("/", "")).toEqual({ kind: "root" });
    expect(parseLocation("", "")).toEqual({ kind: "root" });
  });

  // I4.4 (#181): Index is retired. /index* is redirected (to /samples or, for
  // a slug-bearing /index/<exp>/<sample>, to /sample/:id via IndexSlugRedirect)
  // at the router layer, so parseLocation never legitimately sees it; the parse
  // arm is gone and any /index path now falls through to `stale`.
  it("/index is no longer parsed (falls through to stale)", () => {
    expect(parseLocation("/index", "")).toEqual({
      kind: "stale", raw: "/index",
    });
  });

  it("/index/<exp> is no longer parsed (falls through to stale)", () => {
    expect(parseLocation("/index/lipid-screen", "")).toEqual({
      kind: "stale", raw: "/index/lipid-screen",
    });
  });

  it("/index/<exp>/<sample> is no longer parsed (falls through to stale)", () => {
    expect(parseLocation("/index/lipid-screen/JC001", "")).toEqual({
      kind: "stale", raw: "/index/lipid-screen/JC001",
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

  // I3.6 (#177): Compare is retired. /compare* (both the global /compare/all
  // root and the experiment-scoped /experiments/:eid/compare root, incl.
  // /new, /:id, /:id/edit) is redirected to /series at the router layer, so
  // parseLocation never legitimately sees it; the parse arm is gone and any
  // /compare path now falls through to `stale`.
  it("/compare (legacy list) → stale", () => {
    expect(parseLocation("/compare", "")).toEqual({ kind: "stale", raw: "/compare" });
  });

  it("/compare/all → stale", () => {
    expect(parseLocation("/compare/all", "")).toEqual({ kind: "stale", raw: "/compare/all" });
  });

  it("/compare/all/42 → stale", () => {
    expect(parseLocation("/compare/all/42", "")).toEqual({ kind: "stale", raw: "/compare/all/42" });
  });

  it("/compare/all/42/edit → stale", () => {
    expect(parseLocation("/compare/all/42/edit", "")).toEqual({ kind: "stale", raw: "/compare/all/42/edit" });
  });

  it("/experiments/17/compare → stale", () => {
    expect(parseLocation("/experiments/17/compare", "")).toEqual({ kind: "stale", raw: "/experiments/17/compare" });
  });

  it("/experiments/17/compare/42 → stale", () => {
    expect(parseLocation("/experiments/17/compare/42", "")).toEqual({ kind: "stale", raw: "/experiments/17/compare/42" });
  });

  it("/foo/bar → stale", () => {
    const r = parseLocation("/foo/bar", "?x=1");
    expect(r.kind).toBe("stale");
    if (r.kind === "stale") expect(r.raw).toBe("/foo/bar?x=1");
  });

  it("trailing slash tolerant (/index/ → stale, normalized)", () => {
    expect(parseLocation("/index/", "")).toEqual({
      kind: "stale", raw: "/index/",
    });
  });
});
