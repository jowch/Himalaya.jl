import { describe, it, expect, beforeEach } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { PLOT_VIEWS, usePlotView, viewAxes, viewIntensity } from "../../src/workspace/plotView";

describe("plot views", () => {
  beforeEach(() => window.localStorage.clear());

  it("offers log, linear, log-log and q²I in that order", () => {
    expect(PLOT_VIEWS.map((v) => v.view)).toEqual(["log", "linear", "loglog", "kratky"]);
  });

  it("maps each view to its axis scales", () => {
    expect(viewAxes("log")).toEqual({ xType: "linear", yType: "log" });
    expect(viewAxes("linear")).toEqual({ xType: "linear", yType: "linear" });
    expect(viewAxes("loglog")).toEqual({ xType: "log", yType: "log" });
    expect(viewAxes("kratky")).toEqual({ xType: "linear", yType: "linear" });
  });

  it("multiplies intensity by q² only in the q²I view", () => {
    expect(viewIntensity("kratky", 0.2, 50)).toBeCloseTo(2);
    for (const v of ["log", "linear", "loglog"] as const) expect(viewIntensity(v, 0.2, 50)).toBe(50);
  });

  it("remembers the chosen view and falls back to log for anything unknown", () => {
    const first = renderHook(() => usePlotView());
    expect(first.result.current[0]).toBe("log");
    act(() => first.result.current[1]("kratky"));
    expect(first.result.current[0]).toBe("kratky");
    first.unmount();
    expect(renderHook(() => usePlotView()).result.current[0]).toBe("kratky");
    window.localStorage.setItem("himalaya-ui:workspace-plot-view", "polar");
    expect(renderHook(() => usePlotView()).result.current[0]).toBe("log");
  });
});
