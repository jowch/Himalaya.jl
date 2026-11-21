import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { render, screen } from "@testing-library/react";
import { PLOT_HEIGHT, WorkspacePlot } from "../../src/workspace/ui/WorkspacePlot";
import { PLOT_MARGINS } from "../../src/workspace/plotScale";
import type { PlotView } from "../../src/workspace/plotView";

function fixture(name: string): { q: number[]; I: number[]; sigma: number[] } {
  const rows = readFileSync(resolve(__dirname, "../../../../../test/data", `${name}_tot.dat`), "utf8")
    .trim()
    .split("\n")
    .map((l) => l.trim().split(/\s+/).map(Number));
  return { q: rows.map((r) => r[0]!), I: rows.map((r) => r[1]!), sigma: rows.map((r) => r[2]!) };
}

const plotHeight = PLOT_HEIGHT - PLOT_MARGINS.top - PLOT_MARGINS.bottom;
const VIEWS: PlotView[] = ["log", "linear", "loglog", "kratky"];
const WINDOWS: ([number, number] | null)[] = [null, [0.0075, 0.0085], [0.02, 0.05], [0.08, 0.15], [0.12, 0.16]];

const labels = (role: string): string[] =>
  Array.from(screen.getByTestId("workspace-trace").querySelectorAll(`[data-role="${role}"] text`)).map((t) => t.textContent ?? "");
const ys = (selector: string): number[] =>
  (screen.getByTestId("workspace-trace").querySelector(selector)?.getAttribute("d") ?? "")
    .split(/[ML]/)
    .filter(Boolean)
    .map((p) => p.split(",").map(Number) as [number, number])
    .filter(([x]) => x >= 0 && x <= 640 - PLOT_MARGINS.left - PLOT_MARGINS.right)
    .map(([, y]) => y);

describe("WorkspacePlot on real SSRL traces", () => {
  for (const name of ["cubic", "example", "form-factor"]) {
    it(`${name}: every view and zoom has at least two distinct axis labels on each axis, and the trace inside the plot`, () => {
      const t = fixture(name);
      for (const view of VIEWS) {
        for (const w of WINDOWS) {
          if (w && (w[0] < t.q[0]! || w[1] > t.q[t.q.length - 1]!)) continue;
          const { unmount } = render(<WorkspacePlot trace={t} peaks={[]} combs={[]} view={view} xDomain={w} />);
          for (const role of ["axis-left", "axis-bottom"]) {
            const ls = labels(role);
            expect(ls.length, `${view} ${w} ${role}`).toBeGreaterThanOrEqual(2);
            expect(new Set(ls).size, `${view} ${w} ${role} ${ls.join("|")}`).toBe(ls.length);
          }
          const line = ys('[data-role="trace-line"] path');
          expect(Math.min(...line), `${view} ${w}`).toBeGreaterThan(0);
          expect(Math.max(...line), `${view} ${w}`).toBeLessThan(plotHeight);
          unmount();
        }
      }
    });
  }

  it("keeps a stronger comparison inside the plot in q²I", () => {
    const main = fixture("form-factor");
    const ref = fixture("cubic");
    const compare = { id: 9, name: "cubic", q: ref.q, I: ref.I.map((v) => v * 3), dash: "5 4", emphasis: false };
    render(<WorkspacePlot trace={main} peaks={[]} combs={[]} compares={[compare]} view="kratky" xDomain={[0.12, 0.16]} />);
    const c = ys('[data-compare="9"]');
    expect(Math.min(...c)).toBeGreaterThanOrEqual(0);
  });

  it("does not stretch the y axis for a comparison that ends before the window", () => {
    const main = fixture("example");
    const alone = (): number => {
      const { unmount } = render(<WorkspacePlot trace={main} peaks={[]} combs={[]} view="linear" xDomain={[0.1, 0.15]} />);
      const top = labels("axis-left").at(-1)!;
      unmount();
      return Number(top);
    };
    const base = alone();
    const early = { id: 4, name: "early", q: main.q.filter((q) => q < 0.05), I: main.I.filter((_, i) => main.q[i]! < 0.05).map((v) => v * 10), dash: "5 4", emphasis: false };
    render(<WorkspacePlot trace={main} peaks={[]} combs={[]} compares={[early]} view="linear" xDomain={[0.1, 0.15]} />);
    expect(Number(labels("axis-left").at(-1))).toBe(base);
  });
});
