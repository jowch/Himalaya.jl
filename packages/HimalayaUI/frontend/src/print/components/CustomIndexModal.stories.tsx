import { useMemo, useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { CustomIndexModal, type CustomIndexFit } from "./CustomIndexModal";
import type { CombSeries } from "../comb";
import { phaseColor } from "../../phases";

/* ── Reference physics (the SYMS table + tooth/snap math from the mockup
 *    `2026-05-29-focus-plot.html`). In the real app this lives in the Focus
 *    PAGE — the modal is presentational and only renders the result. The story
 *    reproduces it so the slider actually drives the preview + fit line, the
 *    way the page eventually will. q(hkl) = 2π·√N / a, N = h²+k²+l². ─────── */
const TWO_PI = 2 * Math.PI;
const Q_WINDOW = { min: 0.02, max: 0.23 };

type SymKind = "cubic" | "lamellar" | "hex";
interface Sym {
  kind: SymKind;
  Ms: number[];
  param: string;
  def: number;
  min: number;
  max: number;
}

const SYMS: Record<string, Sym> = {
  Pn3m: { kind: "cubic", Ms: [2, 3, 4, 6, 8, 9], param: "a", def: 197, min: 120, max: 320 },
  Im3m: { kind: "cubic", Ms: [2, 4, 6, 8, 10, 12], param: "a", def: 252, min: 120, max: 360 },
  Ia3d: { kind: "cubic", Ms: [6, 8, 14, 16, 20, 22], param: "a", def: 218, min: 120, max: 360 },
  Lamellar: { kind: "lamellar", Ms: [1, 2, 3, 4, 5], param: "d", def: 60, min: 30, max: 130 },
  Hexagonal: { kind: "hex", Ms: [1, 3, 4, 7, 9, 12], param: "a", def: 70, min: 40, max: 160 },
};

const SYMMETRIES = ["Pn3m", "Im3m", "Ia3d", "Lamellar", "Hexagonal"] as const;

/* Observed peaks for smp_09 (a dominant Pn3m + a subtle Im3m), computed from
 * the real lattices with the mockup's sub-percent fit residuals. */
const qOf = (a: number, n: number): number => (TWO_PI * Math.sqrt(n)) / a;
const OBSERVED: number[] = [
  qOf(197, 2) * 1.0, // Pn3m √2
  qOf(197, 3) * 1.004, // Pn3m √3
  qOf(197, 4) * 0.997, // Pn3m √4
  qOf(252, 4) * 1.006, // Im3m √4 (subtle)
  qOf(252, 6) * 1.014, // Im3m √6 (hand-added, deviates)
];

/** Predicted reflections for one symmetry at a lattice value. */
function customRefls(sym: string, val: number): Array<{ N: number; q: number }> {
  const s = SYMS[sym]!;
  if (s.kind === "lamellar") return s.Ms.map((n) => ({ N: n, q: (TWO_PI * n) / val }));
  if (s.kind === "hex") {
    const q1 = (4 * Math.PI) / (Math.sqrt(3) * val);
    return s.Ms.map((M) => ({ N: M, q: q1 * Math.sqrt(M) }));
  }
  return s.Ms.map((N) => ({ N, q: (TWO_PI * Math.sqrt(N)) / val }));
}

/** Lattice values at which any allowed order lands exactly on an observed peak. */
function snapValues(sym: string): number[] {
  const s = SYMS[sym]!;
  const vals: number[] = [];
  OBSERVED.forEach((pq) =>
    s.Ms.forEach((M) => {
      const v =
        s.kind === "lamellar"
          ? (TWO_PI * M) / pq
          : s.kind === "hex"
            ? (4 * Math.PI * Math.sqrt(M)) / (Math.sqrt(3) * pq)
            : (TWO_PI * Math.sqrt(M)) / pq;
      if (v >= s.min && v <= s.max) vals.push(v);
    }),
  );
  return vals;
}
const isSnapped = (sym: string, v: number): boolean =>
  snapValues(sym).some((s) => Math.abs(s - v) < 0.75);

/** Build the live preview series + fit from a symmetry + lattice value — the
 *  computation the Focus page will own. A tooth is "observed" when a real peak
 *  sits within 2.2% of its predicted q. */
function buildPreview(sym: string, val: number): { series: CombSeries; fit: CustomIndexFit } {
  const teeth = customRefls(sym, val)
    .filter((rf) => rf.q >= Q_WINDOW.min && rf.q <= Q_WINDOW.max)
    .map((rf) => ({
      q: rf.q,
      label: `√${rf.N}`,
      observed: OBSERVED.some((pq) => Math.abs(pq - rf.q) / rf.q < 0.022),
    }));
  const landed = teeth.filter((t) => t.observed).length;
  return {
    series: { phase: sym, color: phaseColor(sym), teeth },
    fit: { landed, total: OBSERVED.length, snapped: isSnapped(sym, val) },
  };
}

const noop = (): void => {};
const initial = buildPreview("Im3m", SYMS["Im3m"]!.def);

const meta = {
  title: "components/CustomIndexModal",
  component: CustomIndexModal,
  // Required props — overridden by the interactive render below.
  args: {
    open: true,
    onClose: noop,
    onCancel: noop,
    onAdd: noop,
    symmetries: SYMMETRIES,
    symmetry: "Im3m",
    paramName: "a",
    paramValue: String(SYMS["Im3m"]!.def),
    paramMin: 120,
    paramMax: 360,
    paramStep: 1,
    onSymmetryChange: noop,
    onParamChange: noop,
    previewSeries: initial.series,
    observed: OBSERVED,
    fit: initial.fit,
  },
} satisfies Meta<typeof CustomIndexModal>;

export default meta;
type Story = StoryObj<typeof meta>;

function OpenDemo(): JSX.Element {
  const [open, setOpen] = useState(true);
  const [symmetry, setSymmetry] = useState<string>("Im3m");
  const [paramValue, setParamValue] = useState<string>(String(SYMS["Im3m"]!.def));

  const cfg = SYMS[symmetry]!;

  // The page-owned computation: predicted teeth + fit recompute on every
  // symmetry / lattice change, so dragging the slider moves the teeth and the
  // "Lands on N of M" line is real, not hard-coded.
  const preview = useMemo(() => buildPreview(symmetry, Number(paramValue)), [symmetry, paramValue]);

  const handleSymmetryChange = (s: string): void => {
    setSymmetry(s);
    setParamValue(String(SYMS[s]!.def));
  };

  return (
    <div>
      <button onClick={() => setOpen(true)}>Open custom index</button>
      <CustomIndexModal
        open={open}
        onClose={() => setOpen(false)}
        onCancel={() => setOpen(false)}
        onAdd={() => {
          setOpen(false);
        }}
        symmetries={SYMMETRIES}
        symmetry={symmetry}
        onSymmetryChange={handleSymmetryChange}
        paramName={cfg.param}
        paramValue={paramValue}
        paramMin={cfg.min}
        paramMax={cfg.max}
        paramStep={1}
        onParamChange={setParamValue}
        unit="Å"
        previewSeries={preview.series}
        observed={OBSERVED}
        fit={preview.fit}
      />
    </div>
  );
}

export const Open: Story = {
  render: () => <OpenDemo />,
};
