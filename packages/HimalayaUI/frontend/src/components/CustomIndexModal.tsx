import { useEffect, useState } from "react";
import { ModalShell } from "./ui/ModalShell";
import { SegmentedControl } from "./ui/SegmentedControl";
import { Button } from "./ui";
import { phaseColor } from "../phases";
import {
  SYMS, customRefls, basisFor, landsOn, snapTo, latticeForFirstOrderOnPeak,
} from "../lib/customIndex";

/**
 * CustomIndexModal (Plan D Task D-9) — build a custom index by symmetry +
 * lattice, with a live-preview comb computed client-side via real physics
 * (`2π√N/a` cubic, `2πn/d` lamellar, `4π√M/(√3·a)` hex), a running
 * "lands on N of M observed peaks" fit, and snap-to-peaks (magnetic slider drag
 * + click an observed peak to snap the first reflection onto it).
 *
 * Commit submits {phase, basis} where `basis` is the first reflection's q (the
 * q₁ slope) — the value `predicted_q_for_phase` reproduces (review finding #4).
 *
 * The lattice numeric input follows the QNumInput focus-gate pattern (finding
 * #6): external snap-drag updates the draft only when the input is NOT focused,
 * so a magnetic snap never interrupts mid-edit.
 *
 * Composes ui/ primitives (ModalShell + SegmentedControl + Button). Appearance
 * threads phaseColor() through SVG attributes (design-guard-safe).
 */

const VB = { l: 24, r: 16, t: 10, b: 26, w: 520, h: 120 };
const PLOT_W = VB.w - VB.l - VB.r;

const SYM_OPTIONS = [
  { value: "Pn3m" as const, label: "Pn3m" },
  { value: "Im3m" as const, label: "Im3m" },
  { value: "Ia3d" as const, label: "Ia3d" },
  { value: "Lamellar" as const, label: "Lamellar" },
  { value: "Hexagonal" as const, label: "Hexagonal" },
];

export interface CustomIndexModalProps {
  open: boolean;
  /** Observed (non-excluded) peak q-values, for the fit + snap targets. */
  peakQs: number[];
  onCommit: (phase: string, basis: number) => void;
  onClose: () => void;
}

export function CustomIndexModal({ open, peakQs, onCommit, onClose }: CustomIndexModalProps): JSX.Element | null {
  const [sym, setSym] = useState<string>("Pn3m");
  const [val, setVal] = useState<number>(SYMS.Pn3m!.def);
  const [draft, setDraft] = useState(String(SYMS.Pn3m!.def));
  const [focused, setFocused] = useState(false);

  // When switching symmetry, seed the lattice at that symmetry's default.
  const setSymmetry = (next: string): void => {
    setSym(next);
    const d = SYMS[next]?.def ?? val;
    setVal(d);
    if (!focused) setDraft(String(d));
  };

  // Focus-gate (finding #6): external (snap-drag) val changes flow into the
  // draft only when the input is not actively being edited.
  useEffect(() => {
    if (!focused) setDraft(String(Math.round(val)));
  }, [val, focused]);

  if (!open) return null;

  const spec = SYMS[sym]!;
  const refls = customRefls(sym, val);
  const color = phaseColor(sym);
  const observed = peakQs.filter((q) => q > 0);
  const lands = landsOn(sym, val, observed);

  // log-q domain over the union of predicted + observed q.
  const qs = [...refls.map((r) => r.q), ...observed].filter((q) => q > 0);
  const qMin = qs.length ? Math.min(...qs) * 0.92 : 0.01;
  const qMax = qs.length ? Math.max(...qs) * 1.08 : 0.2;
  const xOf = (q: number): number => {
    if (q <= 0 || qMin <= 0 || qMax <= qMin) return VB.l;
    return VB.l + ((Math.log(q) - Math.log(qMin)) / (Math.log(qMax) - Math.log(qMin))) * PLOT_W;
  };

  const commitVal = (raw: number): void => {
    const snapped = snapTo(sym, raw, observed);
    const clamped = Math.max(spec.min, Math.min(spec.max, snapped));
    setVal(clamped);
  };

  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="md"
      aria-label="Build a custom index"
      testId="custom-index-modal"
    >
      <div className="flex flex-col gap-4 p-5">
        <div className="text-headline">Custom index</div>

        <SegmentedControl
          options={SYM_OPTIONS}
          value={sym}
          onChange={setSymmetry}
          role="radiogroup"
          size="sm"
          aria-label="Symmetry"
          testId="custom-index-symmetry"
        />

        {/* lattice slider + focus-gated numeric input */}
        <div className="flex items-center gap-3">
          <label className="text-sm text-ink-soft w-4">{spec.param}</label>
          <input
            type="range"
            data-testid="custom-index-slider"
            min={spec.min}
            max={spec.max}
            step={1}
            value={val}
            onChange={(e) => commitVal(parseFloat(e.currentTarget.value))}
            className="flex-1 accent-accent"
            aria-label={`${spec.param} lattice`}
          />
          <input
            type="number"
            data-testid="custom-index-value"
            value={draft}
            step={1}
            onChange={(e) => setDraft(e.currentTarget.value)}
            onFocus={() => setFocused(true)}
            onBlur={(e) => {
              setFocused(false);
              const n = parseFloat(e.currentTarget.value);
              if (Number.isFinite(n)) commitVal(n);
            }}
            onKeyDown={(e) => {
              if (e.key === "Enter") {
                const n = parseFloat((e.currentTarget as HTMLInputElement).value);
                if (Number.isFinite(n)) commitVal(n);
                (e.currentTarget as HTMLInputElement).blur();
              }
            }}
            className="w-[72px] rounded border border-hair-strong bg-paper px-1 py-0.5
                       text-right text-xs tabular-nums text-ink outline-0 focus:border-accent"
            aria-label={`${spec.param} value`}
          />
          <span className="text-xs text-ink-faint">Å</span>
        </div>

        {/* live preview comb */}
        <svg
          data-testid="custom-index-preview"
          viewBox={`0 0 ${VB.w} ${VB.h}`}
          preserveAspectRatio="xMidYMid meet"
          className="w-full rounded-sm border border-hair bg-plate"
          aria-label="Custom-index preview comb"
        >
          {/* observed-peak guides (clickable: snap first order onto the peak) */}
          {observed.map((pq, i) => {
            const x = xOf(pq);
            return (
              <g key={`peak-${i}`} data-testid={`custom-index-peak-${i}`}
                 style={{ cursor: "pointer" }}
                 onClick={() => commitVal(latticeForFirstOrderOnPeak(sym, pq))}>
                <line x1={x} y1={VB.t} x2={x} y2={VB.h - VB.b}
                      stroke="var(--color-ink-faint)" strokeWidth={1} strokeDasharray="1 3" />
                <circle cx={x} cy={VB.h - VB.b} r={3} fill="none"
                        stroke="var(--color-ink-faint)" strokeWidth={1.4} />
                <circle cx={x} cy={(VB.t + VB.h - VB.b) / 2} r={8} fill="transparent" />
              </g>
            );
          })}
          {/* predicted teeth */}
          {refls.map((r) => {
            if (r.q < qMin || r.q > qMax) return null;
            const x = xOf(r.q);
            const onPeak = observed.some((pq) => Math.abs(pq - r.q) / r.q < 0.022);
            return (
              <g key={`refl-${r.N}`} data-testid={`custom-index-tooth-${r.N}`}
                 data-on-peak={onPeak || undefined}>
                <line x1={x} y1={VB.h - VB.b} x2={x} y2={VB.t + 6}
                      stroke={color} strokeWidth={onPeak ? 2.2 : 1.4} opacity={onPeak ? 1 : 0.6} />
                <circle cx={x} cy={VB.t + 6} r={2.4} fill={color} />
              </g>
            );
          })}
          <line x1={VB.l} y1={VB.h - VB.b} x2={VB.w - VB.r} y2={VB.h - VB.b}
                stroke="var(--color-hair-strong)" strokeWidth={1} />
        </svg>

        <div data-testid="custom-index-fit" className="text-xs text-ink-soft">
          Lands on <span className="font-mono font-bold text-ink">{lands}</span> of{" "}
          <span className="font-mono">{observed.length}</span> observed peaks ·{" "}
          {spec.param} = <span className="font-mono font-bold">{Math.round(val)} Å</span>
        </div>

        <div className="flex justify-end gap-2">
          <Button variant="ghost" onClick={onClose}>Cancel</Button>
          <Button
            variant="accent"
            data-testid="custom-index-commit"
            onClick={() => { onCommit(sym, basisFor(sym, val)); onClose(); }}
          >
            Add custom index
          </Button>
        </div>
      </div>
    </ModalShell>
  );
}
