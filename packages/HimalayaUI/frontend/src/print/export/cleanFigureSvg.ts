// cleanFigureSvg.ts — build the exported figure as a standalone SVG string,
// using the SAME greenfield d3 axis math (makeAxis/axisTicks) that the
// on-screen WaterfallChart uses. This is the export idiom: a plain, framed
// GraphPad/journal figure (white ground, Arial, literal hex, L-shaped black
// axes) — deliberately NOT "The Print" branding.
//
// Identity vs colour (the export's deliberate divergence from the plate): a
// trace is named in the KEY BLOCK below the plot (swatch + sample + each phase
// with its lattice a/d and, for bicontinuous cubics, κ). Colour identifies the
// TRACE, not the phase — the series walks a distinguishable categorical palette
// so multiple traces (possibly sharing a phase, or coexisting) stay separable,
// and phase is carried as text. The single-trace Focus figure passes one
// phase-coloured key, so its key block reads as the phase + lattice legend.
import type { WaterfallRow } from "../waterfall/waterfallModel";
import { makeAxis, axisTicks } from "../plot/projection";
import { waterfallQDomain } from "../waterfall/waterfallModel";
import { phaseHex } from "./traceColors";

const ARIAL = "Arial, Helvetica, sans-serif";
const INK = "#111111";
const INK_MUTED = "#666666";

/** One phase a trace carries, pre-formatted for the key block. `detail` is the
 *  ready-to-render lattice (+ κ) string, e.g. "a = 252 Å · κ = 1.70×10⁻⁴ Å⁻²";
 *  empty when no lattice is known. */
export interface FigureKeySegment {
  phase: string;
  detail: string;
}

/** One trace's key-block entry: a colour swatch, a primary label (sample name;
 *  empty on the single-trace figure, where the title already carries it), and
 *  the phases it carries. `note` covers the phaseless states (e.g. "unindexed",
 *  "form factor (no Bragg)"). */
export interface FigureTraceKey {
  color: string;
  label: string;
  segments: FigureKeySegment[];
  note?: string;
}

export interface CleanFigureInput {
  /** Rows low→high (display order); rendered bottom-up, mirroring the plate. */
  rows: WaterfallRow[];
  title: string;
  footer: string;
  /** Per-trace key entries, parallel to `rows` (same index). When present, a
   *  trace is coloured by its key colour and the key block replaces the phase
   *  legend. When absent, traces fall back to phase colour (legacy path). */
  traceKeys?: FigureTraceKey[];
  /** q-axis scale, mirroring the plate's log/linear toggle. Default "log". */
  xType?: "log" | "linear";
  /** Plate trace-offset: scales the inter-row separation. Default 1. */
  offsetScale?: number;
  showPeakTicks?: boolean;
  showPeakLabels?: boolean;
  /** Visible q-domain (the plate's). Falls back to the rows' positive extent. */
  qDomain?: [number, number];
  width?: number;
  height?: number;
}

// Default vertical slot per trace. With the 480px default width and the fixed
// margins below, three traces land at a narrow ~1.85:1 PORTRAIT (a landscape
// waterfall squished the traces flat); the figure grows taller with each trace.
const PER_ROW_H = 248;
const BOTTOM_GAP = 14;
// Key-block line metrics.
const KEY_HEADER_H = 18;
const KEY_SUBLINE_H = 16;
const KEY_TOP_PAD = 20; // gap between the x-axis title and the first key row
const FOOTNOTE_H = 22;
const AXIS_FOOT_H = 40; // tick labels + x-axis title below the baseline

function esc(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

/** Flatten one trace key into its rendered lines. The header line carries the
 *  swatch; with a label it leads, otherwise the first phase segment rides the
 *  header (the single-trace legend register). */
interface KeyLine { swatch: boolean; phase: string; text: string; muted: boolean }
function keyLines(k: FigureTraceKey): KeyLine[] {
  const segLine = (s: FigureKeySegment): Omit<KeyLine, "swatch"> => ({
    phase: s.phase,
    text: s.detail,
    muted: false,
  });
  const lines: KeyLine[] = [];
  if (k.label) {
    lines.push({ swatch: true, phase: "", text: k.label, muted: false });
    for (const s of k.segments) lines.push({ swatch: false, ...segLine(s) });
    if (k.segments.length === 0 && k.note) {
      lines.push({ swatch: false, phase: "", text: k.note, muted: true });
    }
  } else if (k.segments.length > 0) {
    lines.push({ swatch: true, ...segLine(k.segments[0]!) });
    for (const s of k.segments.slice(1)) lines.push({ swatch: false, ...segLine(s) });
  } else {
    lines.push({ swatch: true, phase: "", text: k.note ?? "", muted: true });
  }
  return lines;
}

/** Build the export figure as a standalone `<svg>` markup string. */
export function buildCleanFigureSvg(input: CleanFigureInput): string {
  const {
    rows,
    title,
    footer,
    traceKeys,
    xType = "log",
    offsetScale = 1,
    showPeakTicks = true,
    showPeakLabels = false,
  } = input;
  const width = input.width ?? 480; // narrow portrait (see PER_ROW_H)

  // Resolve each trace's colour: key colour when present, else phase hue.
  const colorOf = (i: number, row: WaterfallRow): string =>
    traceKeys?.[i]?.color ?? phaseHex(row.phase);

  // Lay out the key block first — its height grows with trace + phase count, so
  // it sizes the bottom margin (and thus the figure height). Rendered TOP trace
  // first (reverse of display order) so the key's vertical order matches the
  // bottom-up stacked plot — the eye reads the top key row against the top trace.
  const keyGroups = (traceKeys ?? [])
    .map((k) => ({ color: k.color, lines: keyLines(k) }))
    .reverse();
  const keyBlockH =
    keyGroups.length > 0
      ? KEY_TOP_PAD + keyGroups.length * 4 /* inter-entry gap */ +
        keyGroups.reduce(
          (s, g) => s + KEY_HEADER_H + Math.max(0, g.lines.length - 1) * KEY_SUBLINE_H,
          0,
        )
      : 0;

  const m = {
    l: 58,
    r: 26, // no right gutter — traces are named in the key block
    t: 36,
    b: AXIS_FOOT_H + keyBlockH + FOOTNOTE_H,
  };
  const height =
    input.height ?? m.t + m.b + BOTTOM_GAP + Math.max(1, rows.length) * PER_ROW_H;
  const pw = width - m.l - m.r;
  const ph = height - m.t - m.b;
  const baseY = m.t + ph;
  const scale = Math.max(0.1, offsetScale);

  const qDomain = input.qDomain ?? waterfallQDomain(rows);
  const x = makeAxis(qDomain, [m.l, m.l + pw], xType);
  const inDomain = (q: number): boolean => q >= qDomain[0] && q <= qDomain[1];
  const ticks = axisTicks(x);

  // Band geometry mirrors WaterfallChart: bottom-up (row 0 at the bottom),
  // bandHeight-weighted, the inter-row STEP scaled by the offset, fit to panel.
  const stackH = ph - BOTTOM_GAP;
  const totalWeight = rows.reduce((s, r) => s + Math.max(0, r.bandHeight), 0) || rows.length;
  let cumulative = 0;
  let natural = 0;
  const nat = rows
    .map((row, idx) => ({ row, idx }))
    .reverse()
    .map(({ row, idx }) => {
      const h = (Math.max(0, row.bandHeight) / totalWeight) * stackH;
      const top = cumulative + row.yOffset;
      cumulative += h * scale;
      natural = Math.max(natural, top + h);
      return { row, idx, top, h };
    });
  const fit = natural > 0 ? stackH / natural : 1;
  const bands = nat.map(({ row, idx, top, h }) => {
    const bandTop = m.t + top * fit;
    const bandBottom = bandTop + h * fit;
    return { row, idx, bandTop, bandBottom };
  });

  const parts: string[] = [];
  parts.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" ` +
      `viewBox="0 0 ${width} ${height}" font-family="${ARIAL}">`,
  );
  parts.push(`<rect x="0" y="0" width="${width}" height="${height}" fill="#ffffff"/>`);

  // Title (top-left, like a figure caption header).
  parts.push(
    `<text x="${m.l}" y="24" font-size="16.5" font-weight="700" fill="${INK}">${esc(title)}</text>`,
  );

  // L-shaped black axes.
  parts.push(
    `<line x1="${m.l}" y1="${baseY}" x2="${m.l + pw}" y2="${baseY}" stroke="${INK}" stroke-width="1.4"/>`,
  );
  parts.push(
    `<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${baseY}" stroke="${INK}" stroke-width="1.4"/>`,
  );

  // X ticks + labels (labels only on non-minor ticks).
  for (const t of ticks) {
    if (!inDomain(t.value)) continue;
    const px = x.to(t.value);
    parts.push(
      `<line x1="${px.toFixed(1)}" y1="${baseY}" x2="${px.toFixed(1)}" y2="${baseY + 5}" stroke="${INK}" stroke-width="1"/>`,
    );
    if (t.kind !== "minor") {
      parts.push(
        `<text x="${px.toFixed(1)}" y="${baseY + 18}" text-anchor="middle" font-size="11.5" fill="${INK}">${t.value.toFixed(2)}</text>`,
      );
    }
  }

  // Axis titles.
  parts.push(
    `<text x="${m.l + pw / 2}" y="${baseY + 36}" text-anchor="middle" font-size="13" font-weight="700" fill="${INK}">q (Å⁻¹)</text>`,
  );
  const yMid = m.t + ph / 2;
  parts.push(
    `<text x="16" y="${yMid}" text-anchor="middle" font-size="13" font-weight="700" fill="${INK}" transform="rotate(-90 16 ${yMid})">Intensity (a.u.) + offset</text>`,
  );

  // Traces (log-intensity per band) + peak glyphs.
  for (const { row, idx, bandTop, bandBottom } of bands) {
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = colorOf(idx, row);
    const { q, I } = row.trace;

    // Per-row log-intensity mapping: map log10(I) ∈ [logMin, logMax] → y ∈
    // [bandBottom, bandTop] (low at bottom). Hoisted so peak markers ride the
    // SAME curve as the trace.
    let minI = Infinity;
    let maxI = -Infinity;
    for (const v of I) {
      if (Number.isFinite(v) && v > 0) {
        if (v < minI) minI = v;
        if (v > maxI) maxI = v;
      }
    }
    const hasLog = maxI > 0;
    const logMin = hasLog ? Math.log10(minI) : 0;
    const logRange = hasLog ? Math.max(1e-6, Math.log10(maxI) - logMin) : 1;
    const yOf = (v: number): number =>
      bandBottom - (((v > 0 ? Math.log10(v) : logMin) - logMin) / logRange) * bandH;

    if (hasLog) {
      let d = "";
      let started = false;
      for (let k = 0; k < q.length; k++) {
        const qk = q[k]!;
        if (!inDomain(qk)) continue;
        d += `${started ? "L" : "M"}${x.to(qk).toFixed(1)} ${yOf(I[k] ?? 0).toFixed(1)} `;
        started = true;
      }
      if (d) {
        parts.push(`<path d="${d.trim()}" fill="none" stroke="${color}" stroke-width="2" data-phase="${row.phase ?? "none"}"/>`);
      }
    }

    // Peak anchors — ordinal-numbered (1..n by ascending q). Each marker rides
    // the trace: a small downward triangle pointing at the curve, in the trace
    // colour (we don't attribute a peak to a phase under coexistence).
    if (showPeakTicks && hasLog && row.anchors.length > 0) {
      const sorted = [...row.anchors].filter((a) => inDomain(a.q)).sort((a, b) => a.q - b.q);
      sorted.forEach((a, i) => {
        const px = x.to(a.q);
        const cy = yOf(a.intensity ?? 0);
        const tip = cy - 4;
        const top = tip - 5;
        parts.push(
          `<path d="M${(px - 3.2).toFixed(1)} ${top.toFixed(1)} L${(px + 3.2).toFixed(1)} ${top.toFixed(1)} L${px.toFixed(1)} ${tip.toFixed(1)} Z" fill="${color}"/>`,
        );
        if (showPeakLabels) {
          parts.push(
            `<text x="${px.toFixed(1)}" y="${(top - 2).toFixed(1)}" text-anchor="middle" font-size="9.5" fill="${color}">${i + 1}</text>`,
          );
        }
      });
    }
  }

  // ── Key block (replaces the phase legend + per-row gutter labels) ──────────
  if (keyGroups.length > 0) {
    let ky = baseY + AXIS_FOOT_H + KEY_TOP_PAD;
    for (const g of keyGroups) {
      g.lines.forEach((ln, li) => {
        const isHeader = li === 0;
        if (ln.swatch) {
          parts.push(
            `<rect x="${m.l}" y="${(ky - 9).toFixed(1)}" width="11" height="11" rx="1.5" fill="${g.color}"/>`,
          );
        }
        const tx = m.l + 18;
        const fill = ln.muted ? INK_MUTED : INK;
        if (ln.phase) {
          // "<Phase>   <detail>" — phase name in bold, lattice/κ detail regular.
          parts.push(
            `<text x="${tx}" y="${ky.toFixed(1)}" font-size="11.5" fill="${fill}"><tspan font-weight="700">${esc(ln.phase)}</tspan>${ln.text ? `  ${esc(ln.text)}` : ""}</text>`,
          );
        } else {
          parts.push(
            `<text x="${tx}" y="${ky.toFixed(1)}" font-size="11.5" font-weight="${isHeader ? 700 : 400}" fill="${fill}">${esc(ln.text)}</text>`,
          );
        }
        ky += isHeader ? KEY_HEADER_H : KEY_SUBLINE_H;
      });
      ky += 4; // inter-entry gap
    }
  }

  // Footnote, at the very bottom.
  parts.push(
    `<text x="${(width / 2).toFixed(1)}" y="${(height - 8).toFixed(1)}" text-anchor="middle" font-size="10.5" fill="${INK_MUTED}">${esc(footer)}</text>`,
  );

  parts.push("</svg>");
  return parts.join("");
}
