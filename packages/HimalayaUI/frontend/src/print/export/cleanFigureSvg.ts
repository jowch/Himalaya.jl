// cleanFigureSvg.ts — build the exported figure as a standalone SVG string,
// using the SAME greenfield d3 axis math (makeAxis/axisTicks) that the
// on-screen WaterfallChart uses. This is the export idiom: a plain, framed
// GraphPad/journal figure (white ground, Arial, literal hex, L-shaped black
// axes) — deliberately NOT "The Print" branding.
//
// UNITS: the figure is authored in POINTS (72 pt/in) at a compact physical size
// (~2 in wide), with point-sized fonts (axis titles 10.5 pt, ticks 8 pt). The
// SVG download is therefore the true physical figure; the PNG raster upscales
// (raster.ts, dpi/72) only at export time. So sizes here are read as points.
//
// VERTICAL PACKING: traces are stacked with a fixed amplitude and a SMALLER
// step, so each trace's tall low-q head overlaps a little into the trace above
// (they don't share data — the stack is offset — so the overlap reads fine and
// uses the space instead of wasting it). The plate's offset slider scales the
// step. NOT a fit-to-panel model (that wasted height keeping traces apart).
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

// Font sizes, in points.
const FS_TITLE = 12;
const FS_AXIS = 10.5;
const FS_TICK = 8;
const FS_KEY = 8.5; // sample name / phase name (bold)
const FS_DETAIL = 8; // lattice + κ detail / note
const FS_FOOT = 7;
const FS_ORD = 6.5; // peak ordinal labels

// Vertical packing, in points. A generous amplitude with a step ~75% of it
// keeps the plot PANEL portrait (taller than wide) with only a MILD overlap of
// each trace's tall low-q head into the trace above — enough to use the space,
// not so much that peak labels collide with the trace below.
const AMP = 82; // amplitude of the tallest band (peak-to-floor of one trace)
const STEP_BASE = 62; // baseline-to-baseline step (< AMP ⇒ ~25% overlap)
const STACK_GAP = 6; // lift row 0 off the x-axis
const TOP_HEAD = 6; // headroom above the tallest trace

// Key block + footnote metrics, in points.
const KEY_HEADER_H = 12;
const KEY_SUBLINE_H = 10.5;
const KEY_TOP_PAD = 13; // gap below the x-axis title to the first key row
const KEY_ENTRY_GAP = 4;
const KEY_SWATCH = 8;
const KEY_INDENT = 13; // text inset from the swatch
const FOOTNOTE_H = 14;
const AXIS_FOOT_H = 26; // tick labels + x-axis title below the baseline

const MIN_WIDTH = 132; // ~1.85 in — figure never narrower than this
const CHAR_W = 0.52; // Arial average advance, as a fraction of font size

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

function esc(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v));
}

/** A rendered key line: header carries the swatch; phase lines carry a bold
 *  phase token + a regular detail; a plain line is a label or a muted note. */
interface KeyLine { swatch: boolean; phase: string; text: string; muted: boolean }
function keyLines(k: FigureTraceKey): KeyLine[] {
  const seg = (s: FigureKeySegment): Omit<KeyLine, "swatch"> => ({
    phase: s.phase,
    text: s.detail,
    muted: false,
  });
  const lines: KeyLine[] = [];
  if (k.label) {
    lines.push({ swatch: true, phase: "", text: k.label, muted: false });
    for (const s of k.segments) lines.push({ swatch: false, ...seg(s) });
    if (k.segments.length === 0 && k.note) {
      lines.push({ swatch: false, phase: "", text: k.note, muted: true });
    }
  } else if (k.segments.length > 0) {
    lines.push({ swatch: true, ...seg(k.segments[0]!) });
    for (const s of k.segments.slice(1)) lines.push({ swatch: false, ...seg(s) });
  } else {
    lines.push({ swatch: true, phase: "", text: k.note ?? "", muted: true });
  }
  return lines;
}

/** Estimate a rendered key line's text width (points). */
function keyLineWidth(ln: KeyLine): number {
  const full = ln.phase ? `${ln.phase}  ${ln.text}` : ln.text;
  return full.length * FS_KEY * CHAR_W;
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

  // Resolve each trace's colour: key colour when present, else phase hue.
  const colorOf = (i: number, row: WaterfallRow): string =>
    traceKeys?.[i]?.color ?? phaseHex(row.phase);

  // Key block: render TOP trace first (reverse of display order) so the key's
  // vertical order matches the bottom-up stacked plot.
  const keyGroups = (traceKeys ?? [])
    .map((k) => ({ color: k.color, lines: keyLines(k) }))
    .reverse();
  const keyBlockH =
    keyGroups.length > 0
      ? KEY_TOP_PAD + keyGroups.length * KEY_ENTRY_GAP +
        keyGroups.reduce(
          (s, g) => s + KEY_HEADER_H + Math.max(0, g.lines.length - 1) * KEY_SUBLINE_H,
          0,
        )
      : 0;

  const m = {
    l: 38,
    r: 8, // no right gutter — traces are named in the key block
    t: 20,
    b: AXIS_FOOT_H + keyBlockH + FOOTNOTE_H,
  };

  // Adaptive width: wide enough for the longest key line and the title, but at
  // least MIN_WIDTH. The key (with its long "a = … · κ = … Å⁻²" detail) is
  // usually the widest element.
  const keyW = keyGroups.reduce(
    (mx, g) => Math.max(mx, ...g.lines.map((ln) => m.l + KEY_INDENT + keyLineWidth(ln))),
    0,
  );
  const titleW = m.l + title.length * FS_TITLE * CHAR_W;
  const width = input.width ?? Math.ceil(Math.max(MIN_WIDTH, keyW + m.r, titleW + m.r));

  // ── Vertical geometry: fixed amplitude + smaller step (overlap), NOT fit ────
  const maxBand = Math.max(1e-6, ...rows.map((r) => Math.max(0, r.bandHeight)));
  const ampOf = (r: WaterfallRow): number =>
    rows.length === 0 ? AMP : (Math.max(0, r.bandHeight) / maxBand) * AMP;
  const step = clamp(STEP_BASE * Math.max(0.1, offsetScale), 18, AMP * 1.6);
  // Trace i (0 = bottom): its band rises `STACK_GAP + i*step` above the axis,
  // up to its amplitude. The plot height is the tallest trace's reach.
  const ph =
    Math.max(
      AMP,
      ...rows.map((r, i) => STACK_GAP + i * step + ampOf(r)),
    ) + TOP_HEAD;

  const height = input.height ?? m.t + ph + m.b;
  const pw = width - m.l - m.r;
  const baseY = m.t + ph;

  const qDomain = input.qDomain ?? waterfallQDomain(rows);
  const x = makeAxis(qDomain, [m.l, m.l + pw], xType);
  const inDomain = (q: number): boolean => q >= qDomain[0] && q <= qDomain[1];
  const ticks = axisTicks(x);

  const bands = rows.map((row, idx) => {
    const bandBottom = baseY - STACK_GAP - idx * step;
    return { row, idx, bandTop: bandBottom - ampOf(row), bandBottom };
  });

  const parts: string[] = [];
  parts.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" ` +
      `viewBox="0 0 ${width} ${height}" font-family="${ARIAL}">`,
  );
  parts.push(`<rect x="0" y="0" width="${width}" height="${height}" fill="#ffffff"/>`);

  // Title (top-left, like a figure caption header).
  parts.push(
    `<text x="${m.l}" y="15" font-size="${FS_TITLE}" font-weight="700" fill="${INK}">${esc(title)}</text>`,
  );

  // L-shaped black axes.
  parts.push(
    `<line x1="${m.l}" y1="${baseY}" x2="${m.l + pw}" y2="${baseY}" stroke="${INK}" stroke-width="1"/>`,
  );
  parts.push(
    `<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${baseY}" stroke="${INK}" stroke-width="1"/>`,
  );

  // X ticks + labels (labels only on non-minor ticks).
  for (const t of ticks) {
    if (!inDomain(t.value)) continue;
    const px = x.to(t.value);
    parts.push(
      `<line x1="${px.toFixed(1)}" y1="${baseY}" x2="${px.toFixed(1)}" y2="${baseY + 3.5}" stroke="${INK}" stroke-width="0.8"/>`,
    );
    if (t.kind !== "minor") {
      parts.push(
        `<text x="${px.toFixed(1)}" y="${baseY + 13}" text-anchor="middle" font-size="${FS_TICK}" fill="${INK}">${t.value.toFixed(2)}</text>`,
      );
    }
  }

  // Axis titles.
  parts.push(
    `<text x="${m.l + pw / 2}" y="${baseY + 25}" text-anchor="middle" font-size="${FS_AXIS}" font-weight="700" fill="${INK}">q (Å⁻¹)</text>`,
  );
  const yMid = m.t + ph / 2;
  parts.push(
    `<text x="11" y="${yMid}" text-anchor="middle" font-size="${FS_AXIS}" font-weight="700" fill="${INK}" transform="rotate(-90 11 ${yMid})">Intensity (a.u.) + offset</text>`,
  );

  // Traces (log-intensity per band) + peak glyphs.
  for (const { row, idx, bandTop, bandBottom } of bands) {
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = colorOf(idx, row);
    const { q, I } = row.trace;

    // Per-row log-intensity mapping: log10(I) ∈ [logMin, logMax] → y ∈
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
        parts.push(`<path d="${d.trim()}" fill="none" stroke="${color}" stroke-width="1.2" data-phase="${row.phase ?? "none"}"/>`);
      }
    }

    // Peak anchors — ordinal-numbered (1..n by ascending q), riding the curve in
    // the trace colour (we don't attribute a peak to a phase under coexistence).
    if (showPeakTicks && hasLog && row.anchors.length > 0) {
      const sorted = [...row.anchors].filter((a) => inDomain(a.q)).sort((a, b) => a.q - b.q);
      sorted.forEach((a, i) => {
        const px = x.to(a.q);
        const cy = yOf(a.intensity ?? 0);
        const tip = cy - 2.5;
        const top = tip - 3.5;
        parts.push(
          `<path d="M${(px - 2.2).toFixed(1)} ${top.toFixed(1)} L${(px + 2.2).toFixed(1)} ${top.toFixed(1)} L${px.toFixed(1)} ${tip.toFixed(1)} Z" fill="${color}"/>`,
        );
        if (showPeakLabels) {
          parts.push(
            `<text x="${px.toFixed(1)}" y="${(top - 1.5).toFixed(1)}" text-anchor="middle" font-size="${FS_ORD}" fill="${color}">${i + 1}</text>`,
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
            `<rect x="${m.l}" y="${(ky - KEY_SWATCH + 0.5).toFixed(1)}" width="${KEY_SWATCH}" height="${KEY_SWATCH}" rx="1" fill="${g.color}"/>`,
          );
        }
        const tx = m.l + KEY_INDENT;
        const fill = ln.muted ? INK_MUTED : INK;
        if (ln.phase) {
          // "<Phase>  <detail>" — phase name bold, lattice/κ detail regular.
          parts.push(
            `<text x="${tx}" y="${ky.toFixed(1)}" font-size="${FS_DETAIL}" fill="${fill}"><tspan font-weight="700">${esc(ln.phase)}</tspan>${ln.text ? `  ${esc(ln.text)}` : ""}</text>`,
          );
        } else {
          parts.push(
            `<text x="${tx}" y="${ky.toFixed(1)}" font-size="${isHeader ? FS_KEY : FS_DETAIL}" font-weight="${isHeader ? 700 : 400}" fill="${fill}">${esc(ln.text)}</text>`,
          );
        }
        ky += isHeader ? KEY_HEADER_H : KEY_SUBLINE_H;
      });
      ky += KEY_ENTRY_GAP;
    }
  }

  // Footnote, at the very bottom.
  parts.push(
    `<text x="${(width / 2).toFixed(1)}" y="${(height - 5).toFixed(1)}" text-anchor="middle" font-size="${FS_FOOT}" fill="${INK_MUTED}">${esc(footer)}</text>`,
  );

  parts.push("</svg>");
  return parts.join("");
}
