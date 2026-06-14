// cleanFigureSvg.ts — build the exported series figure as a standalone SVG
// string, using the SAME greenfield d3 axis math (makeAxis/axisTicks) that the
// on-screen WaterfallChart uses. This is the export idiom: a plain, framed
// GraphPad/journal figure (white ground, Arial, literal hex, L-shaped black
// axes) — deliberately NOT "The Print" branding. It renders the SAME data the
// plate composes (traces, indexed-peak anchors, phase colours, offset, q-scale).
import type { WaterfallRow } from "../waterfall/waterfallModel";
import { makeAxis, axisTicks } from "../plot/projection";
import { waterfallQDomain } from "../waterfall/waterfallModel";

const ARIAL = "Arial, Helvetica, sans-serif";
const INK = "#111111";
const INK_MUTED = "#666666";
const TRACE_FALLBACK = "#777777";

/** Phase → literal sRGB hex, converted from phases.ts PHASE_PALETTE (OKLCH) so
 *  the export hue matches the on-screen trace colour. The export idiom is plain
 *  hex (no var(--…) / no OKLCH — the file is a self-contained figure). */
const PHASE_HEX: Record<string, string> = {
  Pn3m: "#b65b00",
  Im3m: "#007d53",
  Ia3d: "#7555a8",
  Fm3m: "#b54952",
  Fd3m: "#855095",
  Hexagonal: "#a44b79",
  Lamellar: "#375fb9",
  Square: "#548126",
};

export interface CleanFigureInput {
  /** Rows low→high (display order); rendered bottom-up, mirroring the plate. */
  rows: WaterfallRow[];
  title: string;
  footer: string;
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

const PER_ROW_H = 96; // default vertical slot per trace (drives default height)
const BOTTOM_GAP = 14;

function phaseHex(phase: string | null): string {
  return (phase != null && PHASE_HEX[phase]) || TRACE_FALLBACK;
}

function esc(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

/** Build the export figure as a standalone `<svg>` markup string. */
export function buildCleanFigureSvg(input: CleanFigureInput): string {
  const {
    rows,
    title,
    footer,
    xType = "log",
    offsetScale = 1,
    showPeakTicks = true,
    showPeakLabels = false,
  } = input;
  const width = input.width ?? 760;
  // Bottom margin stacks four rows below the axis: tick labels, x-axis title,
  // the phase legend, and the footnote — each on its own line so they never
  // collide (the legend used to overprint the x-axis title).
  const m = { l: 56, r: 70, t: 34, b: 96 };
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
  // bandHeight-weighted, the inter-row STEP scaled by the offset, then fit to
  // the panel so the stack always fills it.
  const stackH = ph - BOTTOM_GAP;
  const totalWeight = rows.reduce((s, r) => s + Math.max(0, r.bandHeight), 0) || rows.length;
  let cumulative = 0;
  let natural = 0;
  const nat = [...rows].reverse().map((row) => {
    const h = (Math.max(0, row.bandHeight) / totalWeight) * stackH;
    const top = cumulative + row.yOffset;
    cumulative += h * scale;
    natural = Math.max(natural, top + h);
    return { row, top, h };
  });
  const fit = natural > 0 ? stackH / natural : 1;
  // bandTop/bandBottom in absolute SVG y. Bottom gap lifts row 0 off the axis.
  const bands = nat.map(({ row, top, h }) => {
    const bandTop = m.t + top * fit;
    const bandBottom = bandTop + h * fit;
    return { row, bandTop, bandBottom };
  });

  const parts: string[] = [];
  parts.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" ` +
      `viewBox="0 0 ${width} ${height}" font-family="${ARIAL}">`,
  );
  parts.push(`<rect x="0" y="0" width="${width}" height="${height}" fill="#ffffff"/>`);

  // Title (top-left, like a figure caption header).
  parts.push(
    `<text x="${m.l}" y="22" font-size="15" font-weight="700" fill="${INK}">${esc(title)}</text>`,
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
        `<text x="${px.toFixed(1)}" y="${baseY + 17}" text-anchor="middle" font-size="10" fill="${INK}">${t.value.toFixed(2)}</text>`,
      );
    }
  }

  // Axis titles.
  parts.push(
    `<text x="${m.l + pw / 2}" y="${baseY + 34}" text-anchor="middle" font-size="11.5" font-weight="700" fill="${INK}">q (Å⁻¹)</text>`,
  );
  const yMid = m.t + ph / 2;
  parts.push(
    `<text x="16" y="${yMid}" text-anchor="middle" font-size="11.5" font-weight="700" fill="${INK}" transform="rotate(-90 16 ${yMid})">Intensity (a.u.) + offset</text>`,
  );

  // Traces (log-intensity per band) + peak glyphs + per-row labels.
  for (const { row, bandTop, bandBottom } of bands) {
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = phaseHex(row.phase);
    const { q, I } = row.trace;

    // Per-row log-intensity mapping (matches multiTraceExportMarks): map
    // log10(I) ∈ [logMin, logMax] → y ∈ [bandBottom, bandTop] (low at bottom).
    let minI = Infinity;
    let maxI = -Infinity;
    for (const v of I) {
      if (Number.isFinite(v) && v > 0) {
        if (v < minI) minI = v;
        if (v > maxI) maxI = v;
      }
    }
    if (maxI > 0) {
      const logMin = Math.log10(minI);
      const logRange = Math.max(1e-6, Math.log10(maxI) - logMin);
      let d = "";
      let started = false;
      for (let k = 0; k < q.length; k++) {
        const qk = q[k]!;
        if (!inDomain(qk)) continue;
        const v = I[k] ?? 0;
        const li = v > 0 ? Math.log10(v) : logMin;
        const py = bandBottom - ((li - logMin) / logRange) * bandH;
        d += `${started ? "L" : "M"}${x.to(qk).toFixed(1)} ${py.toFixed(1)} `;
        started = true;
      }
      if (d) {
        parts.push(`<path d="${d.trim()}" fill="none" stroke="${color}" stroke-width="2" data-phase="${row.phase ?? "none"}"/>`);
      }
    }

    // Peak anchors — ordinal-numbered (1..n by ascending q), a small filled
    // triangle just above the band top; ordinal label above when enabled.
    if (showPeakTicks && row.anchors.length > 0) {
      const sorted = [...row.anchors].filter((a) => inDomain(a.q)).sort((a, b) => a.q - b.q);
      sorted.forEach((a, i) => {
        const px = x.to(a.q);
        const ty = bandTop - 2;
        parts.push(
          `<path d="M${(px - 3).toFixed(1)} ${ty.toFixed(1)} L${(px + 3).toFixed(1)} ${ty.toFixed(1)} L${px.toFixed(1)} ${(ty - 5).toFixed(1)} Z" fill="${color}"/>`,
        );
        if (showPeakLabels) {
          parts.push(
            `<text x="${px.toFixed(1)}" y="${(ty - 7).toFixed(1)}" text-anchor="middle" font-size="9" fill="${color}">${i + 1}</text>`,
          );
        }
      });
    }

    // Per-row label in the right gutter, at the band's vertical centre.
    parts.push(
      `<text x="${m.l + pw + 6}" y="${(bandTop + bandH * 0.5).toFixed(1)}" font-size="10" fill="${INK}" dominant-baseline="middle">${esc(row.label)}</text>`,
    );
  }

  // Phase legend (distinct phases present), centred under the axis title.
  const legendPhases: string[] = [];
  for (const r of rows) {
    const p = r.phase ?? "unphased";
    if (!legendPhases.includes(p)) legendPhases.push(p);
  }
  if (legendPhases.length > 0) {
    const legendY = baseY + 58;
    // Size each item from the RENDERED label (not the phase key — "unphased"
    // displays as the longer "unphased / unbound"), so items can't overlap.
    const GAP = 18;
    const items = legendPhases.map((p) => {
      const label = p === "unphased" ? "unphased / unbound" : p;
      const hex = p === "unphased" ? TRACE_FALLBACK : phaseHex(p);
      return { label, hex, w: 16 + label.length * 6.1 };
    });
    const totalW = items.reduce((s, it) => s + it.w + GAP, 0) - GAP;
    let cx = Math.max(m.l, (width - totalW) / 2);
    for (const it of items) {
      parts.push(`<rect x="${cx.toFixed(1)}" y="${legendY - 8}" width="10" height="10" fill="${it.hex}"/>`);
      parts.push(`<text x="${(cx + 15).toFixed(1)}" y="${legendY}" font-size="10" fill="${INK}">${esc(it.label)}</text>`);
      cx += it.w + GAP;
    }
  }

  // Footnote, centred under the legend.
  parts.push(
    `<text x="${(width / 2).toFixed(1)}" y="${baseY + 80}" text-anchor="middle" font-size="10" fill="${INK_MUTED}">${esc(footer)}</text>`,
  );

  parts.push("</svg>");
  return parts.join("");
}
