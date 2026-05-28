import { useCallback, useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../api";
import { phaseColor } from "../phases";
import { prettifyUnits } from "../lib/units";
import { invertQ } from "../lib/plot/invertQ";
import { formatAxis } from "../lib/plot/formatAxis";

export interface TraceViewerProps {
	trace: Trace;
	peaks: Peak[];
	activeGroupIndices: IndexEntry[];
	hoveredIndex: IndexEntry | undefined;
	/** When set, the matching peak triangle gets a highlight ring (drives the chat-mention peak chip hover effect). */
	hoveredPeakId?: number | undefined;
	/** q-link (#180): the hovered q-value. The peak triangle whose q is within
	 *  tolerance gets the same hover halo as `hoveredPeakId`. Optional — the
	 *  legacy Index path passes neither this nor `onHoverQ`. */
	hoveredQ?: number | undefined;
	/** q-link (#180): fired on cursor move with the (peak-snapped) q, and with
	 *  `null` on leave. Drives the ephemeral `hoveredQ` Zustand field. */
	onHoverQ?: (q: number | null) => void;
	/** R4 candidate-swap preview (L-11): peak ids the hovered candidate would
	 *  orphan if committed (it claims an overlapping peak of an active phase).
	 *  Those triangles dim so the cost of the swap is visible before the click. */
	losingPeakIds?: Set<number> | undefined;
	onAddPeak: (q: number) => void;
	onRemovePeak: (peakId: number) => void;
	onTogglePeakExclusion: (peakId: number, excluded: boolean) => void;
	/** Visible q-range. null = auto (full trace). */
	xDomain: [number, number] | null;
	onXDomain: (d: [number, number] | null) => void;
	/** Visible intensity range. null = auto (full data range). */
	yDomain?: [number, number] | null;
	/** X-axis scale. Defaults to "log" (SAXS convention). */
	xType?: "log" | "linear";
	/** Called on plot dblclick. Defaults to `() => onXDomain(null)`. */
	onReset?: () => void;
	/** Q-axis units label. Defaults to "Å⁻¹" if omitted. */
	qUnits?: string;
}

// ── constants ──────────────────────────────────────────────────────────────

/** Click within this many pixels of a peak triangle to act on it. */
const PEAK_HIT_PX = 10;

/** q-link (#180): relative tolerance for matching `hoveredQ` to a peak's q
 *  when lighting the trace halo (a hovered detector ring carries a peak's
 *  exact q, so this only needs to absorb float noise). */
const Q_LINK_REL_TOL = 0.01;

/** Pixel offset of a peak triangle above the trace. ~half of the previous log-space lift. */
const PEAK_OFFSET_PX = 7;

/** Plot margins (kept in sync with the effect so overlay can clip / position). */
const MARGIN_LEFT = 50;
const MARGIN_RIGHT = 14;
const MARGIN_TOP = 36;
const MARGIN_BOTTOM = 40;

/** Triangle marker geometry. */
const TRIANGLE_HALF_W = 4;
const TRIANGLE_H = 7;

/** When a predicted-q tick has no nearby detected peak, terminate this many pixels above the trace. */
const TICK_END_OFFSET_PX = 14;

/** Pixel tolerance to consider a predicted q "matched" by a detected peak. */
const TICK_MATCH_PX = 12;

/** Predicted-q track row above the plot (carries the phase-colour swatches). */
const TRACK_Y_TOP = 6;
const TRACK_Y_BOTTOM = 18;

// ── helpers ────────────────────────────────────────────────────────────────

export function findNearestPeak(
	peaks: Peak[],
	qClick: number,
	tolerance: number,
): Peak | null {
	let best: Peak | null = null;
	let bestDist = Infinity;
	for (const p of peaks) {
		const d = Math.abs(p.q - qClick);
		if (d < bestDist) {
			best = p;
			bestDist = d;
		}
	}
	return best && bestDist <= tolerance ? best : null;
}

/**
 * Return the click-target peak nearest the cursor pixel, skipping optimistic
 * placeholders (id < 0). Optimistic peaks have no server-side row yet —
 * routing a remove click to one would DELETE /api/peaks/-N, hit a 404,
 * trigger rollback that re-inserts the placeholder, and then the still-pending
 * add would resolve into a real peak the user believed they removed (the
 * "ghost peak" desync). Until the add confirms the click is a no-op; the
 * placeholder renders as an outlined triangle so the unclickability is visible.
 */
export function nearestClickablePeak(
	peaks: Peak[],
	clickX: number,
	toPx: (q: number) => number,
	tolerancePx: number,
): Peak | null {
	let best: Peak | null = null;
	let bestDist = tolerancePx;
	for (const p of peaks) {
		if (p.id < 0) continue;
		const px = toPx(p.q);
		if (!Number.isFinite(px)) continue;
		const d = Math.abs(px - clickX);
		if (d <= bestDist) {
			best = p;
			bestDist = d;
		}
	}
	return best;
}

export function interpolateI(q: number, trace: Trace): number {
	const qs = trace.q;
	let nearest = 0;
	for (let i = 1; i < qs.length; i++) {
		if (Math.abs(qs[i]! - q) < Math.abs(qs[nearest]! - q)) nearest = i;
	}
	return trace.I[nearest]!;
}

type Scale =
	| { invert?: (v: number) => number; apply?: (v: number) => number }
	| undefined;

interface IndexTick {
	q: number;
	color: string;
	indexId: number;
	speculative: boolean;
}
function indexTicks(indices: IndexEntry[]): IndexTick[] {
	return indices.flatMap((ix) =>
		ix.predicted_q.map((q) => ({
			q,
			color: phaseColor(ix.phase),
			indexId: ix.id,
			speculative: ix.kind === "speculative",
		})),
	);
}

// ── component ──────────────────────────────────────────────────────────────

export function TraceViewer({
	trace,
	peaks,
	activeGroupIndices,
	hoveredIndex,
	hoveredPeakId,
	hoveredQ,
	onHoverQ,
	losingPeakIds,
	onAddPeak,
	onRemovePeak,
	onTogglePeakExclusion,
	xDomain,
	onXDomain,
	yDomain = null,
	xType = "log",
	onReset,
	qUnits,
}: TraceViewerProps): JSX.Element {
	const hostRef = useRef<HTMLDivElement>(null);
	const plotContainer = useRef<HTMLDivElement>(null);
	const overlayRef = useRef<SVGSVGElement>(null);
	const plotElRef = useRef<HTMLElement | SVGElement | null>(null);

	const [_resizeKey, setResizeKey] = useState(0);

	useEffect(() => {
		const el = plotContainer.current;
		if (!el) return;
		const obs = new ResizeObserver(() => setResizeKey((k) => k + 1));
		obs.observe(el);
		return () => obs.disconnect();
	}, []);

	// Re-render trigger for the overlay (needed because our deps include peaks
	// and indices; effects close over those values).
	useEffect(() => {
		const host = hostRef.current;
		const container = plotContainer.current;
		if (!host || !container) return;

		const bandData = trace.q.map((q, i) => ({
			q,
			I: trace.I[i]!,
			lo: Math.max(1e-12, trace.I[i]! - trace.sigma[i]!),
			hi: trace.I[i]! + trace.sigma[i]!,
		}));

		const el = Plot.plot({
			width: container.clientWidth || 400,
			height: container.clientHeight || 300,
			marginLeft: MARGIN_LEFT,
			marginRight: MARGIN_RIGHT,
			marginTop: MARGIN_TOP,
			marginBottom: MARGIN_BOTTOM,
			style: {
				fontFamily: "var(--font-sans)",
				color: "var(--color-fg-muted)",
				background: "transparent",
				overflow: "visible",
			},
			x: {
				type: xType,
				label: `q (${prettifyUnits(qUnits ?? "A-1")})`,
				// Plain decimal tick labels — Plot's default SI-suffix formatter
				// renders 0.040 as "40m" which is unhelpful for SAXS q values.
				tickFormat: (d: number) => formatAxis(d),
				...(xDomain ? { domain: xDomain } : {}),
			},
			y: {
				type: "log",
				label: "I (a.u.)",
				tickFormat: (d: number) => formatAxis(d),
				...(yDomain ? { domain: yDomain } : {}),
			},
			// The Plot only renders the data trace + sigma band. Everything else
			// (peak triangles, predicted-q lines, cursor) lives in the overlay so
			// we have full control over geometry, fade-on-hover, and cursor follow.
			marks: [
				Plot.areaY(bandData, {
					x: "q",
					y1: "lo",
					y2: "hi",
					fill: "var(--color-accent)",
					fillOpacity: 0.12,
				}),
				Plot.line(bandData, {
					x: "q",
					y: "I",
					stroke: "var(--color-fg)",
					strokeWidth: 1,
				}),
			],
		});

		container.replaceChildren(el);
		plotElRef.current = el as unknown as HTMLElement;

		// ── click: add / remove / toggle-exclude based on what's near the cursor ─
		function handleClick(ev: Event): void {
			const me = ev as MouseEvent;
			const xScale: Scale = (
				plotElRef.current as unknown as { scale: (n: string) => Scale }
			)?.scale("x");
			if (!xScale?.invert || !xScale.apply) return;
			const rect = container!.getBoundingClientRect();
			const clickX = me.clientX - rect.left;
			const clickY = me.clientY - rect.top;
			if (!insideInterior(clickX, clickY, rect.width, rect.height)) return;

			const bestPeak = nearestClickablePeak(
				peaks,
				clickX,
				xScale.apply!,
				PEAK_HIT_PX,
			);

			if (bestPeak) {
				if (bestPeak.source === "manual") onRemovePeak(bestPeak.id);
				else onTogglePeakExclusion(bestPeak.id, !bestPeak.excluded);
				return;
			}

			// Empty area → add a manual peak at the exact clicked q (no snap).
			const q = xScale.invert(clickX);
			if (Number.isFinite(q) && q > 0) onAddPeak(q);
		}
		(el as unknown as EventTarget).addEventListener("click", handleClick);

		// ── wheel: zoom x-domain around cursor ───────────────────────────────
		function handleWheel(evRaw: Event): void {
			const ev = evRaw as WheelEvent;
			ev.preventDefault();
			const rect = container!.getBoundingClientRect();
			const cursorQ = invertQ(plotElRef.current, ev.clientX - rect.left);
			if (cursorQ === null) return;
			const curMin = xDomain ? xDomain[0] : trace.q[0]!;
			const curMax = xDomain ? xDomain[1] : trace.q[trace.q.length - 1]!;
			const factor = Math.exp(ev.deltaY * 0.001);
			const q0 = trace.q[0]!;
			const qN = trace.q[trace.q.length - 1]!;
			let newMin: number;
			let newMax: number;
			if (xType === "log") {
				const logMin = Math.log(curMin);
				const logMax = Math.log(curMax);
				const logCur = Math.log(Math.max(cursorQ, 1e-6));
				newMin = Math.max(q0, Math.exp(logCur - (logCur - logMin) * factor));
				newMax = Math.min(qN, Math.exp(logCur + (logMax - logCur) * factor));
			} else {
				// Linear x: zoom symmetrically in linear space around the cursor.
				newMin = Math.max(q0, cursorQ - (cursorQ - curMin) * factor);
				newMax = Math.min(qN, cursorQ + (curMax - cursorQ) * factor);
			}
			if (newMax - newMin < (qN - q0) * 1e-4) return;
			onXDomain([newMin, newMax]);
		}
		(el as unknown as EventTarget).addEventListener("wheel", handleWheel, {
			passive: false,
		} as AddEventListenerOptions);

		// ── dblclick: reset domain (delegates to onReset if provided) ───────
		function handleDblClick(): void {
			if (onReset) onReset();
			else onXDomain(null);
		}
		(el as unknown as EventTarget).addEventListener("dblclick", handleDblClick);

		renderOverlay();

		return () => {
			(el as unknown as EventTarget).removeEventListener("click", handleClick);
			(el as unknown as EventTarget).removeEventListener("wheel", handleWheel);
			(el as unknown as EventTarget).removeEventListener(
				"dblclick",
				handleDblClick,
			);
			container.replaceChildren();
			plotElRef.current = null;
		};
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [
		trace,
		peaks,
		activeGroupIndices,
		hoveredIndex,
		xDomain,
		yDomain,
		xType,
		qUnits,
		onAddPeak,
		onRemovePeak,
		onTogglePeakExclusion,
		onXDomain,
		onReset,
		_resizeKey,
	]);

	// ── overlay renderer (peaks + predicted-q lines + cursor) ───────────────
	const renderOverlay = useCallback((): void => {
		const host = hostRef.current;
		const overlay = overlayRef.current;
		if (!host || !overlay) return;
		const plotEl = plotElRef.current as unknown as SVGElement | null;
		if (!plotEl) return;
		const xScale: Scale = (
			plotEl as unknown as { scale: (n: string) => Scale }
		).scale("x");
		const yScale: Scale = (
			plotEl as unknown as { scale: (n: string) => Scale }
		).scale("y");
		if (!xScale?.apply || !yScale?.apply) return;

		const bbox = (plotEl as Element).getBoundingClientRect();
		overlay.setAttribute("width", String(bbox.width));
		overlay.setAttribute("height", String(bbox.height));

		// Wipe and redraw the peak, tick and track layers.
		const peakRoot = overlay.querySelector<SVGGElement>(
			"[data-role=peak-root]",
		)!;
		const tickRoot = overlay.querySelector<SVGGElement>(
			"[data-role=tick-root]",
		)!;
		const trackRoot = overlay.querySelector<SVGGElement>(
			"[data-role=track-root]",
		)!;
		while (peakRoot.firstChild) peakRoot.removeChild(peakRoot.firstChild);
		while (tickRoot.firstChild) tickRoot.removeChild(tickRoot.firstChild);
		while (trackRoot.firstChild) trackRoot.removeChild(trackRoot.firstChild);

		const dimOthers = hoveredIndex !== undefined;
		const insideX = (px: number): boolean =>
			bbox.width === 0 ||
			(px >= MARGIN_LEFT && px <= bbox.width - MARGIN_RIGHT);

		// ── 1. Peak triangles. Compute pixel positions; cache for tick lookup ──
		interface PeakDraw {
			peak: Peak;
			px: number;
			py: number;
		}
		const peakDraws: PeakDraw[] = [];
		for (const p of peaks) {
			const px = xScale.apply!(p.q);
			if (!Number.isFinite(px) || !insideX(px)) continue;
			const I = interpolateI(p.q, trace);
			const py = yScale.apply!(I) - PEAK_OFFSET_PX;
			if (!Number.isFinite(py)) continue;
			peakDraws.push({ peak: p, px, py });
		}

		for (const { peak, px, py } of peakDraws) {
			const isAuto = peak.source === "auto";
			// Bright/neon for "active workflow": auto = ice blue, manual = magenta.
			const baseColor = isAuto
				? "var(--color-accent)"
				: "var(--color-peak-manual)";
			// Auto peaks: filled triangle. Excluded auto peaks: same color but ~30% opacity.
			// Manual peaks: filled magenta triangle (always full opacity when not faded).
			let fill: string;
			let opacity: number;
			const excludedAuto = isAuto && peak.excluded;
			// R4 (L-11): a previewed candidate would orphan this peak on a swap.
			// Dim it (and desaturate) so the cost of committing the swap is
			// visible before the click. Takes precedence over the dimOthers fade.
			const losing = losingPeakIds?.has(peak.id) ?? false;
			if (losing) {
				fill = "var(--color-fg-dim)";
				opacity = 0.3;
			} else if (excludedAuto) {
				// Excluded auto peaks keep their identity (ice blue, ghosted) — they
				// are user curation, not "context to fade away."
				fill = baseColor;
				opacity = 0.3;
			} else if (dimOthers) {
				// Hovering an index: peaks are not the focus → desaturate to gray
				// rather than just thin the alpha. Removing the color signal entirely
				// is what makes the hovered phase pop.
				fill = "var(--color-fg-dim)";
				opacity = 0.5;
			} else {
				fill = baseColor;
				opacity = isAuto ? 0.95 : 1;
			}

			// Hover ring: a faint pulse halo behind the triangle when this peak is
			// hovered elsewhere (chat mention chip). Drawn first so the triangle
			// sits on top.
			// Hover halo: lit by the chat-mention peak chip (`hoveredPeakId`) OR
			// by the q-link (#180) when `hoveredQ` lands on this peak's q. Both
			// channels share the same visual.
			const litByPeakId = hoveredPeakId === peak.id;
			const litByQ = hoveredQ !== undefined
				&& Math.abs(peak.q - hoveredQ) <= peak.q * Q_LINK_REL_TOL;
			if (litByPeakId || litByQ) {
				const halo = document.createElementNS(
					"http://www.w3.org/2000/svg",
					"circle",
				);
				halo.setAttribute("cx", String(px));
				halo.setAttribute("cy", String(py - TRIANGLE_H / 2));
				halo.setAttribute("r", String(TRIANGLE_H + 4));
				halo.setAttribute("fill", "none");
				halo.setAttribute("stroke", baseColor);
				halo.setAttribute("stroke-width", "1.5");
				halo.setAttribute("stroke-opacity", "0.7");
				peakRoot.appendChild(halo);
			}

			const tri = document.createElementNS(
				"http://www.w3.org/2000/svg",
				"polygon",
			);
			tri.setAttribute(
				"points",
				`${px - TRIANGLE_HALF_W},${py - TRIANGLE_H} ` +
					`${px + TRIANGLE_HALF_W},${py - TRIANGLE_H} ` +
					`${px},${py}`,
			);
			// Optimistic placeholders (negative id) render outlined to signal
			// "in flight, not yet interactable". The empty triangle fills in once
			// the server confirms and the row swaps to a positive id.
			const isOptimistic = peak.id < 0;
			if (isOptimistic) {
				tri.setAttribute("fill", "none");
				tri.setAttribute("stroke", baseColor);
				tri.setAttribute("stroke-width", "1.25");
				tri.setAttribute("stroke-opacity", String(opacity));
				tri.setAttribute("data-optimistic", "true");
			} else {
				tri.setAttribute("fill", fill);
				tri.setAttribute("fill-opacity", String(opacity));
				tri.setAttribute("stroke", "var(--color-bg)");
				tri.setAttribute("stroke-width", "0.75");
			}
			tri.setAttribute("data-peak-id", String(peak.id));
			if (losing) tri.setAttribute("data-losing", "true");
			peakRoot.appendChild(tri);
		}

		// Returns true if a detected peak lies within TICK_MATCH_PX of the given q.
		function isPeakMatched(q: number): boolean {
			const px = xScale!.apply!(q);
			return peakDraws.some((d) => Math.abs(d.px - px) <= TICK_MATCH_PX);
		}

		// ── 2. Predicted-q vlines. Adaptive endpoint: terminate above a matched
		//       peak triangle if any, else just above the trace. Matched ticks
		//       render at full opacity; unmatched (predicted but not found) are
		//       noticeably dimmer so the user can see which positions were confirmed.
		function drawTickLine(
			t: IndexTick,
			opts: { strong: boolean; faded: boolean; matched: boolean },
		): void {
			const px = xScale!.apply!(t.q);
			if (!Number.isFinite(px) || !insideX(px)) return;

			// Find a matching peak triangle within TICK_MATCH_PX.
			let matchedPy: number | null = null;
			for (const draw of peakDraws) {
				if (Math.abs(draw.px - px) <= TICK_MATCH_PX) {
					matchedPy = draw.py - TRIANGLE_H - 7;
					break;
				}
			}
			const traceY =
				yScale!.apply!(interpolateI(t.q, trace)) - TICK_END_OFFSET_PX;
			const y2 = matchedPy ?? traceY;

			// Plot vlines stay neutral (gray) until a specific index is hovered —
			// colour now lives in the track row above the plot.
			const stroke = opts.strong ? t.color : "var(--color-fg-dim)";
			const strokeWidth = opts.strong ? "1.5" : "1";
			const strokeOpacity = opts.faded
				? opts.matched
					? "0.3"
					: "0.12"
				: opts.strong
					? opts.matched
						? "1"
						: "0.45"
					: opts.matched
						? "0.35"
						: "0.15";

			const line = document.createElementNS(
				"http://www.w3.org/2000/svg",
				"line",
			);
			line.setAttribute("x1", String(px));
			line.setAttribute("x2", String(px));
			line.setAttribute("y1", String(MARGIN_TOP - 2));
			line.setAttribute("y2", String(y2));
			line.setAttribute("stroke", stroke);
			line.setAttribute("stroke-width", strokeWidth);
			line.setAttribute("stroke-opacity", strokeOpacity);
			line.setAttribute("stroke-linecap", "round");
			// Speculative ticks always dashed — visual warning that the index is
			// hand-built below the auto-analysis minpeaks bar.
			if (t.speculative) line.setAttribute("stroke-dasharray", "2,3");
			tickRoot.appendChild(line);
		}

		// ── 3. Predicted-q track row above the plot. Each tick keeps its phase
		//       colour at all times so the row reads as a quiet legend; default
		//       state matches the old vline opacity, hover lights it solid, and
		//       the off-hover indices fade to gray (same scheme as plot vlines).
		function drawTrackTick(
			t: IndexTick,
			opts: { strong: boolean; faded: boolean; matched: boolean },
		): void {
			const px = xScale!.apply!(t.q);
			if (!Number.isFinite(px) || !insideX(px)) return;

			const stroke = opts.faded ? "var(--color-fg-dim)" : t.color;
			const strokeWidth = opts.strong ? "1.75" : "1.25";
			const strokeOpacity = opts.faded
				? opts.matched
					? "0.3"
					: "0.12"
				: opts.strong
					? opts.matched
						? "1"
						: "0.45"
					: opts.matched
						? "0.55"
						: "0.2";

			const line = document.createElementNS(
				"http://www.w3.org/2000/svg",
				"line",
			);
			line.setAttribute("x1", String(px));
			line.setAttribute("x2", String(px));
			line.setAttribute("y1", String(TRACK_Y_TOP));
			line.setAttribute("y2", String(TRACK_Y_BOTTOM));
			line.setAttribute("stroke", stroke);
			line.setAttribute("stroke-width", strokeWidth);
			line.setAttribute("stroke-opacity", strokeOpacity);
			line.setAttribute("stroke-linecap", "round");
			if (t.speculative) line.setAttribute("stroke-dasharray", "1,2");
			trackRoot.appendChild(line);
		}

		const baseIndices = hoveredIndex
			? activeGroupIndices.filter((ix) => ix.id !== hoveredIndex.id)
			: activeGroupIndices;
		for (const t of indexTicks(baseIndices)) {
			const matched = isPeakMatched(t.q);
			drawTickLine(t, { strong: false, faded: dimOthers, matched });
			drawTrackTick(t, { strong: false, faded: dimOthers, matched });
		}
		if (hoveredIndex) {
			for (const t of indexTicks([hoveredIndex])) {
				const matched = isPeakMatched(t.q);
				drawTickLine(t, { strong: true, faded: false, matched });
				drawTrackTick(t, { strong: true, faded: false, matched });
			}
		}
	}, [peaks, trace, hoveredIndex, hoveredPeakId, hoveredQ, losingPeakIds, activeGroupIndices, xDomain]);

	// Re-render overlay whenever anything that affects it changes.
	useEffect(() => {
		renderOverlay();
	}, [renderOverlay]);

	// ── cursor crosshair (separate, doesn't rebuild the plot) ──────────────
	useEffect(() => {
		const host = hostRef.current;
		const overlay = overlayRef.current;
		if (!host || !overlay) return;

		let rafId = 0;
		function drawCursor(mx: number, my: number): void {
			const plotEl = plotElRef.current as unknown as SVGElement | null;
			if (!plotEl) return;
			const yScale: Scale = (
				plotEl as unknown as { scale: (n: string) => Scale }
			).scale("y");
			if (!yScale?.apply) return;

			const bbox = (plotEl as Element).getBoundingClientRect();
			const relX = mx - bbox.left;
			const relY = my - bbox.top;

			const line = overlay!.querySelector<SVGLineElement>(
				"[data-role=cursor-line]",
			)!;
			const dot = overlay!.querySelector<SVGCircleElement>(
				"[data-role=cursor-dot]",
			)!;
			const label = overlay!.querySelector<SVGTextElement>(
				"[data-role=cursor-label]",
			)!

			if (!insideInterior(relX, relY, bbox.width, bbox.height)) {
				line.setAttribute("opacity", "0");
				dot.setAttribute("opacity", "0");
				label.setAttribute("opacity", "0");
				return;
			}

			const q = invertQ(plotEl, relX);
			if (q === null) return;

			// q-link (#180): emit the hovered q. Peak-snap — if the cursor is
			// within PEAK_HIT_PX of a peak triangle, emit that peak's exact q;
			// otherwise the raw cursor q. nearestClickablePeak takes a PIXEL
			// clickX + a q->px mapper (the x-scale's apply), not an inverted q.
			if (onHoverQ) {
				const xScale: Scale = (
					plotEl as unknown as { scale: (n: string) => Scale }
				).scale("x");
				const snapped = xScale?.apply
					? nearestClickablePeak(peaks, relX, xScale.apply, PEAK_HIT_PX)
					: null;
				onHoverQ(snapped ? snapped.q : q);
			}

			const Iv = interpolateI(q, trace);
			const py = yScale.apply(Iv);

			line.setAttribute("x1", String(relX));
			line.setAttribute("x2", String(relX));
			line.setAttribute("y1", String(MARGIN_TOP));
			line.setAttribute("y2", String(bbox.height - MARGIN_BOTTOM));
			line.setAttribute("opacity", "1");
			dot.setAttribute("cx", String(relX));
			dot.setAttribute("cy", String(py));
			dot.setAttribute("opacity", Number.isFinite(py) ? "1" : "0");

			label.setAttribute("x", String(relX));
			label.setAttribute("y", String(bbox.height - MARGIN_BOTTOM + 14));
			label.textContent = formatAxis(q);
			label.setAttribute("opacity", "1");
		}

		function onMove(ev: MouseEvent): void {
			const mx = ev.clientX,
				my = ev.clientY;
			if (rafId) cancelAnimationFrame(rafId);
			rafId = requestAnimationFrame(() => drawCursor(mx, my));
		}
		function onLeave(): void {
			if (rafId) cancelAnimationFrame(rafId);
			onHoverQ?.(null); // q-link: clear the cross-highlight on leave
			const line = overlay!.querySelector<SVGLineElement>(
				"[data-role=cursor-line]",
			)!;
			const dot = overlay!.querySelector<SVGCircleElement>(
				"[data-role=cursor-dot]",
			)!;
			const label = overlay!.querySelector<SVGTextElement>(
				"[data-role=cursor-label]",
			)!
			line.setAttribute("opacity", "0");
			dot.setAttribute("opacity", "0");
			label.setAttribute("opacity", "0");
		}

		host.addEventListener("mousemove", onMove);
		host.addEventListener("mouseleave", onLeave);
		return () => {
			host.removeEventListener("mousemove", onMove);
			host.removeEventListener("mouseleave", onLeave);
			if (rafId) cancelAnimationFrame(rafId);
		};
	}, [trace, peaks, onHoverQ]);

	return (
		<div
			ref={hostRef}
			className="w-full h-full relative anim-overlay"
			data-testid="trace-viewer"
		>
			<div ref={plotContainer} className="w-full h-full" />
			<svg
				ref={overlayRef}
				className="absolute inset-0 pointer-events-none"
				aria-hidden="true"
			>
				<g data-role="track-root" />
				<g data-role="tick-root" />
				<g data-role="peak-root" />
				<line
					data-role="cursor-line"
					stroke="var(--color-fg-dim)"
					strokeWidth={1}
					strokeOpacity={0.7}
					opacity={0}
				/>
				<circle
					data-role="cursor-dot"
					r={3.5}
					fill="none"
					stroke="var(--color-accent)"
					strokeWidth={1.5}
					opacity={0}
				/>
				<text
					data-role="cursor-label"
					textAnchor="middle"
					dominantBaseline="hanging"
					fill="var(--color-fg-muted)"
					fontSize="10"
					fontFamily="var(--font-mono)"
					opacity={0}
				/>
			</svg>
		</div>
	);
}

function insideInterior(x: number, y: number, w: number, h: number): boolean {
	return (
		x >= MARGIN_LEFT &&
		x <= w - MARGIN_RIGHT &&
		y >= MARGIN_TOP &&
		y <= h - MARGIN_BOTTOM
	);
}
