/**
 * MemberMetaGutter — composes the per-trace metadata rows + (in edit mode)
 * the band-resize dividers and HTML5 drag-and-drop reorder (Plan §Phase 7,
 * Tasks 7.1, 7.2, 7.3).
 *
 * Vertical positioning sources from the SAME `computeYBands` math the plot
 * uses, so reorder + resize stay aligned without any duplicate state.
 *
 * Reorder uses the platform's HTML5 drag-and-drop API — `@dnd-kit` was
 * not in the codebase as of Phase 7 (no existing reuse path despite the
 * plan's reference to ThumbnailGallery, which doesn't have reorder). The
 * grip on `MemberMetaRow` carries `draggable`; this component listens at
 * the row container level for `dragover` / `drop`, decodes the source
 * member id from `dataTransfer`, and dispatches `reorderMembers(newOrder)`.
 *
 * Test selector contract:
 *   - `data-testid="member-meta-gutter"` for the container
 *   - rows / dividers expose their own testids (see MemberMetaRow,
 *     BandResizeDivider).
 */
import { useCallback, useEffect, useRef, useState } from "react";
import type { ComparisonMember } from "../api";
import { useAppState } from "../state";
import { computeYBands } from "../lib/comparison/yBands";
import { useActiveBand } from "./ActiveBandContext";
import { MemberMetaRow } from "./MemberMetaRow";
import { BandResizeGap } from "./BandResizeDivider";

export interface MemberMetaGutterProps {
  /** Members in render order (already sorted by display_order). */
  members: ComparisonMember[];
  /** Panel height in pixels — drives the y-band envelopes. */
  panelHeight: number;
  mode: "review" | "edit";
  /**
   * Per-member display label resolved by the parent via
   * `resolveDisplayLabels` (`lib/comparison/labels.ts`) — the single
   * source of truth for the fallback chain (issues #52, #69, #73). The
   * resolver always populates an entry for every member id, so this
   * component does no fallback reconstruction of its own; doing so would
   * re-introduce the drift seam #73 closed.
   */
  displayLabelByMemberId: Map<number, string>;
  /**
   * Compare UX E-4 — optional observer of an inter-row resize drag.
   * Receives `{ dy }` (cumulative pointer delta) on every pointerMove and
   * a final call on pointerUp carrying the total dy. The band-height
   * mutation itself goes through the existing `resizeBands` Zustand action
   * inside `BandResizeGap`; this callback is purely for parents that want
   * to observe the gesture. Genuinely optional under
   * `exactOptionalPropertyTypes` — most call sites omit it.
   */
  onResize?: (delta: { dy: number }) => void;
}

const DRAG_MIME = "application/x-himalaya-member-id";

export function MemberMetaGutter(props: MemberMetaGutterProps): JSX.Element {
  const { members, panelHeight, mode, displayLabelByMemberId, onResize } = props;
  const reorderMembers = useAppState((s) => s.reorderMembers);
  const { setDropTargetMemberId } = useActiveBand();
  const dragSourceIdxRef = useRef<number | null>(null);
  // Compare UX E-5 fix — reactive mirror of `dragSourceIdxRef`. The ref is
  // kept because the native HTML5 `drop` handler must read the source index
  // SYNCHRONOUSLY at event time (a prior Critical fix made `drop` the single
  // reorder-commit path — see `handleDrop`). But refs are not reactive, so
  // the render-time `dropIndicatorY` edge decision reads THIS state instead.
  // Ref and state are set/cleared together at every site so they never
  // diverge: primed in `handleDragStart` / `handleRowDragStart`, cleared in
  // `endReorderRef` (which `handleDrop` / `handleDragEnd` call).
  const [dragSourceIdx, setDragSourceIdx] = useState<number | null>(null);
  const [hoverIdx, setHoverIdx] = useState<number | null>(null);
  // Compare UX E-5 — reorder drop indicator. While a pointer-drag reorder
  // is active, `dropTargetIdx` is the band the dragged row would land in
  // (also the band mirrored on the plot via `ActiveBandContext`). `null`
  // means no reorder drag in progress, so no indicator renders.
  const [dropTargetIdx, setDropTargetIdx] = useState<number | null>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  // Compare UX E-2 — the single-expanded-at-a-time invariant lives here.
  // Exactly one member id may be expanded; `null` means all collapsed.
  // Threaded into each `MemberMetaRow` as a controlled `expanded` prop.
  const [expandedMemberId, setExpandedMemberId] = useState<number | null>(null);

  const ratios = members.map((m) => m.band_height || 1);
  const yBands = computeYBands(ratios, panelHeight);

  const handleDragStart = useCallback(
    (e: React.DragEvent, fromIdx: number) => {
      dragSourceIdxRef.current = fromIdx;
      setDragSourceIdx(fromIdx);
      // Some browsers reject the drag without setData — empty payload + a
      // custom MIME satisfies the API while keeping the source-of-truth in
      // the ref (faster, no string parsing).
      try {
        e.dataTransfer.setData(DRAG_MIME, String(members[fromIdx]?.id ?? ""));
        e.dataTransfer.effectAllowed = "move";
      } catch {
        // Older JSDOM versions throw on setData — we still drive the
        // reorder via the ref, so this is non-fatal.
      }
    },
    [members],
  );

  // Compare UX E-3 — grab-anywhere reorder. `MemberMetaRow` calls this the
  // moment a pointer gesture crosses the 4px drag threshold (off the grip
  // OR anywhere on the row body). It primes the same `dragSourceIdxRef`
  // the grip's native `dragstart` sets, so the existing `dragover`/`drop`
  // listeners below complete the reorder unchanged — no second drag path.
  const handleRowDragStart = useCallback(
    (memberId: number) => {
      const fromIdx = members.findIndex((m) => m.id === memberId);
      if (fromIdx >= 0) {
        dragSourceIdxRef.current = fromIdx;
        setDragSourceIdx(fromIdx);
      }
    },
    [members],
  );

  // ── Compare UX E-5 — drop indicator + plot mirror ────────────────────
  // While a pointer-drag reorder is active (`dragSourceIdxRef` primed by
  // E-3's threshold gate), each pointer-move maps the pointer's y onto the
  // `computeYBands` envelopes to find the band the dragged row is over.
  // That band index drives the inline-positioned 2px drop indicator and is
  // published to `ActiveBandContext` so the plot's matching per-band
  // overlay mirrors it with `data-drop-target`.
  //
  // The row body holds `setPointerCapture` (E-3), so its move/up/cancel
  // events still bubble to this container — these handlers ride that same
  // bubbling, no second drag path. Outside a reorder drag they no-op.
  //
  // `publishedDropRef` mirrors what was last pushed to the context so the
  // unmount-cleanup effect can clear a leaked drop-target without stomping
  // unrelated state (the exact class of bug an E-4 review caught).
  const publishedDropRef = useRef<number | null>(null);
  const publishDropTarget = useCallback(
    (idx: number | null) => {
      const memberId = idx !== null ? members[idx]?.id ?? null : null;
      publishedDropRef.current = memberId;
      setDropTargetMemberId(memberId);
    },
    [members, setDropTargetMemberId],
  );

  // Compare UX E-5 fix — visual-only teardown. Clears the drop indicator
  // (`dropTargetIdx`/`hoverIdx`) and the plot-mirror context, but does NOT
  // touch `dragSourceIdxRef`. The native HTML5 `drop` is the SINGLE
  // reorder-commit path and must still be able to read `dragSourceIdxRef`
  // after `pointercancel` fires (a real browser fires `pointercancel` the
  // instant a `draggable` element starts a native drag, BEFORE `drop`).
  // Every gesture-end caller runs this; `dragSourceIdxRef` teardown is
  // owned exclusively by `drop`/`dragend` (see `endReorderRef`).
  const endReorderVisuals = useCallback(() => {
    setDropTargetIdx(null);
    setHoverIdx(null);
    publishDropTarget(null);
  }, [publishDropTarget]);

  // Compare UX E-5 fix — `dragSourceIdxRef` teardown. Owned by the native
  // drag lifecycle (`drop` / `dragend`) only, so the ref survives the
  // `pointercancel` that the browser fires at native-dragstart and the
  // commit in `handleDrop` can still read it.
  const endReorderRef = useCallback(() => {
    dragSourceIdxRef.current = null;
    setDragSourceIdx(null);
  }, []);

  const handleReorderPointerMove = useCallback(
    (e: React.PointerEvent) => {
      if (dragSourceIdxRef.current === null) return;
      const container = containerRef.current;
      if (!container) return;
      const localY = e.clientY - container.getBoundingClientRect().top;
      // Clamp into [0, panelHeight] then locate the containing band.
      const y = Math.max(0, Math.min(panelHeight, localY));
      let idx = yBands.findIndex(([top, bottom]) => y >= top && y < bottom);
      if (idx < 0) idx = members.length - 1; // pointer past the last band
      setDropTargetIdx(idx);
      publishDropTarget(idx);
    },
    [panelHeight, yBands, members.length, publishDropTarget],
  );

  // Unmount cleanup — if the gutter is torn down mid-drag, pointerup/cancel
  // never fire, so the published drop-target would leak into the context.
  // Clear it on unmount iff this gutter is the live publisher.
  useEffect(() => {
    return () => {
      if (publishedDropRef.current !== null) setDropTargetMemberId(null);
    };
  }, [setDropTargetMemberId]);

  const handleDragOver = useCallback((e: React.DragEvent, overIdx: number) => {
    if (dragSourceIdxRef.current === null) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = "move";
    setHoverIdx(overIdx);
  }, []);

  // Compare UX E-5 fix — the native HTML5 `drop` is the SINGLE
  // reorder-commit path. It reads `dragSourceIdxRef` (primed at
  // `dragstart`), commits the reorder, then tears down BOTH the visuals
  // and the ref. `pointercancel` (fired by the browser at native
  // dragstart) deliberately leaves `dragSourceIdxRef` intact so this
  // handler can still see the source index.
  const handleDrop = useCallback(
    (e: React.DragEvent, toIdx: number) => {
      e.preventDefault();
      const fromIdx = dragSourceIdxRef.current;
      endReorderRef();
      endReorderVisuals();
      if (fromIdx === null) return;
      if (fromIdx === toIdx) return;
      // Build newOrder by lifting fromIdx out and inserting at toIdx.
      const order: number[] = [];
      for (let i = 0; i < members.length; i++) order.push(i);
      order.splice(fromIdx, 1);
      order.splice(toIdx, 0, fromIdx);
      reorderMembers(order);
    },
    [members.length, reorderMembers, endReorderRef, endReorderVisuals],
  );

  // Compare UX E-5 fix — `dragend` is the native-drag completion event;
  // it fires whether or not a `drop` landed (e.g. drag aborted off-target,
  // or pressing Esc mid-drag). It commits nothing — `drop` owns the
  // commit — but it MUST clear `dragSourceIdxRef` so a cancelled drag
  // doesn't leak a stale source index into the next gesture. Visuals are
  // cleared too in case the drag ended without a `drop`.
  const handleDragEnd = useCallback(() => {
    endReorderRef();
    endReorderVisuals();
  }, [endReorderRef, endReorderVisuals]);

  // Compare UX E-5 fix — `pointercancel` fires the instant a `draggable`
  // element begins a native drag (the real-browser reorder sequence:
  // pointerdown → pointermove → dragstart → pointercancel → drop). It is
  // VISUAL-ONLY teardown: it must NOT null `dragSourceIdxRef`, or the
  // native `drop` that follows would read `null` and the reorder would
  // never commit. The native `drop`/`dragend` own the ref's teardown.
  const handleReorderPointerCancel = useCallback(() => {
    endReorderVisuals();
  }, [endReorderVisuals]);

  // Compare UX E-5 — y of the drop indicator. The dragged row drops INTO
  // `dropTargetIdx`; draw the 2px rule at that band's near edge — the top
  // edge when dropping above the source, the bottom edge when below — so
  // the indicator reads as "the row lands here". The edge decision reads
  // the reactive `dragSourceIdx` STATE (not `dragSourceIdxRef`) so this
  // render-path branch genuinely re-evaluates whenever the source index
  // changes — no reliance on an incidental co-occurring render.
  const dropBand = dropTargetIdx !== null ? yBands[dropTargetIdx] : undefined;
  const dropIndicatorY =
    dropBand && dragSourceIdx !== null
      ? dropTargetIdx! <= dragSourceIdx
        ? dropBand[0]
        : dropBand[1]
      : null;

  return (
    <div
      ref={containerRef}
      data-testid="member-meta-gutter"
      className="relative h-full w-full"
      style={{ height: `${panelHeight}px` }}
      {...(mode === "edit"
        ? {
            // Compare UX E-5 fix — `onPointerMove` tracks the indicator;
            // `onPointerCancel` does visual-only teardown; `onDragEnd`
            // (native) clears `dragSourceIdxRef`. There is intentionally
            // NO `onPointerUp` reorder handler: the native HTML5 `drop`
            // is the SINGLE reorder-commit path. A row-body native drag
            // consumes the pointer, so no `pointerup` is delivered for a
            // reorder gesture anyway.
            onPointerMove: handleReorderPointerMove,
            onPointerCancel: handleReorderPointerCancel,
            onDragEnd: handleDragEnd,
          }
        : {})}
    >
      {members.map((m, i) => {
        const band = yBands[i] ?? [0, 0];
        const top = band[0];
        const height = band[1] - band[0];
        const isHover = hoverIdx === i && dragSourceIdxRef.current !== null;
        return (
          <div
            key={m.id}
            // The drop target wraps the row so dragover/drop fire across the
            // entire band envelope.
            onDragOver={mode === "edit" ? (e) => handleDragOver(e, i) : undefined}
            onDrop={mode === "edit" ? (e) => handleDrop(e, i) : undefined}
            style={{
              position: "absolute",
              left: 0,
              right: 0,
              top: `${top}px`,
              height: `${height}px`,
            }}
            className={isHover ? "ring-1 ring-accent/50 rounded" : undefined}
          >
            <MemberMetaRow
              member={m}
              top={0}
              height={height}
              mode={mode}
              memberIndex={i}
              // `resolveDisplayLabels` always populates an entry per member;
              // the `?? ""` is a TypeScript-narrowing safety net, never a
              // user-facing fallback (those all live in `labels.ts`, #73).
              displayLabel={displayLabelByMemberId.get(m.id) ?? ""}
              // Compare UX E-2 — controlled collapse/expand. Toggling a row
              // expands it (clearing any other), or collapses it if it was
              // the expanded one.
              expanded={expandedMemberId === m.id}
              onToggleExpand={() =>
                setExpandedMemberId((cur) => (cur === m.id ? null : m.id))
              }
              {...(mode === "edit"
                ? {
                    onGripDragStart: (e: React.DragEvent) => handleDragStart(e, i),
                    onDragStart: handleRowDragStart,
                  }
                : {})}
            />
          </div>
        );
      })}
      {mode === "edit" &&
        members.slice(0, -1).map((m, i) => {
          const band = yBands[i];
          if (!band) return null;
          const boundaryY = band[1]; // bottom of band[i] = top of band[i+1]
          const next = members[i + 1]!;
          return (
            <BandResizeGap
              key={`gap-${m.id}-${next.id}`}
              aboveId={m.id}
              belowId={next.id}
              memberIndex={i}
              top={boundaryY}
              totalHeightPx={panelHeight}
              {...(onResize ? { onResize } : {})}
            />
          );
        })}
      {/* Compare UX E-5 — 2px reorder drop indicator. Rendered only while a
          pointer-drag reorder is active; height + y come from inline style
          so JSDOM tests can assert them without a layout pass. */}
      {dropIndicatorY !== null && (
        <div
          data-testid="drop-indicator"
          aria-hidden="true"
          className="absolute left-0 right-0 z-20 bg-accent rounded-full pointer-events-none"
          style={{
            top: `${dropIndicatorY - 1}px`,
            height: "2px",
          }}
        />
      )}
    </div>
  );
}
