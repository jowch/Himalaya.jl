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
import { useCallback, useRef, useState } from "react";
import type { ComparisonMember } from "../api";
import { useAppState } from "../state";
import { computeYBands } from "./MultiTracePlot";
import { MemberMetaRow } from "./MemberMetaRow";
import { BandResizeDivider } from "./BandResizeDivider";

export interface MemberMetaGutterProps {
  /** Members in render order (already sorted by display_order). */
  members: ComparisonMember[];
  /** Panel height in pixels — drives the y-band envelopes. */
  panelHeight: number;
  mode: "review" | "edit";
  /**
   * Optional per-member display label resolved by the parent (typically
   * `${sample.label||sample.name} · ${exposure.filename}`, with
   * `member.label_override` honored first). When absent, MemberMetaRow
   * falls back to its internal default. Issue #52.
   */
  displayLabelByMemberId?: Map<number, string>;
}

const DRAG_MIME = "application/x-himalaya-member-id";

export function MemberMetaGutter(props: MemberMetaGutterProps): JSX.Element {
  const { members, panelHeight, mode, displayLabelByMemberId } = props;
  const reorderMembers = useAppState((s) => s.reorderMembers);
  const dragSourceIdxRef = useRef<number | null>(null);
  const [hoverIdx, setHoverIdx] = useState<number | null>(null);

  const ratios = members.map((m) => m.band_height || 1);
  const yBands = computeYBands(ratios, panelHeight);

  const handleDragStart = useCallback(
    (e: React.DragEvent, fromIdx: number) => {
      dragSourceIdxRef.current = fromIdx;
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

  const handleDragOver = useCallback((e: React.DragEvent, overIdx: number) => {
    if (dragSourceIdxRef.current === null) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = "move";
    setHoverIdx(overIdx);
  }, []);

  const handleDrop = useCallback(
    (e: React.DragEvent, toIdx: number) => {
      e.preventDefault();
      const fromIdx = dragSourceIdxRef.current;
      dragSourceIdxRef.current = null;
      setHoverIdx(null);
      if (fromIdx === null) return;
      if (fromIdx === toIdx) return;
      // Build newOrder by lifting fromIdx out and inserting at toIdx.
      const order: number[] = [];
      for (let i = 0; i < members.length; i++) order.push(i);
      order.splice(fromIdx, 1);
      order.splice(toIdx, 0, fromIdx);
      reorderMembers(order);
    },
    [members.length, reorderMembers],
  );

  return (
    <div
      data-testid="member-meta-gutter"
      className="relative h-full w-full"
      style={{ height: `${panelHeight}px` }}
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
              {...(displayLabelByMemberId?.get(m.id) !== undefined
                ? { displayLabel: displayLabelByMemberId.get(m.id)! }
                : {})}
              {...(mode === "edit"
                ? { onGripDragStart: (e: React.DragEvent) => handleDragStart(e, i) }
                : {})}
            />
          </div>
        );
      })}
      {mode === "edit" &&
        members.slice(0, -1).map((m, i) => {
          const band = yBands[i];
          if (!band) return null;
          const dividerY = band[1]; // bottom of band[i] = top of band[i+1]
          const next = members[i + 1]!;
          return (
            <BandResizeDivider
              key={`divider-${m.id}-${next.id}`}
              aboveId={m.id}
              belowId={next.id}
              memberIndex={i}
              top={dividerY}
              totalHeightPx={panelHeight}
            />
          );
        })}
    </div>
  );
}
