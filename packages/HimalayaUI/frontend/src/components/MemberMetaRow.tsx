/**
 * MemberMetaRow — per-trace metadata gutter row (Plan §Phase 7, Task 7.1).
 *
 * Two render paths share the same vertical-positioning math (`top`/`height`
 * from `computeYBands` in `MultiTracePlot`) so reorders + resizes stay in
 * lockstep with the plot.
 *
 * Review-mode (read-only): one line — label / phase chip / d / R² / NGC for
 * cubics. Hover or click expands a detail card (peak count + secondary
 * metadata).
 *
 * Edit-mode: drag-handle (left edge) for reorder + inline controls (label
 * override, normalization, q-window). The numeric q-window inputs follow
 * the focus-gated `QNumInput` pattern from `PlotCard.tsx` so external
 * state changes don't interrupt mid-edit.
 *
 * Stale members surface a present-when-true `data-stale` attribute (matches
 * the codebase convention used by `data-active`, `data-pinned`, etc.) plus
 * an inline warning icon.
 */
import { useCallback, useEffect, useState } from "react";
import type { ComparisonMember } from "../api";
import { useAppState } from "../state";
import { phaseColor, CUBIC_PHASES } from "../phases";
import { COMPARE_PALETTE } from "../lib/comparison/coloring";
import type { DraftMemberNormalization } from "../lib/comparison/draft";

export interface MemberMetaRowProps {
  member: ComparisonMember;
  /** y-band envelope top in pixels (derived from `computeYBands`). */
  top: number;
  /** y-band envelope height in pixels. */
  height: number;
  mode: "review" | "edit";
  /**
   * Index into `activeDraft.members`. Required in edit mode so the row's
   * controls dispatch through `updateMember(index, partial)`.
   */
  memberIndex?: number;
  /**
   * Forwarded by `MemberMetaGutter` so the parent can orchestrate HTML5
   * drag-and-drop reorder across rows. Hooks into the grip's `onDragStart`.
   */
  onGripDragStart?: (e: React.DragEvent) => void;
  /**
   * Pre-resolved display label (typically `${sample.label||sample.name} ·
   * ${exposure.filename}`). When provided, it wins over `defaultLabel`. The
   * Compare review page assembles this from the per-member exposure +
   * sample caches (issue #52). `member.label_override` is honored upstream.
   */
  displayLabel?: string;
}

const NORMALIZATION_OPTIONS: DraftMemberNormalization[] = [
  "none", "max", "area", "qwindow",
];

function formatLattice(d: number | null | undefined): string {
  return d != null ? d.toFixed(2) : "—";
}
function formatR2(r: number | null | undefined): string {
  return r != null ? r.toFixed(3) : "—";
}
function formatKappa(k: number | null | undefined): string {
  return k != null ? k.toFixed(2) : "—";
}

function defaultLabel(member: ComparisonMember): string {
  if (member.label_override) return member.label_override;
  if (member.exposure_id == null) return "(deleted exposure)";
  return `Exposure #${member.exposure_id}`;
}

export function MemberMetaRow(props: MemberMetaRowProps): JSX.Element {
  const { member, top, height, mode, memberIndex, onGripDragStart, displayLabel } = props;
  const [expanded, setExpanded] = useState(false);

  const updateMember = useAppState((s) => s.updateMember);
  const setHighlight = useAppState((s) => s.setHighlightedCompareMemberId);
  const idx = memberIndex ?? -1;

  // Phase 9.5 — only members with a confirmed index participate in
  // hover-driven phase coloring. Hover/focus/click on rows without a
  // confirmed_index are inert (no highlight, not in tab order).
  const canHighlight = member.snapshot?.confirmed_index != null;
  // Pin lifecycle is tracked LOCALLY — the Zustand `highlightedCompareMemberId`
  // is the render target (a single id), but it can't tell us whether the
  // current target was hovered vs. pinned. Tracking pin state on the row
  // lets `onMouseLeave` know whether to clear or keep.
  const [isPinned, setIsPinned] = useState(false);

  const onLabelCommit = useCallback(
    (val: string) => {
      if (idx < 0) return;
      updateMember(idx, { label_override: val === "" ? undefined : val });
    },
    [updateMember, idx],
  );

  const onNormChange = useCallback(
    (val: DraftMemberNormalization) => {
      if (idx < 0) return;
      updateMember(idx, { normalization: val });
    },
    [updateMember, idx],
  );

  const onQWindowMin = useCallback(
    (val: number | undefined) => {
      if (idx < 0) return;
      updateMember(idx, { q_window_min: val });
    },
    [updateMember, idx],
  );
  const onQWindowMax = useCallback(
    (val: number | undefined) => {
      if (idx < 0) return;
      updateMember(idx, { q_window_max: val });
    },
    [updateMember, idx],
  );

  const onResetColor = useCallback(() => {
    if (idx < 0) return;
    updateMember(idx, { color_override: undefined });
  }, [updateMember, idx]);

  const onPickColor = useCallback(
    (color: string) => {
      if (idx < 0) return;
      updateMember(idx, { color_override: color });
    },
    [updateMember, idx],
  );

  // ── hover-driven phase coloring handlers (Phase 9.5) ─────────────────
  // Hover/focus → set highlight; leave/blur → clear UNLESS pinned. Click
  // toggles the pin (sets to this id, or clears if already pinned to this
  // id). Keyboard parity: Enter pins; Esc clears.
  const onMouseEnter = useCallback(() => {
    if (!canHighlight) return;
    setHighlight(member.id);
  }, [canHighlight, member.id, setHighlight]);
  const onMouseLeave = useCallback(() => {
    if (!canHighlight) return;
    if (isPinned) return;
    setHighlight(undefined);
  }, [canHighlight, isPinned, setHighlight]);
  const onFocus = useCallback(() => {
    if (!canHighlight) return;
    setHighlight(member.id);
  }, [canHighlight, member.id, setHighlight]);
  const onBlur = useCallback(() => {
    if (!canHighlight) return;
    if (isPinned) return;
    setHighlight(undefined);
  }, [canHighlight, isPinned, setHighlight]);
  const onKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (!canHighlight) return;
      if (e.key === "Enter") {
        // Pin: track locally + reassert the highlight target.
        setIsPinned(true);
        setHighlight(member.id);
      } else if (e.key === "Escape") {
        setIsPinned(false);
        setHighlight(undefined);
      }
    },
    [canHighlight, member.id, setHighlight],
  );

  const ci = member.snapshot?.confirmed_index ?? null;
  const isCubic = ci !== null && CUBIC_PHASES.has(ci.phase);

  return (
    <div
      data-testid="member-meta-row"
      data-member-id={String(member.id)}
      {...(member.is_stale ? { "data-stale": "" } : {})}
      {...(isPinned ? { "data-highlighted": "" } : {})}
      // Tab into the row only when there's a confirmed index to highlight —
      // otherwise pressing Tab past it is dead air.
      {...(canHighlight ? { tabIndex: 0 } : {})}
      onClick={() => {
        // Click toggles both the expansion (always) and the pin (when
        // hoverable). Pin lifecycle: click pins, click-again unpins. The
        // expanded panel mirrors that — second click closes both.
        setExpanded((e) => !e);
        if (canHighlight) {
          if (isPinned) {
            setIsPinned(false);
            setHighlight(undefined);
          } else {
            setIsPinned(true);
            setHighlight(member.id);
          }
        }
      }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
      onFocus={onFocus}
      onBlur={onBlur}
      onKeyDown={onKeyDown}
      style={{
        position: "absolute",
        left: 0,
        right: 0,
        top: `${top}px`,
        height: `${height}px`,
      }}
      className="flex flex-col gap-0.5 px-2 py-1 text-xs hover:bg-bg-elevated/40 cursor-pointer
                 outline-0 focus-visible:ring-1 focus-visible:ring-accent"
    >
      {/* Primary single line */}
      <div className="flex items-center gap-2 min-w-0">
        {mode === "edit" && (
          <span
            data-testid="member-reorder-grip"
            data-member-id={String(member.id)}
            draggable
            // The actual drag-orchestration lives in the parent gutter
            // (Task 7.3) — this is the visible affordance.
            className="cursor-grab select-none px-1 text-fg-dim hover:text-fg"
            title="Drag to reorder"
            onClick={(e) => e.stopPropagation()}
            {...(onGripDragStart ? { onDragStart: onGripDragStart } : {})}
          >
            ⋮⋮
          </span>
        )}
        <span
          data-testid="member-meta-label"
          className="truncate min-w-0 max-w-[16ch] text-fg"
          title={displayLabel ?? defaultLabel(member)}
        >
          {displayLabel ?? defaultLabel(member)}
        </span>
        {ci !== null ? (
          <>
            <span
              data-testid="member-meta-phase-chip"
              className="text-data-strong px-1.5 py-0.5 rounded-sm border shrink-0"
              style={{
                color: phaseColor(ci.phase),
                background: `color-mix(in oklab, ${phaseColor(ci.phase)} 10%, transparent)`,
                borderColor: `color-mix(in oklab, ${phaseColor(ci.phase)} 35%, transparent)`,
              }}
            >
              {ci.phase}
            </span>
            <span className="text-fg-muted tabular-nums shrink-0">
              <span className="text-fg-dim">d=</span>
              {formatLattice(ci.lattice_d)}
            </span>
            <span className="text-fg-muted tabular-nums shrink-0">
              <span className="text-fg-dim">R²=</span>
              {formatR2(ci.r_squared)}
            </span>
            {isCubic && (
              <span
                data-testid="member-meta-ngc"
                className="text-fg-muted tabular-nums shrink-0"
              >
                <span className="text-fg-dim">κ=</span>
                {formatKappa(ci.ngc)}
              </span>
            )}
          </>
        ) : (
          <span className="text-fg-dim italic">no index</span>
        )}
        {member.is_stale && (
          <span
            data-testid="member-meta-stale-icon"
            title="Stale snapshot — re-submit to refresh"
            className="text-warning shrink-0"
            aria-label="stale"
          >
            ⚠
          </span>
        )}
      </div>

      {/* Edit-mode controls row */}
      {mode === "edit" && (
        <div
          className="flex items-center gap-1 flex-wrap"
          onClick={(e) => e.stopPropagation()}
        >
          <input
            type="text"
            data-testid="member-meta-label-input"
            placeholder="Label override"
            defaultValue={member.label_override ?? ""}
            onBlur={(e) => onLabelCommit(e.currentTarget.value)}
            className="bg-bg border border-border rounded px-1 py-0.5 text-fg text-xs
                       outline-0 focus:border-accent w-[12ch]"
          />
          <select
            data-testid="member-meta-normalization"
            value={member.normalization}
            onChange={(e) =>
              onNormChange(e.currentTarget.value as DraftMemberNormalization)
            }
            className="bg-bg border border-border rounded px-1 py-0.5 text-fg text-xs"
          >
            {NORMALIZATION_OPTIONS.map((n) => (
              <option key={n} value={n}>
                {n}
              </option>
            ))}
          </select>
          <span className="text-fg-dim">q∈[</span>
          <QWindowInput
            value={member.q_window_min ?? null}
            testId="member-meta-qwindow-min"
            onCommit={onQWindowMin}
          />
          <span className="text-fg-dim">,</span>
          <QWindowInput
            value={member.q_window_max ?? null}
            testId="member-meta-qwindow-max"
            onCommit={onQWindowMax}
          />
          <span className="text-fg-dim">]</span>
          {member.color_override != null && (
            <button
              type="button"
              data-testid="member-meta-reset-color"
              onClick={onResetColor}
              className="text-fg-dim hover:text-fg text-xs underline"
            >
              Reset color
            </button>
          )}
        </div>
      )}

      {/* Expand-on-click detail card overlay (review + edit) */}
      {expanded && (
        <div
          data-testid="member-meta-detail"
          className="absolute z-10 left-2 top-full mt-1 bg-bg-elevated border border-border
                     rounded shadow-md p-2 text-xs flex flex-col gap-2 min-w-[180px]"
          onClick={(e) => e.stopPropagation()}
        >
          <div>
            <span className="text-fg-dim">Peaks: </span>
            <span className="tabular-nums">
              {member.snapshot?.effective_peaks?.length ?? 0}
            </span>
          </div>
          {ci !== null && (
            <div>
              <span className="text-fg-dim">Index id: </span>
              <span className="tabular-nums">{ci.id}</span>
            </div>
          )}
          {member.exposure_id !== null && (
            <div>
              <span className="text-fg-dim">Exposure: </span>
              <span className="tabular-nums">#{member.exposure_id}</span>
            </div>
          )}
          {mode === "edit" && (
            <ColorPickerSwatchGrid
              activeColor={member.color_override ?? null}
              onPick={onPickColor}
            />
          )}
        </div>
      )}
    </div>
  );
}

/**
 * 12-color swatch grid for per-member color override (Phase 9.4). Sources
 * its palette from `COMPARE_PALETTE` so the override matches the
 * grouping-mode default palette one-to-one — keeps figures visually
 * consistent within a comparison. Custom hex input is deferred to v2.
 */
interface ColorPickerSwatchGridProps {
  activeColor: string | null;
  onPick: (color: string) => void;
}

function ColorPickerSwatchGrid(
  { activeColor, onPick }: ColorPickerSwatchGridProps,
): JSX.Element {
  return (
    <div
      data-testid="member-color-picker-grid"
      role="radiogroup"
      aria-label="Trace color override"
      className="grid grid-cols-6 gap-1"
    >
      {COMPARE_PALETTE.map((c, i) => {
        const active = c === activeColor;
        return (
          <button
            key={i}
            type="button"
            role="radio"
            aria-checked={active}
            data-testid={`member-color-picker-swatch-${i}`}
            data-color={c}
            data-active={active ? "true" : "false"}
            title={c}
            onClick={(e) => {
              e.stopPropagation();
              onPick(c);
            }}
            style={{ background: c }}
            className={
              "h-5 w-5 rounded-sm border " +
              (active ? "border-fg ring-1 ring-fg" : "border-border hover:border-fg-muted")
            }
          />
        );
      })}
    </div>
  );
}

/**
 * Focus-gated numeric input for q-window (mirrors `QNumInput` from
 * `PlotCard.tsx`). External `value` changes are synced into the draft only
 * when the input is not focused — wheel-zoom or other live updates can't
 * interrupt mid-edit.
 *
 * Empty input commits `undefined` (clears the q-window bound).
 */
interface QWindowInputProps {
  value: number | null;
  testId: string;
  onCommit: (val: number | undefined) => void;
}

function QWindowInput({ value, testId, onCommit }: QWindowInputProps): JSX.Element {
  const initial = value !== null ? value.toFixed(3) : "";
  const [draft, setDraft] = useState(initial);
  const [focused, setFocused] = useState(false);

  useEffect(() => {
    if (!focused) setDraft(value !== null ? value.toFixed(3) : "");
  }, [value, focused]);

  return (
    <input
      type="number"
      step="0.001"
      data-testid={testId}
      value={draft}
      onChange={(e) => setDraft(e.currentTarget.value)}
      onFocus={() => setFocused(true)}
      onBlur={(e) => {
        setFocused(false);
        const raw = e.currentTarget.value.trim();
        if (raw === "") {
          onCommit(undefined);
          return;
        }
        const n = parseFloat(raw);
        if (Number.isFinite(n)) onCommit(n);
      }}
      onKeyDown={(e) => {
        if (e.key === "Enter") {
          (e.currentTarget as HTMLInputElement).blur();
        }
      }}
      className="w-[60px] bg-bg border border-border rounded px-1 py-0.5
                 text-fg text-xs tabular-nums text-right
                 outline-0 focus:border-accent"
    />
  );
}
