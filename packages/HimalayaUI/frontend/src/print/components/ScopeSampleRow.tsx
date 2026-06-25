import { useEffect, useRef, useState } from "react";
import { FlagButton, GripHandle, IconButton, Input } from "../ui";
import { Sparkline } from "../plot/Sparkline";
import { cx } from "../../lib/cx";


export interface ScopeSampleRowProps {
  name: string;
  sampleId: string;
  /** Sparkline data. */
  trace: { q: number[]; I: number[] };
  /** Dominant phase → spark hue; null/undefined → unindexed. */
  phase?: string | null;
  /** The parsed ordering value, e.g. "1 : 0.25". */
  value: string;
  /** Skipped read → FlagButton flagged look + faint accent row wash. */
  flagged?: boolean;
  /** Forwarded to FlagButton onClick. */
  onToggleFlag?: () => void;
  /**
   * Keyboard-reorder contract (SC-KBD, WCAG 2.1.1). When provided, the grip
   * renders as a real button named "Reorder {name}" that handles ArrowUp
   * (`onMoveBy(-1)`) and ArrowDown (`onMoveBy(1)`), preventing default so the
   * page does not scroll. When absent, the grip stays the aria-hidden visual
   * affordance — no dead control. The page owns boundary handling and
   * announcements (it sees the row's position; this row does not).
   */
  onMoveBy?: (delta: -1 | 1) => void;
  /**
   * Inline value correction (SC-VALUECORRECT). When provided, a pencil control
   * appears beside the skip toggle; clicking it swaps the value into a mono
   * Input. Commit (Enter / the check button / blur) calls `onEditValue` with the
   * new value ONLY when it changed AND is non-empty (the write must never carry
   * `value:""` — same guard as the commit gate's `isKept`); Escape cancels. When
   * absent the value stays display + skip only (the cold path and the loose
   * candidate rows have no correctable read). Edit and skip are DISTINCT
   * affordances — correcting a misread is not the same gesture as dropping it.
   */
  onEditValue?: (value: string) => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * One row of the series-scoping worksheet (mockup `.srow`): a drag grip, a small
 * trace sparkline, the sample's name + id, and a trailing parsed-value control.
 *
 * When `flagged` (the user skipped this read from the batch write), the whole
 * row takes a faint accent wash and the value control marks the exclusion. The row
 * carries `group` so the grip brightens on row hover; it also draws its own
 * bottom hairline (the parent strips the last). Pointer drag-reorder is
 * page-deferred (the page wraps rows in `useDragReorder` props); the keyboard
 * reorder path lives HERE via `onMoveBy` — provide it wherever rows are
 * reorderable so the grip is not pointer-only.
 */
export function ScopeSampleRow({
  name,
  sampleId,
  trace,
  phase,
  value,
  flagged,
  onToggleFlag,
  onMoveBy,
  onEditValue,
  className,
}: ScopeSampleRowProps): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);

  // Focus + select on entering edit so the misread value can be retyped or
  // corrected without first clicking into the field.
  useEffect(() => {
    if (editing) {
      inputRef.current?.focus();
      inputRef.current?.select();
    }
  }, [editing]);

  function startEdit(): void {
    setDraft(value);
    setEditing(true);
  }
  // Commit only a CHANGED, non-empty value: an unchanged commit must not push a
  // history entry, and `value:""` would corrupt the sample on the batch write
  // (the commit gate filters it, but never let it leave this control either).
  function commit(): void {
    const next = draft.trim();
    if (next !== "" && next !== value && onEditValue) onEditValue(next);
    setEditing(false);
  }
  function cancel(): void {
    setDraft(value);
    setEditing(false);
  }

  return (
    <div
      data-testid="scope-sample-row"
      data-flagged={flagged ? "true" : "false"}
      className={cx(
        "group flex items-center gap-3 px-2 py-2.5 border-b border-hair",
        flagged && "bg-accent/5",
        className,
      )}
    >
      {onMoveBy ? (
        /* Real, focusable grip: the keyboard path the aria-hidden glyph cannot
           be (FOL-KBD precedent: native button, focus-visible token outline,
           never a dead control). The glyph look is unchanged — the button is a
           transparent wrapper. */
        <button
          type="button"
          aria-keyshortcuts="ArrowUp ArrowDown"
          data-hit-area="24"
          onKeyDown={(e) => {
            // Modified arrows stay native (Cmd+ArrowDown = macOS scroll-to-end,
            // Alt+Arrow = word-nav): hijacking those would be its own trap.
            if (e.metaKey || e.ctrlKey || e.altKey) return;
            if (e.key === "ArrowUp") {
              e.preventDefault();
              onMoveBy(-1);
            } else if (e.key === "ArrowDown") {
              e.preventDefault();
              onMoveBy(1);
            }
          }}
          // Padding grows the hit box to >=24x24 (WCAG 2.5.8) while the matching
          // negative margin cancels it in layout so the glyph and row geometry
          // are unchanged (Checkbox precedent). The glyph is only ~9px wide, so
          // the horizontal padding is px-2 (-> ~25px) while py-1.5 covers height;
          // a symmetric p-1.5 would have left the width at ~21px. inline-flex
          // keeps the padded box centred on the glyph. title surfaces the
          // arrow-key reorder to sighted keyboard users (the shortcut is
          // otherwise only in aria-keyshortcuts).
          title="Drag to reorder, or focus and use ↑ ↓"
          className="inline-flex items-center justify-center flex-shrink-0 py-1.5 -my-1.5 px-2 -mx-2 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
        >
          <span className="sr-only">Reorder {name}</span>
          <GripHandle />
        </button>
      ) : (
        <GripHandle />
      )}
      <Sparkline trace={trace} {...(phase != null ? { phase } : {})} />
      <div className="flex-1 min-w-0">
        <div className="text-body font-semibold text-ink">{name}</div>
        <div className="text-caption font-mono text-ink-soft">{sampleId}</div>
      </div>
      {editing ? (
        /* draggable={false} releases this cluster from the row wrapper's
           draggable=true (useDragReorder) so the value text stays
           click-selectable instead of starting a row drag. */
        <span
          draggable={false}
          className="flex items-center gap-1 flex-shrink-0"
        >
          <Input
            inputRef={inputRef}
            value={draft}
            onValueChange={setDraft}
            mono
            inputSize="sm"
            aria-label={`Corrected value for ${name}`}
            testId="value-edit-input"
            className="w-28"
            onKeyDown={(e: React.KeyboardEvent) => {
              if (e.key === "Enter") {
                e.preventDefault();
                commit();
              } else if (e.key === "Escape") {
                e.preventDefault();
                cancel();
              }
            }}
            onBlur={commit}
          />
          {/* preventDefault on mousedown keeps focus in the field, so the click
              commits without the blur firing a redundant first commit. */}
          <IconButton
            label={`Save value for ${name}`}
            tone="accent"
            onMouseDown={(e) => e.preventDefault()}
            onClick={commit}
          >
            ✓
          </IconButton>
        </span>
      ) : (
        // SC-PENCIL-FAINT: widen the gap so the edit pencil is not crowded
        // against the value's skip toggle (a fat-finger skip when the user meant
        // to edit). The pencil is a quiet secondary affordance: hidden at rest,
        // revealed on row hover OR when it takes keyboard focus, so it reads as
        // a clear "correct this value" control rather than a permanent faint glyph.
        <span className="flex items-center gap-2.5 flex-shrink-0">
          <FlagButton
            value={value}
            {...(flagged ? { flagged: true } : {})}
            {...(onToggleFlag ? { onClick: onToggleFlag } : {})}
          />
          {onEditValue && (
            <IconButton
              label={`Edit value for ${name}`}
              tone="ghost"
              onClick={startEdit}
              className="opacity-0 transition-opacity group-hover:opacity-100 focus-visible:opacity-100"
            >
              ✎
            </IconButton>
          )}
        </span>
      )}
    </div>
  );
}
