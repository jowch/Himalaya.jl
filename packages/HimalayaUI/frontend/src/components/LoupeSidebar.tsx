import { useCallback, useState } from "react";
import type { CorpusSample, Exposure } from "../api";
import { useAddCorpusSampleTag, useRemoveCorpusSampleTag } from "../queries";

interface Props {
  /** The active exposure. */
  exposure: Exposure;
  /** All exposures for the sample — drives the "frame N of M" position. */
  exposures: Exposure[];
  /** The sample, for the tags section + corpus-tag mutations (#159). */
  sample: CorpusSample;
  /**
   * Signal-strength level, 0–5, for the meta-list meter (M-8). Derived by the
   * parent from the active exposure's detected-peak count — an honest proxy for
   * ordering-signal strength (a strongly-ordered phase resolves more Bragg
   * reflections). Clamped to 0..5 here.
   */
  signalLevel: number;
  /** Drop the exposure if kept, restore it if dropped. */
  onDropToggle: () => void;
  /** Mark the active exposure as the indexing representative. */
  onSetRepresentative: () => void;
}

function MetaRow(
  { label, value, testid }: { label: string; value: string; testid: string },
): JSX.Element {
  return (
    <div className="flex justify-between font-mono text-[11.5px]">
      <span className="text-ink-faint">{label}</span>
      <span data-testid={testid} className="text-ink">{value}</span>
    </div>
  );
}

/**
 * SignalMeter — five fixed bars, the leftmost `level` filled (mockup
 * `.signal-bars`). `level` is clamped to 0..5. Used as a meta-row value.
 */
function SignalMeter({ level }: { level: number }): JSX.Element {
  const on = Math.max(0, Math.min(5, Math.round(level)));
  return (
    <span data-testid="loupe-meta-signal" className="inline-flex gap-0.5">
      {[0, 1, 2, 3, 4].map((i) => (
        <i
          key={i}
          data-bar={i < on ? "on" : "off"}
          className={[
            "h-[11px] w-[5px] rounded-[1px]",
            i < on ? "bg-ink-soft" : "bg-hair-strong",
          ].join(" ")}
        />
      ))}
    </span>
  );
}

function SectionHeading({ children }: { children: string }): JSX.Element {
  return (
    <div className="mb-2 text-[10px] font-bold uppercase tracking-wide text-ink-faint">
      {children}
    </div>
  );
}

/**
 * LoupeTagsEditor — the read-write tag chip rail for the loupe sidebar
 * (#159 / #207). Each chip carries an inline remove ✕ (aria-label
 * "Remove <key> tag"), and the trailing `+ tag` button reveals a tiny
 * inline form (placeholder "key" + "value", Add). Mutations route through
 * `useAddCorpusSampleTag` / `useRemoveCorpusSampleTag` — the same queue
 * mutators the contact sheet wires to. The corpus cache row is the single
 * source of truth, so an add here is visible on the next contact-sheet
 * paint without a refetch round-trip.
 */
function LoupeTagsEditor({ sample }: { sample: CorpusSample }): JSX.Element {
  const add = useAddCorpusSampleTag(sample.id);
  const remove = useRemoveCorpusSampleTag(sample.id);
  const [editing, setEditing] = useState(false);
  const [keyDraft, setKeyDraft] = useState("");
  const [valDraft, setValDraft] = useState("");

  const reset = useCallback(() => {
    setKeyDraft("");
    setValDraft("");
    setEditing(false);
  }, []);

  const submit = useCallback(() => {
    const v = valDraft.trim();
    // The value carries the user-visible token; key is optional metadata
    // (`""` when absent, mirroring SampleMetadataCard's contract).
    if (v === "") return;
    add.mutate({ key: keyDraft.trim(), value: v });
    reset();
  }, [add, keyDraft, valDraft, reset]);

  return (
    <div data-testid="loupe-tags" className="flex flex-wrap items-center gap-1.5">
      {sample.tags.length === 0 && !editing ? (
        <button
          type="button"
          data-testid="loupe-tag-add"
          onClick={() => setEditing(true)}
          className="rounded-full border border-dashed border-hair-strong px-2 py-0.5
                     text-[10.5px] font-semibold text-ink-faint
                     hover:border-print-accent hover:text-print-accent"
        >
          + tag
        </button>
      ) : (
        <>
          {sample.tags.map((tag) => (
            <span
              key={tag.id}
              className="inline-flex items-center gap-1 rounded-full border border-hair
                         bg-plate px-2 py-0.5 text-[10.5px] font-semibold text-ink-soft"
            >
              {tag.value}
              <button
                type="button"
                aria-label={`Remove ${tag.key || tag.value} tag`}
                onClick={() => remove.mutate(tag.id)}
                className="text-ink-faint hover:text-print-accent"
              >
                ×
              </button>
            </span>
          ))}
          {editing ? (
            <span
              data-testid="loupe-tag-form"
              className="inline-flex items-center gap-1 rounded-full border border-hair-strong
                         bg-plate px-1.5 py-0.5"
            >
              <input
                aria-label="tag key"
                placeholder="key"
                value={keyDraft}
                onChange={(e) => setKeyDraft(e.target.value)}
                className="w-12 bg-transparent text-[10.5px] text-ink outline-none
                           placeholder:text-ink-faint"
              />
              <span className="text-ink-faint">:</span>
              <input
                aria-label="tag value"
                placeholder="value"
                value={valDraft}
                onChange={(e) => setValDraft(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") submit();
                  else if (e.key === "Escape") reset();
                }}
                autoFocus
                className="w-16 bg-transparent text-[10.5px] text-ink outline-none
                           placeholder:text-ink-faint"
              />
              <button
                type="button"
                onClick={submit}
                className="text-[10.5px] font-semibold text-print-accent
                           hover:underline"
              >
                Add
              </button>
              <button
                type="button"
                aria-label="cancel adding tag"
                onClick={reset}
                className="text-ink-faint hover:text-print-accent"
              >
                ×
              </button>
            </span>
          ) : (
            <button
              type="button"
              data-testid="loupe-tag-add"
              onClick={() => setEditing(true)}
              title="add a tag"
              className="rounded-full border border-dashed border-hair-strong px-2 py-0.5
                         text-[10.5px] font-semibold text-ink-faint
                         hover:border-print-accent hover:text-print-accent"
            >
              +
            </button>
          )}
        </>
      )}
    </div>
  );
}

/**
 * LoupeSidebar — the right column of the loupe. Per-exposure facts, a
 * keep/drop verdict, the indexing-representative control, the sample-tags
 * read+write rail (#159 / #207), and a keyboard legend.
 *
 * Tag editing lives here because single-sample tagging is loupe-scoped in the
 * mockup (`sample-table.html:914-917`). The corpus contact sheet keeps the
 * inert `+ tag` placeholder — its `loupe-tag-add` analogue is intentionally
 * the only on-ramp until the contact-sheet tag rail gets its own affordance
 * (out of scope for #207).
 */
export function LoupeSidebar({
  exposure,
  exposures,
  sample,
  signalLevel,
  onDropToggle,
  onSetRepresentative,
}: Props): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  const frameIndex = exposures.findIndex((e) => e.id === exposure.id);
  const frameLabel =
    frameIndex >= 0 ? `${frameIndex + 1} of ${exposures.length}` : "—";
  // Display an existing rejection reason if one was authored elsewhere (the
  // Inspect card / a future culling surface). The loupe never *authors* a
  // reason — its drop is a plain status toggle (spec §6).
  const rejectionReason = exposure.tags.find(
    (t) => t.key === "rejection_reason",
  )?.value;

  return (
    <aside data-testid="loupe-sidebar" className="flex flex-col gap-5">
      {/* This exposure — R3-S02 (#256): the meta-list is instrument-facing,
          mirroring the mockup's `frame / integration / collected / signal`
          (lowercase). The schema nouns "Filename" + "Kind"
          (averaged/background_subtracted) are gone — they leaked SQLite at a
          user who thinks in "2s exposure at 14:23". `integration` + `collected`
          stub to "—" until the backend fields are plumbed (mockup ships mock
          data here too). R2-M13: the redundant Status row stays dropped — the
          verdict card below carries that fact. */}
      <section>
        <SectionHeading>This exposure</SectionHeading>
        <div className="flex flex-col gap-1.5">
          <MetaRow label="frame" value={frameLabel} testid="loupe-meta-frame" />
          <MetaRow label="integration" value="—" testid="loupe-meta-integration" />
          <MetaRow label="collected" value="—" testid="loupe-meta-collected" />
          {/* M-8: signal-strength meter — peak-count proxy (see Props). */}
          <div className="flex items-center justify-between font-mono text-[11.5px]">
            <span className="text-ink-faint">signal</span>
            <SignalMeter level={signalLevel} />
          </div>
        </div>
      </section>

      {/* Verdict */}
      <section
        data-testid="loupe-verdict"
        className="flex items-center gap-3 rounded border border-hair-strong bg-paper-sunk p-3"
      >
        <span
          data-testid="loupe-kept-dot"
          className={[
            "h-2.5 w-2.5 shrink-0 rounded-full",
            // T-4/T-5: kept = sage success token; dropped = terracotta accent.
            isRejected ? "bg-print-accent" : "bg-success",
          ].join(" ")}
        />
        <div className="flex-1">
          <div data-testid="loupe-verdict-state" className="text-[13px] font-bold text-ink">
            {isRejected ? "Dropped" : "Kept"}
          </div>
          <div className="text-[10.5px] text-ink-faint">
            {isRejected
              ? (rejectionReason ?? "Dropped from this sample.")
              : "Everything is kept until you drop it."}
          </div>
        </div>
        <button
          data-testid="loupe-drop-toggle"
          onClick={onDropToggle}
          className="rounded border border-hair-strong bg-paper px-2.5 py-1.5
                     text-[11.5px] font-semibold text-ink hover:bg-paper-sunk"
        >
          {isRejected ? "Restore" : "Drop"}
          {/* R3-S03 (#256): the mono X keycap surfaces the drop/restore
              keyboard shortcut from the right rail, matching CullBar's
              "Drop X" and the footer-legend keycap idiom. */}
          <span className="ml-1 font-mono text-[10px] opacity-60">X</span>
        </button>
      </section>

      {/* Representative (R2-M12: when active, swap the neutral hair border for
          a terracotta-tinged border per the mockup's `.rep-box.is-rep`). */}
      <section
        data-testid="loupe-rep"
        data-is-rep={isRepresentative ? "true" : undefined}
        className={[
          "rounded border p-3",
          isRepresentative ? "border-print-accent/40" : "border-hair-strong",
        ].join(" ")}
      >
        {isRepresentative ? (
          <div className="flex items-center gap-2 text-xs font-bold text-print-accent">
            <span className="h-2 w-2 rounded-full bg-print-accent" />
            Representative for indexing
          </div>
        ) : (
          <>
            <div className="text-[11.5px] text-ink-soft">
              One exposure per sample carries forward to the Index stage.
              Pick the cleanest, strongest frame.
            </div>
            <button
              data-testid="loupe-set-representative"
              onClick={onSetRepresentative}
              className="mt-2 rounded border border-hair-strong bg-paper px-2.5 py-1.5
                         text-[11.5px] font-semibold text-ink hover:bg-paper-sunk"
            >
              Set as representative
            </button>
          </>
        )}
      </section>

      {/* Sample tags — read+write, routes through the corpus tag mutators. */}
      <section>
        <SectionHeading>Sample tags</SectionHeading>
        <LoupeTagsEditor sample={sample} />
      </section>

      {/* Keyboard legend */}
      <section>
        <SectionHeading>Keys</SectionHeading>
        <div className="flex flex-col gap-1 text-[11px] text-ink-faint">
          <div>
            <kbd className="font-mono">←</kbd> <kbd className="font-mono">→</kbd> flip frames
          </div>
          <div><kbd className="font-mono">X</kbd> drop / restore</div>
          <div><kbd className="font-mono">R</kbd> set representative</div>
          <div><kbd className="font-mono">Esc</kbd> back to the sheet</div>
        </div>
      </section>
    </aside>
  );
}
