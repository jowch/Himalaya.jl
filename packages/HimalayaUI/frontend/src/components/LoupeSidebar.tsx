import type { CorpusSample, Exposure } from "../api";

interface Props {
  /** The active exposure. */
  exposure: Exposure;
  /** All exposures for the sample — drives the "frame N of M" position. */
  exposures: Exposure[];
  /** The sample, for the read-only tags section. */
  sample: CorpusSample;
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

function SectionHeading({ children }: { children: string }): JSX.Element {
  return (
    <div className="mb-2 text-[10px] font-bold uppercase tracking-wide text-ink-faint">
      {children}
    </div>
  );
}

/**
 * LoupeSidebar — the right column of the loupe. Per-exposure facts, a
 * keep/drop verdict, the indexing-representative control, a read-only
 * sample-tags section, and a keyboard legend.
 *
 * The tags section is intentionally read-only: the corpus sample-tag
 * add/remove round-trip is #159 (I1.3), out of scope for #161.
 */
export function LoupeSidebar({
  exposure,
  exposures,
  sample,
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
      {/* This exposure */}
      <section>
        <SectionHeading>This exposure</SectionHeading>
        <div className="flex flex-col gap-1.5">
          <MetaRow
            label="Filename"
            value={exposure.filename ?? "—"}
            testid="loupe-meta-filename"
          />
          <MetaRow label="Kind" value={exposure.kind} testid="loupe-meta-kind" />
          <MetaRow label="Frame" value={frameLabel} testid="loupe-meta-frame" />
          <MetaRow
            label="Status"
            value={exposure.status ?? "pending"}
            testid="loupe-meta-status"
          />
        </div>
      </section>

      {/* Verdict */}
      <section
        data-testid="loupe-verdict"
        className="flex items-center gap-3 rounded border border-border bg-bg-subtle p-3"
      >
        <span
          className={[
            "h-2.5 w-2.5 shrink-0 rounded-full",
            isRejected ? "bg-accent" : "bg-success",
          ].join(" ")}
        />
        <div className="flex-1">
          <div className="text-[13px] font-bold text-ink">
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
          className="rounded border border-border bg-paper px-2.5 py-1.5
                     text-[11.5px] font-semibold text-ink hover:bg-bg-subtle"
        >
          {isRejected ? "Restore" : "Drop"}
        </button>
      </section>

      {/* Representative */}
      <section
        data-testid="loupe-rep"
        className="rounded border border-border p-3"
      >
        {isRepresentative ? (
          <div className="flex items-center gap-2 text-xs font-bold text-accent">
            <span className="h-2 w-2 rounded-full bg-accent" />
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
              className="mt-2 rounded border border-border bg-paper px-2.5 py-1.5
                         text-[11.5px] font-semibold text-ink hover:bg-bg-subtle"
            >
              Set as representative
            </button>
          </>
        )}
      </section>

      {/* Sample tags — read-only (editing is #159) */}
      <section>
        <SectionHeading>Sample tags</SectionHeading>
        <div data-testid="loupe-tags" className="flex flex-wrap gap-1">
          {sample.tags.length === 0 ? (
            <span className="text-[11.5px] text-ink-faint">No tags yet</span>
          ) : (
            sample.tags.map((tag) => (
              <span
                key={tag.id}
                className="inline-flex items-center rounded-full border border-border
                           bg-bg-subtle px-2 py-0.5 text-xs text-ink-soft"
              >
                {tag.key}: {tag.value}
              </span>
            ))
          )}
        </div>
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
