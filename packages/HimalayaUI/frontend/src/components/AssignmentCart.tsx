import { phaseColor } from "../phases";
import { Card, IconButton, ScoreBar, Button } from "./ui";
import { SegmentedControl } from "./ui/SegmentedControl";
import { seriesRatio } from "../lib/seriesRatio";
import type { Assignment, AssignmentState, IndexEntry } from "../api";

/**
 * AssignmentCart (Plan D Task D-4) — the redesigned phase-assignment cart.
 *
 * Renders the durable 3-state assignment:
 *   - "indexed" with N members → one phase block per member (× remove in the
 *     block header, serif name + mono score + lattice/reflections + score bar),
 *     or the empty "form factor / unindexed" declaration when 0 members.
 *   - "form_factor" / "null" → a single declaration chip.
 *
 * A SegmentedControl toggles the 3-state. A contextual Bonnet note appears only
 * when substantive: (a) a bicontinuous cubic is in the cart, no second cubic,
 * and peaks remain unexplained → suggest a coexisting cubic at the Bonnet
 * lattice; (b) two cubics coexist → report the lattice-ratio ≈ Bonnet
 * consistency. A `+ custom index…` footer Button opens the modal.
 *
 * Visual contract: docs/redesign-mockups/2026-05-29-focus-plot.html (.assign).
 * Appearance is threaded via phaseColor() through style attributes (the
 * design-guard-safe path, matching PhaseCallBlock); className is placement-only.
 */

const CUBICS = new Set(["Pn3m", "Im3m", "Ia3d"]);

function formatScore(s: number | null | undefined): string {
  return s != null ? s.toFixed(2) : "—";
}

interface AssignmentBlockProps {
  index: IndexEntry;
  latticeUnit: string;
  onRemove: () => void;
}

function AssignmentBlock({ index, latticeUnit, onRemove }: AssignmentBlockProps): JSX.Element {
  const color = phaseColor(index.phase);
  const ratio = seriesRatio(index.phase, index.peaks.map((p) => p.ratio_position));
  return (
    <div data-testid={`assignment-block-${index.id}`} className="px-4 py-3">
      <div className="flex items-baseline gap-2">
        <span
          className="flex-1 min-w-0 font-serif font-medium leading-none tracking-tight text-xl"
          style={{ color }}
        >
          {index.phase}
        </span>
        <span className="font-mono text-xs font-bold tabular-nums" style={{ color }}>
          {formatScore(index.score)}
        </span>
        <IconButton
          label={`remove ${index.phase}`}
          tone="ghost"
          className="self-center"
          onClick={onRemove}
        >
          <span aria-hidden className="text-base leading-none">×</span>
        </IconButton>
      </div>
      <div className="mt-1.5 font-mono text-sm text-ink-soft">
        {index.lattice_d != null && (
          <>{`a = ${index.lattice_d.toFixed(0)} ${latticeUnit}`}&nbsp; ·&nbsp; </>
        )}
        {index.peaks.length} reflections
      </div>
      {index.score != null && (
        <ScoreBar value={index.score} phase={index.phase} size="bar" className="mt-2" />
      )}
      {ratio && (
        <div className="mt-1.5 text-xs text-ink-faint">
          series&nbsp;&nbsp;
          <span className="font-mono font-semibold text-ink-soft">{ratio}</span>
        </div>
      )}
    </div>
  );
}

/**
 * Build the contextual Bonnet note text, or null when not substantive.
 *  - one cubic in, no second cubic, peaks unexplained → suggestion
 *  - two cubics in → consistency (lattice ratio ≈ Bonnet ⭙ 1.279)
 */
function bonnetNote(active: IndexEntry[], peakCount: number): JSX.Element | null {
  const cubics = active.filter((ix) => CUBICS.has(ix.phase));
  if (cubics.length === 1) {
    const claimed = new Set<number>();
    for (const ix of active) for (const p of ix.peaks) claimed.add(p.peak_id);
    const leftover = Math.max(0, peakCount - claimed.size);
    if (leftover > 0) {
      return (
        <>
          <b className="text-ink-soft font-semibold">{leftover} peak{leftover > 1 ? "s" : ""}</b>
          {" "}unexplained. A coexisting bicontinuous cubic at the Gauss–Bonnet
          lattice would account for them.
        </>
      );
    }
    return null;
  }
  if (cubics.length >= 2) {
    const [a, b] = cubics;
    const ratio = (a!.lattice_d != null && b!.lattice_d != null && a!.lattice_d > 0)
      ? (Math.max(a!.lattice_d, b!.lattice_d) / Math.min(a!.lattice_d, b!.lattice_d))
      : null;
    return (
      <>
        {cubics.map((c) => c.phase).join(" + ")} coexistence.
        {ratio != null && (
          <> Lattice ratio{" "}
            <span className="font-mono text-ink-soft">{ratio.toFixed(3)}</span>
            {" "}≈ Bonnet{" "}
            <span className="font-mono text-ink-soft">1.279</span>. Consistent.
          </>
        )}
      </>
    );
  }
  return null;
}

const STATE_OPTIONS = [
  { value: "indexed" as const, label: "Indexed" },
  { value: "form_factor" as const, label: "Form factor" },
  { value: "null" as const, label: "No scattering" },
];

export interface AssignmentCartProps {
  assignment: Assignment | undefined;
  indices: IndexEntry[];
  /** Total observed-peak count (for the Bonnet leftover heuristic). */
  peakCount: number;
  latticeUnit: string;
  onRemovePhase: (indexId: number) => void;
  onSetState: (state: AssignmentState) => void;
  onOpenCustom: () => void;
}

export function AssignmentCart({
  assignment, indices, peakCount, latticeUnit,
  onRemovePhase, onSetState, onOpenCustom,
}: AssignmentCartProps): JSX.Element {
  const state: AssignmentState = assignment?.state ?? "indexed";
  const memberIds = new Set(assignment?.members ?? []);
  const active = indices
    .filter((ix) => memberIds.has(ix.id))
    .sort((a, b) => (b.score ?? 0) - (a.score ?? 0));

  const note = state === "indexed" ? bonnetNote(active, peakCount) : null;

  return (
    <div className="flex flex-col gap-2.5">
      <SegmentedControl
        options={STATE_OPTIONS}
        value={state}
        onChange={onSetState}
        role="radiogroup"
        aria-label="Assignment state"
        testId="assignment-state"
      />

      <Card className="overflow-hidden">
        {state !== "indexed" ? (
          <div
            data-testid="assignment-declaration"
            className="px-4 py-4 text-sm text-ink-soft"
          >
            {state === "form_factor"
              ? "Form factor — diffuse scattering, no lattice to index."
              : "No scattering — this exposure has no usable signal."}
          </div>
        ) : active.length === 0 ? (
          <div
            data-testid="assignment-empty"
            className="px-4 py-4 text-xs text-ink-faint leading-relaxed"
          >
            No phase yet — this is form factor / unindexed. Evaluate the
            candidates and add the ones that make sense.
          </div>
        ) : (
          <div className="divide-y divide-hair">
            {active.map((ix) => (
              <AssignmentBlock
                key={ix.id}
                index={ix}
                latticeUnit={latticeUnit}
                onRemove={() => onRemovePhase(ix.id)}
              />
            ))}
          </div>
        )}

        {/* + custom index… footer — a real Button so it reads as one on hover
            (fixes the mockup-noted hover gap). Border-top via a placement
            divider class; appearance comes from the Button primitive. */}
        <div className="border-t border-dashed border-hair-strong">
          <Button
            variant="ghost"
            data-testid="custom-index-footer"
            className="w-full justify-start text-xs font-semibold"
            onClick={onOpenCustom}
          >
            <span aria-hidden className="text-accent font-bold">+ </span>custom index…
          </Button>
        </div>
      </Card>

      {note && (
        <p
          data-testid="bonnet-note"
          className="text-xs leading-[1.55] text-ink-faint"
        >
          {note}
        </p>
      )}
    </div>
  );
}
