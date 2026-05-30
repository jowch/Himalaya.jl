import type { ReactNode } from "react";
import { Skeleton } from "boneyard-js/react";
import { useIndices, useAssignment, useAddAssignmentPhase, useRemoveAssignmentPhase, useDeleteIndex, useExperiment } from "../queries";
import { useAppState } from "../state";
import { phaseColor } from "../phases";
import { Card, HintText, IconButton, ScoreBar, Kicker } from "./ui";
import { StaleIndicesBanner } from "./StaleIndicesBanner";
import { SpeculativeBuilder } from "./SpeculativeBuilder";
import { latticeUnitFromQUnits } from "../lib/units";
import { seriesRatio } from "../lib/seriesRatio";
import type { IndexEntry } from "../api";
import { deriveActiveIndices } from "../lib/assignment";

function formatScore(s: number | null | undefined): string {
  return s != null ? s.toFixed(2) : "—";
}

function ratioOf(index: IndexEntry): string {
  return seriesRatio(index.phase, index.peaks.map((p) => p.ratio_position));
}

// ── Phase-call output block (R4 L-9) ─────────────────────────────────────────
// One card per phase currently in the active set: serif phase name + mono
// score, lattice + peak count, a score bar tinted in the phase colour, and the
// √N ratio series. Mirrors mockup `.phasecall` / `.pc-block`.

interface PhaseCallBlockProps {
  index: IndexEntry;
  latticeUnit: string;
}

function PhaseCallBlock({ index, latticeUnit }: PhaseCallBlockProps): JSX.Element {
  const color = phaseColor(index.phase);
  const score = index.score ?? 0;
  const ratio = ratioOf(index);
  return (
    <div data-testid={`phase-call-block-${index.id}`} className="px-4 py-3">
      <div className="flex items-baseline justify-between gap-2">
        <span
          className="font-serif font-medium leading-none tracking-tight text-xl"
          style={{ color }}
        >
          {index.phase}
        </span>
        <span className="font-mono text-xs font-bold tabular-nums" style={{ color }}>
          {formatScore(index.score)}
        </span>
      </div>
      <div className="mt-1.5 font-mono text-sm text-ink-soft">
        {index.lattice_d != null && (
          <>{`a = ${index.lattice_d.toFixed(0)} ${latticeUnit}`}&nbsp; ·&nbsp; </>
        )}
        {index.peaks.length} peaks
      </div>
      <ScoreBar value={score} phase={index.phase} size="bar" className="mt-2" />
      {ratio && (
        <div className="mt-1.5 text-xs text-ink-faint">
          series&nbsp;&nbsp;
          <span className="font-mono font-semibold text-ink-soft">{ratio}</span>
        </div>
      )}
    </div>
  );
}

// ── Candidate row (R4 L-10) ──────────────────────────────────────────────────
// A multi-select checkbox row: the active set is a SET, not a single pick.
// Clicking toggles membership; hovering previews on the trace (phase colour,
// L-11). Mirrors mockup `.cand` / `.c-mark`.

interface CandidateRowProps {
  index: IndexEntry;
  inCall: boolean;
  onToggle: () => void;
  onHover: () => void;
  onLeave: () => void;
  onDelete?: () => void;
}

function CandidateRow({ index, inCall, onToggle, onHover, onLeave, onDelete }: CandidateRowProps): JSX.Element {
  const color = phaseColor(index.phase);
  const ratio = ratioOf(index);
  const score = index.score ?? 0;
  return (
    <div
      data-testid={`candidate-row-${index.id}`}
      data-alternative-id={index.id}
      data-active={inCall || undefined}
      className="group flex items-center gap-3 rounded-lg border px-3 py-2.5 transition-colors cursor-pointer focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
      style={{
        // Hover-preview is phase-coloured (L-11): tint the border in the
        // phase's own hue on hover; membership uses the terracotta accent.
        borderColor: inCall
          ? `color-mix(in oklab, var(--color-accent) 42%, var(--color-hair))`
          : `var(--color-hair)`,
        background: inCall
          ? `color-mix(in oklab, var(--color-accent) 6%, var(--color-plate))`
          : `var(--color-plate)`,
        ["--pc" as string]: color,
      }}
      role="checkbox"
      aria-checked={inCall}
      aria-label={index.phase}
      tabIndex={0}
      onClick={onToggle}
      onKeyDown={(e) => {
        if (e.key === " " || e.key === "Enter") { e.preventDefault(); onToggle(); }
      }}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
      onFocus={onHover}
      onBlur={onLeave}
    >
      {/* checkbox mark */}
      <span
        aria-hidden
        className="flex h-4 w-4 shrink-0 items-center justify-center rounded border-[1.5px] text-paper"
        style={{
          borderColor: inCall ? "var(--color-accent)" : "var(--color-hair-strong)",
          background: inCall ? "var(--color-accent)" : "transparent",
        }}
      >
        {inCall && (
          <svg viewBox="0 0 16 16" className="h-3 w-3">
            <path d="M4 8.2l2.7 2.7 5.2-5.6" fill="none" stroke="currentColor"
                  strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        )}
      </span>

      {/* body */}
      <div className="min-w-0 flex-1">
        <div className="flex items-center gap-1.5">
          <span className="font-mono text-base font-bold" style={{ color }}>
            {index.phase}
          </span>
          {/* Plan D F3: Bonnet badge — a coexisting cubic whose lattice matches
              the Gauss–Bonnet ratio (IndexEntry.bonnet.consistent). The single
              automation surfaced in the candidate list. */}
          {index.bonnet?.consistent && (
            <span
              data-testid={`bonnet-badge-${index.id}`}
              className="inline-flex items-center gap-0.5 rounded-full border border-accent/40 bg-accent/10 px-1.5 text-xs font-bold uppercase tracking-wide text-accent"
              title={`Coexisting cubic at the Gauss–Bonnet lattice (predicted a ≈ ${index.bonnet.predicted_a.toFixed(0)})`}
            >
              <span aria-hidden>⭙</span> Bonnet
            </span>
          )}
        </div>
        <div className="mt-0.5 text-xs text-ink-faint">
          explains {index.peaks.length} peaks{inCall ? " · in the call" : ""}
          {ratio && <span className="font-mono"> · {ratio}</span>}
        </div>
      </div>

      {/* score */}
      <div className="text-right font-mono">
        <div className="text-base font-bold text-ink tabular-nums">{formatScore(index.score)}</div>
        <ScoreBar value={score} phase={index.phase} size="compact" className="mt-1" />
      </div>

      {/* speculative delete */}
      {onDelete && (
        <IconButton
          label={`Delete speculative index ${index.id}`}
          tone="danger"
          title="Delete this speculative index"
          data-testid={`spec-delete-${index.id}`}
          onClick={(e) => { e.stopPropagation(); onDelete(); }}
          className="shrink-0"
        >
          <svg viewBox="0 0 16 16" className="h-3.5 w-3.5" aria-hidden="true">
            <path
              d="M3.5 4.5h9M6 4.5V3.2a1 1 0 011-1h2a1 1 0 011 1v1.3M5 4.5l.55 8a1 1 0 001 .93h2.9a1 1 0 001-.93l.55-8"
              fill="none" stroke="currentColor" strokeWidth="1.3"
              strokeLinecap="round" strokeLinejoin="round"
            />
          </svg>
        </IconButton>
      )}
    </div>
  );
}

// ── Rail section heading ─────────────────────────────────────────────────────

function RailHead({ children }: { children: ReactNode }): JSX.Element {
  return <Kicker tone="faint">{children}</Kicker>;
}

// ── Panel ────────────────────────────────────────────────────────────────────

const FIXTURE_INDICES: IndexEntry[] = [
  { id:1, exposure_id:0, phase:"Pn3m",     basis:0.15, score:0.91, r_squared:0.995,
    lattice_d:197, ngc:-1.51, status:"candidate", kind:"auto", inputs_hash:null,
    peaks:[{ peak_id:1, ratio_position:1, residual:0.001, q_observed:0.045 },
           { peak_id:2, ratio_position:2, residual:0.002, q_observed:0.055 },
           { peak_id:3, ratio_position:3, residual:0.002, q_observed:0.064 }],
    predicted_q:[0.045,0.055,0.064] },
  { id:2, exposure_id:0, phase:"Lamellar", basis:0.12, score:0.55, r_squared:0.960,
    lattice_d:61, ngc:null, status:"candidate", kind:"auto", inputs_hash:null,
    peaks:[{ peak_id:4, ratio_position:1, residual:0.004, q_observed:0.103 }],
    predicted_q:[0.103,0.206] },
];

const PHASE_PANEL_FIXTURE = (
  <div className="flex-1 overflow-y-auto p-3 flex flex-col gap-5">
    <div className="flex flex-col gap-2.5">
      <RailHead>Phase call</RailHead>
      <Card className="overflow-hidden">
        <PhaseCallBlock index={FIXTURE_INDICES[0]!} latticeUnit="Å" />
      </Card>
    </div>
    <div className="flex flex-col gap-2.5">
      <RailHead>Candidate indexings</RailHead>
      <div className="flex flex-col gap-[7px]">
        <CandidateRow index={FIXTURE_INDICES[0]!} inCall onToggle={() => {}} onHover={() => {}} onLeave={() => {}} />
        <CandidateRow index={FIXTURE_INDICES[1]!} inCall={false} onToggle={() => {}} onHover={() => {}} onLeave={() => {}} />
      </div>
    </div>
  </div>
);

export interface PhasePanelProps {
  exposureId: number | undefined;
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const indicesQ = useIndices(exposureId);
  const assignmentQ = useAssignment(exposureId);
  const experimentQ = useExperiment(activeExperimentId ?? 0);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const setPreviewIndex = useAppState((s) => s.setPreviewIndex);
  // Plan D-8: candidate toggles drive the assignment cart natively (the legacy
  // group dual-write remains live on the backend until D-10, but the frontend
  // no longer touches /groups).
  const addMember    = useAddAssignmentPhase(exposureId ?? 0);
  const removeMember = useRemoveAssignmentPhase(exposureId ?? 0);
  const deleteIndex  = useDeleteIndex(exposureId ?? 0);
  const builder      = useAppState((s) => s.speculativeBuilder);
  const openBuilder  = useAppState((s) => s.openSpeculativeBuilder);
  const closeBuilder = useAppState((s) => s.closeSpeculativeBuilder);

  const qUnits        = experimentQ.data?.q_units ?? null;
  const latticeUnit   = latticeUnitFromQUnits(qUnits);

  if (exposureId === undefined) {
    return (
      <div className="p-4">
        <HintText>No exposure selected.</HintText>
      </div>
    );
  }

  const indices = (indicesQ.data ?? []).slice().sort(
    (a, b) => (b.score ?? 0) - (a.score ?? 0),
  );
  // Active member set sourced from the durable assignment cart.
  const activeIndices   = deriveActiveIndices(assignmentQ.data, indices);
  const memberIds       = new Set(activeIndices.map((ix) => ix.id));
  const speculatives    = indices.filter((ix) => ix.kind === "speculative");
  const auto            = indices.filter((ix) => ix.kind !== "speculative");
  // Phase call = the active set, ordered by lowest claimed q (mockup ordering).
  const inCall = indices.filter((ix) => memberIds.has(ix.id)).sort((a, b) => {
    const aq = Math.min(...a.peaks.map((p) => p.q_observed), Infinity);
    const bq = Math.min(...b.peaks.map((p) => p.q_observed), Infinity);
    return aq - bq;
  });

  const toggle = (ix: IndexEntry): void => {
    if (memberIds.has(ix.id)) removeMember.mutate(ix.id);
    else addMember.mutate(ix.id);
  };

  // Plan D-7: hover a candidate → preview on the PLOT only (ghost combs + trace
  // highlight); the cart is left untouched. NEVER a mutator/SSE event. Clear on
  // leave/blur so a stale ghost never masks the real cart.
  const previewOn  = (ix: IndexEntry): void => { setHoveredIndex(ix.id); setPreviewIndex(ix.id); };
  const previewOff = (): void => { setHoveredIndex(undefined); setPreviewIndex(undefined); };

  return (
    <div className="flex flex-col h-full min-h-0">

      {/* ── Sticky header ── */}
      {/* R4-N1 (#209): the explanatory sentence "Check every phase that is
          present…" appeared here AND in the rail-note paragraph below the
          candidates ~3 lines apart. Mockup `focus-workspace.html` has only the
          label "Index choices" in the rail head, with the explanatory note
          below the candidates — so we keep the rail-note (which has more
          context) and drop the subtitle here. */}
      <div className="card-header">
        <div className="flex flex-col justify-center min-w-0">
          <div className="text-title tracking-tight">Index choices</div>
        </div>
      </div>

      {/* ── Stale-indices banner (shown when re-analysis is needed) ── */}
      <div className="px-3 pt-2">
        <StaleIndicesBanner exposureId={exposureId} />
      </div>

      {/* ── Scrollable rail ── */}
      <Skeleton
        name="phase-panel"
        className="flex-1 min-h-0 flex flex-col"
        loading={indicesQ.isLoading || assignmentQ.isLoading}
        stagger={50}
        transition={200}
        fixture={PHASE_PANEL_FIXTURE}
        fallback={<div className="p-4"><HintText>Loading phase assignments…</HintText></div>}
      >
        <div className="flex-1 overflow-y-auto p-3 flex flex-col gap-5">

          {/* Phase call — the output */}
          <div className="flex flex-col gap-2.5">
            <RailHead>Phase call</RailHead>
            <Card className="overflow-hidden">
              {inCall.length === 0 ? (
                <div data-testid="phase-call-empty" className="px-4 py-4 text-xs text-ink-faint">
                  No phase assigned; every peak is unindexed. Check a candidate below.
                </div>
              ) : (
                <>
                  {inCall.length > 1 && (
                    <Kicker
                      tone="faint"
                      data-testid="coexistence-tag"
                      className="px-4 pt-2.5"
                    >
                      Coexistence · {inCall.length} phases
                    </Kicker>
                  )}
                  <div className="divide-y divide-hair">
                    {inCall.map((ix) => (
                      <PhaseCallBlock key={ix.id} index={ix} latticeUnit={latticeUnit} />
                    ))}
                  </div>
                </>
              )}
            </Card>
          </div>

          {/* Candidate indexings — the multi-select active set */}
          <div className="flex flex-col gap-2.5">
            <RailHead>Candidate indexings</RailHead>
            {auto.length === 0 ? (
              <HintText>No candidate indexings.</HintText>
            ) : (
              <div className="flex flex-col gap-[7px]">
                {auto.map((ix) => (
                  <CandidateRow
                    key={ix.id}
                    index={ix}
                    inCall={memberIds.has(ix.id)}
                    onToggle={() => toggle(ix)}
                    onHover={() => previewOn(ix)}
                    onLeave={previewOff}
                  />
                ))}
              </div>
            )}
            <p className="text-sm leading-[1.55] text-ink-faint">
              Check every phase that is present; a sample can hold more than one.
              Candidates that explain the same peaks swap; independent phases coexist.
            </p>
          </div>

          {/* Speculative — user-built sub-minpeaks indices. R3-F02: tucked
              below-the-fold behind a native <details> disclosure so the rail
              reads as the calm two-section output (Phase call + Candidates).
              Auto-open when speculatives exist so the user's own builds stay
              visible; collapsed by default otherwise. The "+ Add speculative"
              CTA lives inside the disclosure body, so it only surfaces when the
              section is open (the browser hides closed-<details> content). */}
          <details
            data-testid="speculative-disclosure"
            open={speculatives.length > 0}
            className="flex flex-col gap-2.5"
          >
            <summary className="cursor-pointer list-none">
              <RailHead>Speculative</RailHead>
            </summary>
            {speculatives.length > 0 && (
              <div className="mt-2.5 flex flex-col gap-[7px]">
                {speculatives.map((ix) => (
                  <CandidateRow
                    key={ix.id}
                    index={ix}
                    inCall={memberIds.has(ix.id)}
                    onToggle={() => toggle(ix)}
                    onHover={() => previewOn(ix)}
                    onLeave={previewOff}
                    onDelete={() => deleteIndex.mutate(ix.id)}
                  />
                ))}
              </div>
            )}
            <button
              type="button"
              data-testid="add-speculative-button"
              className="mt-2.5 w-full text-xs text-ink-faint border border-dashed border-hair rounded-md py-1.5 hover:text-ink hover:bg-paper-sunk transition-colors"
              onClick={() => openBuilder(exposureId)}
            >
              + Add speculative
            </button>
          </details>

        </div>
      </Skeleton>

      {builder && builder.exposureId === exposureId && (
        <SpeculativeBuilder exposureId={exposureId} onClose={closeBuilder} />
      )}
    </div>
  );
}
