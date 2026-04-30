import { Skeleton } from "boneyard-js/react";
import { useIndices, useGroups, useAddIndexToGroup, useRemoveIndexFromGroup, useDeleteIndex, useExperiment } from "../queries";
import { useAppState } from "../state";
import { phaseColor, CUBIC_PHASES } from "../phases";
import { HintText } from "./ui";
import { StaleIndicesBanner } from "./StaleIndicesBanner";
import { SpeculativeBuilder } from "./SpeculativeBuilder";
import { latticeUnitFromQUnits, inverseSquareUnits } from "../lib/units";
import type { GroupEntry, IndexEntry } from "../api";

const R2_THRESHOLD = 0.98;

function activeGroup(groups: GroupEntry[]): GroupEntry | undefined {
  return groups.find((g) => g.active);
}

function formatLattice(d: number | null): string {
  return d != null ? d.toFixed(2) : "—";
}

function formatR2(r: number | null): string {
  return r != null ? r.toFixed(3) : "—";
}

function formatScore(s: number | null | undefined): string {
  return s != null ? s.toFixed(2) : "—";
}

/**
 * κ ranges over many orders of magnitude depending on the experiment's
 * length unit (≈10⁻² nm⁻² ≈ 10⁻⁴ Å⁻²). Toggle to scientific notation when
 * a fixed-point representation would lose all the leading sig figs.
 */
function formatKappa(k: number): string {
  const abs = Math.abs(k);
  if (abs === 0) return "0";
  if (abs < 0.01) return k.toExponential(2);
  return k.toFixed(3);
}

// ── Individual index card ────────────────────────────────────────────────────

interface IndexCardProps {
  index: IndexEntry;
  isActive: boolean;
  onAction: () => void;
  onDelete?: () => void;
  onHover?: () => void;
  onLeave?: () => void;
  /** Lattice unit derived from the experiment's q_units (e.g. "Å", "nm"). */
  latticeUnit: string;
  /** Inverse-square form for κ display (e.g. "Å⁻²", "nm⁻²"). */
  curvatureUnit: string;
  /** Forwarded to the li for E2E selectors */
  "data-alternative-id"?: number;
}

function IndexCard({ index, isActive, onAction, onDelete, onHover, onLeave, latticeUnit, curvatureUnit, "data-alternative-id": altId }: IndexCardProps): JSX.Element {
  const color = phaseColor(index.phase);
  const isSpeculative = index.kind === "speculative";
  // R²-gate dimming: speculative indices bypass the gate entirely (a 2-peak fit
  // is R²=1 by construction, so the gate is meaningless — the "speculative"
  // label is the warning instead).
  const lowR2 = !isActive && !isSpeculative && index.r_squared != null && index.r_squared < R2_THRESHOLD;

  return (
    <li
      data-index-id={index.id}
      data-alternative-id={altId}
      data-active={isActive || undefined}
      data-low-r2={lowR2 || undefined}
      data-speculative={isSpeculative || undefined}
      className={[
        "grid items-stretch rounded-lg overflow-hidden transition-all",
        isSpeculative
          ? "border border-dashed border-border-soft"
          : "border border-border-soft",
        lowR2 ? "opacity-40" : "",
      ].join(" ")}
      style={{
        gridTemplateColumns: "3px 1fr auto",
        background: isActive
          ? `color-mix(in oklab, ${color} 6%, transparent)`
          : undefined,
        borderColor: isActive ? `color-mix(in oklab, ${color} 28%, var(--color-border-soft))` : undefined,
      }}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
    >
      {/* Left color bar */}
      <div style={{ background: color }} />

      {/* Main content */}
      <div className="px-2.5 py-2 flex flex-col gap-1 min-w-0">
        {/* Primary row: phase chip + lattice param */}
        <div className="flex items-center gap-2 min-w-0">
          <span
            className="text-data-strong px-1.5 py-0.5 rounded-sm border shrink-0"
            style={{
              color,
              background: `color-mix(in oklab, ${color} 10%, transparent)`,
              borderColor: `color-mix(in oklab, ${color} 35%, transparent)`,
            }}
          >
            {index.phase}
          </span>
          {index.lattice_d != null && (
            <span className="text-data truncate min-w-0">
              <span className="text-fg-dim">a =</span>{" "}
              {formatLattice(index.lattice_d)}{" "}
              <span className="text-fg-dim text-xs">{latticeUnit}</span>
            </span>
          )}
        </div>

        {/* Secondary row: score bar + R² + peak count */}
        <div className="flex items-center gap-3 font-mono text-xs text-fg-dim">
          <span className="flex items-center gap-1.5">
            <span>score</span>
            <span className="inline-block w-12 h-1.5 bg-bg-hover rounded-full overflow-hidden">
              <span
                data-score-bar
                className="block h-full"
                style={{
                  width: `${Math.round((index.score ?? 0) * 100)}%`,
                  background: color,
                }}
              />
            </span>
            <span className="text-fg-muted tabular-nums">{formatScore(index.score)}</span>
          </span>
          <span>
            R²{" "}
            <span className={index.r_squared != null && index.r_squared >= R2_THRESHOLD
              ? "text-fg-muted" : "text-fg-dim"}>
              {formatR2(index.r_squared)}
            </span>
          </span>
          {CUBIC_PHASES.has(index.phase) && index.ngc != null && (
            <span data-testid={`ngc-${index.id}`}>
              κ{" "}
              <span className="text-fg-muted tabular-nums">
                {formatKappa(index.ngc)}
              </span>{" "}
              <span className="text-fg-dim">{curvatureUnit}</span>
            </span>
          )}
          <span className="ml-auto px-1.5 py-0.5 border border-border-soft rounded-full text-xs text-fg-dim">
            {index.peaks.length} peaks
          </span>
        </div>
      </div>

      {/* Add / remove button + (speculative) delete */}
      <div className="flex flex-col border-l border-border-soft">
        <button
          className="flex-1 w-[34px] bg-transparent text-fg-dim hover:text-fg hover:bg-bg-hover transition-colors text-base font-semibold"
          onClick={onAction}
          aria-label={isActive ? `Remove index ${index.id}` : `Add index ${index.id}`}
        >
          {isActive ? "−" : "+"}
        </button>
        {onDelete && (
          <button
            className="flex-1 w-[34px] bg-transparent text-fg-dim hover:text-error hover:bg-bg-hover transition-colors text-xs border-t border-border-soft"
            onClick={onDelete}
            aria-label={`Delete speculative index ${index.id}`}
            data-testid={`spec-delete-${index.id}`}
            title="Delete this speculative index"
          >
            🗑
          </button>
        )}
      </div>
    </li>
  );
}

// ── Group heading ────────────────────────────────────────────────────────────

function GroupHead({ label, count }: { label: string; count: number }): JSX.Element {
  return (
    <div className="flex items-center justify-between mb-2 px-1">
      <span className="text-xs uppercase tracking-widest text-fg-dim font-semibold">
        {label}
      </span>
      <span className="text-xs text-fg-dim">{count}</span>
    </div>
  );
}

// ── Panel ────────────────────────────────────────────────────────────────────

const FIXTURE_INDICES: IndexEntry[] = [
  { id:1, exposure_id:0, phase:"Pn3m",     basis:0.15, score:0.91, r_squared:0.995,
    lattice_d:64.2, ngc:-1.51, status:"candidate", kind:"auto",
    peaks:[{ peak_id:1, ratio_position:1, residual:0.001, q_observed:0.15 },
           { peak_id:2, ratio_position:2, residual:0.002, q_observed:0.21 }],
    predicted_q:[0.15,0.21,0.26] },
  { id:2, exposure_id:0, phase:"Im3m",     basis:0.14, score:0.72, r_squared:0.981,
    lattice_d:57.1, ngc:-2.06, status:"candidate", kind:"auto",
    peaks:[{ peak_id:1, ratio_position:1, residual:0.003, q_observed:0.15 }],
    predicted_q:[0.15,0.22] },
  { id:3, exposure_id:0, phase:"Lamellar", basis:0.12, score:0.55, r_squared:0.960,
    lattice_d:52.4, ngc:null, status:"candidate", kind:"auto",
    peaks:[{ peak_id:2, ratio_position:1, residual:0.004, q_observed:0.21 }],
    predicted_q:[0.21,0.42] },
];

const PHASE_PANEL_FIXTURE = (
  <div className="flex-1 overflow-y-auto p-3 flex flex-col gap-4">
    <div>
      <GroupHead label="Active set" count={1} />
      <ul className="flex flex-col gap-1.5">
        <IndexCard key={1} index={FIXTURE_INDICES[0]!} isActive
                   latticeUnit="Å" curvatureUnit="Å⁻²" onAction={() => {}} />
      </ul>
    </div>
    <div>
      <GroupHead label="Candidates" count={2} />
      <ul className="flex flex-col gap-1.5">
        <IndexCard key={2} index={FIXTURE_INDICES[1]!} isActive={false}
                   latticeUnit="Å" curvatureUnit="Å⁻²" onAction={() => {}} />
        <IndexCard key={3} index={FIXTURE_INDICES[2]!} isActive={false}
                   latticeUnit="Å" curvatureUnit="Å⁻²" onAction={() => {}} />
      </ul>
    </div>
  </div>
);

export interface PhasePanelProps {
  exposureId: number | undefined;
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const indicesQ = useIndices(exposureId);
  const groupsQ  = useGroups(exposureId);
  const experimentQ = useExperiment(activeExperimentId ?? 0);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const active = (groupsQ.data && activeGroup(groupsQ.data)) ?? undefined;
  const addMember    = useAddIndexToGroup(exposureId ?? 0, active?.id ?? 0);
  const removeMember = useRemoveIndexFromGroup(exposureId ?? 0, active?.id ?? 0);
  const deleteIndex  = useDeleteIndex(exposureId ?? 0);
  const builder      = useAppState((s) => s.speculativeBuilder);
  const openBuilder  = useAppState((s) => s.openSpeculativeBuilder);
  const closeBuilder = useAppState((s) => s.closeSpeculativeBuilder);

  const qUnits        = experimentQ.data?.q_units ?? null;
  const latticeUnit   = latticeUnitFromQUnits(qUnits);
  const curvatureUnit = inverseSquareUnits(qUnits);

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
  const memberIds       = new Set(active?.members ?? []);
  const speculatives    = indices.filter((ix) => ix.kind === "speculative");
  const auto            = indices.filter((ix) => ix.kind !== "speculative");
  const activeMembers   = indices.filter((ix) =>  memberIds.has(ix.id));
  const alternatives    = auto.filter((ix) => !memberIds.has(ix.id));

  return (
    <div className="flex flex-col h-full min-h-0">

      {/* ── Sticky header ── */}
      <div className="card-header">
        <div className="flex flex-col justify-center min-w-0">
          <div className="text-title tracking-tight">Index choices</div>
          <div className="text-xs text-fg-dim leading-tight">
            Hover a candidate to preview peaks
          </div>
        </div>
      </div>

      {/* ── Stale-indices banner (shown when re-analysis is needed) ── */}
      <div className="px-3 pt-2">
        <StaleIndicesBanner exposureId={exposureId} />
      </div>

      {/* ── Scrollable list ── */}
      <Skeleton
        name="phase-panel"
        className="flex-1 min-h-0 flex flex-col"
        loading={indicesQ.isLoading || groupsQ.isLoading}
        stagger={50}
        transition={200}
        fixture={PHASE_PANEL_FIXTURE}
        fallback={<div className="p-4"><HintText>Loading phase assignments…</HintText></div>}
      >
        <div className="flex-1 overflow-y-auto p-3 flex flex-col gap-4">

          {/* Active set */}
          <div>
            <GroupHead label="Active set" count={activeMembers.length} />
            {activeMembers.length === 0 ? (
              <HintText>No indices in the active set.</HintText>
            ) : (
              <ul className="flex flex-col gap-1.5">
                {activeMembers.map((ix) => (
                  <IndexCard
                    key={ix.id}
                    index={ix}
                    isActive
                    latticeUnit={latticeUnit}
                    curvatureUnit={curvatureUnit}
                    onAction={() => { if (active) removeMember.mutate(ix.id); }}
                    onHover={() => setHoveredIndex(ix.id)}
                    onLeave={() => setHoveredIndex(undefined)}
                  />
                ))}
              </ul>
            )}
          </div>

          {/* Candidates */}
          <div>
            <GroupHead label="Candidates" count={alternatives.length} />
            {alternatives.length === 0 ? (
              <HintText>No alternatives.</HintText>
            ) : (
              <ul className="flex flex-col gap-1.5">
                {alternatives.map((ix) => (
                  <IndexCard
                    key={ix.id}
                    index={ix}
                    isActive={false}
                    latticeUnit={latticeUnit}
                    curvatureUnit={curvatureUnit}
                    data-alternative-id={ix.id}
                    onAction={() => addMember.mutate(ix.id)}
                    onHover={() => setHoveredIndex(ix.id)}
                    onLeave={() => setHoveredIndex(undefined)}
                  />
                ))}
              </ul>
            )}
          </div>

          {/* Speculative — user-built sub-minpeaks indices */}
          <div>
            <GroupHead label="Speculative" count={speculatives.length} />
            {speculatives.length === 0 ? (
              <HintText>No speculative indices yet.</HintText>
            ) : (
              <ul className="flex flex-col gap-1.5">
                {speculatives.map((ix) => {
                  const inActive = memberIds.has(ix.id);
                  return (
                    <IndexCard
                      key={ix.id}
                      index={ix}
                      isActive={inActive}
                      latticeUnit={latticeUnit}
                      curvatureUnit={curvatureUnit}
                      onAction={() =>
                        inActive
                          ? active && removeMember.mutate(ix.id)
                          : active && addMember.mutate(ix.id)
                      }
                      onDelete={() => deleteIndex.mutate(ix.id)}
                      onHover={() => setHoveredIndex(ix.id)}
                      onLeave={() => setHoveredIndex(undefined)}
                    />
                  );
                })}
              </ul>
            )}
            <button
              type="button"
              data-testid="add-speculative-button"
              className="mt-2 w-full text-xs text-fg-dim border border-dashed border-border-soft rounded-md py-1.5 hover:text-fg hover:bg-bg-hover transition-colors"
              onClick={() => openBuilder(exposureId)}
            >
              + Add speculative
            </button>
          </div>

        </div>
      </Skeleton>

      {builder && builder.exposureId === exposureId && (
        <SpeculativeBuilder exposureId={exposureId} onClose={closeBuilder} />
      )}
    </div>
  );
}
