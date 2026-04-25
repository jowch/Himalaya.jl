import { useIndices, useGroups, useAddIndexToGroup, useRemoveIndexFromGroup } from "../queries";
import { useAppState } from "../state";
import { phaseColor } from "../phases";
import { HintText } from "./ui";
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

// ── Individual index card ────────────────────────────────────────────────────

interface IndexCardProps {
  index: IndexEntry;
  isActive: boolean;
  onAction: () => void;
  onHover?: () => void;
  onLeave?: () => void;
  /** Forwarded to the li for E2E selectors */
  "data-alternative-id"?: number;
}

function IndexCard({ index, isActive, onAction, onHover, onLeave, "data-alternative-id": altId }: IndexCardProps): JSX.Element {
  const color = phaseColor(index.phase);
  const lowR2 = !isActive && index.r_squared != null && index.r_squared < R2_THRESHOLD;

  return (
    <li
      data-index-id={index.id}
      data-alternative-id={altId}
      data-active={isActive || undefined}
      data-low-r2={lowR2 || undefined}
      className={[
        "flex items-start gap-2.5 rounded-xl border px-3 py-2.5 transition-all",
        isActive ? "" : "border-transparent hover:bg-bg-hover cursor-default",
        lowR2 ? "opacity-40" : "",
      ].join(" ")}
      style={isActive ? {
        background: `${color}12`,
        borderColor: `${color}38`,
        boxShadow: `0 4px 16px ${color}18`,
      } : undefined}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
    >
      {/* Phase color dot — glows on active items */}
      <span
        className="mt-[5px] w-2 h-2 rounded-full shrink-0"
        aria-label={index.phase}
        role="img"
        style={{
          background: color,
          boxShadow: isActive ? `0 0 6px ${color}99` : undefined,
        }}
      />

      {/* Text block */}
      <div className="flex-1 min-w-0">
        {/* Line 1: phase name + lattice parameter */}
        <div className="flex items-baseline gap-2 min-w-0">
          <span className="font-semibold text-[13px] shrink-0">
            {index.phase}
          </span>
          <span className="text-[12px] text-fg-muted truncate">
            {formatLattice(index.lattice_d)}
            <span className="text-[10px] ml-0.5">nm</span>
          </span>
          {lowR2 && (
            <span className="text-[10px] text-warning uppercase tracking-wider">
              low R²
            </span>
          )}
        </div>
        {/* Line 2: score · R² as quiet metadata */}
        <div className="flex gap-3 mt-0.5 text-[11px] text-fg-dim items-center">
          <span>
            R²{" "}
            <span className={index.r_squared != null && index.r_squared >= R2_THRESHOLD
              ? "text-fg-muted" : "text-fg-dim"}>
              {formatR2(index.r_squared)}
            </span>
          </span>
          <span className="flex items-center gap-1.5">
            score{" "}
            <span className="text-fg-muted">{formatScore(index.score)}</span>
            {!isActive && index.score != null && (
              <span className="inline-block w-8 h-1 bg-bg-hover rounded overflow-hidden">
                <span
                  data-score-bar
                  className="block h-full bg-accent"
                  style={{ width: `${Math.round(index.score * 100)}%` }}
                />
              </span>
            )}
          </span>
        </div>
      </div>

      {/* Add / remove button */}
      <button
        className="w-5 h-5 shrink-0 mt-0.5 rounded-md flex items-center justify-center text-[14px] text-fg-dim hover:text-fg hover:bg-bg-hover transition-colors"
        onClick={onAction}
        aria-label={isActive ? `Remove index ${index.id}` : `Add index ${index.id}`}
      >
        {isActive ? "−" : "+"}
      </button>
    </li>
  );
}

// ── Group heading ────────────────────────────────────────────────────────────

function GroupHead({ label, count }: { label: string; count: number }): JSX.Element {
  return (
    <div className="flex items-center justify-between mb-2 px-1">
      <span className="text-[10px] uppercase tracking-widest text-fg-dim font-semibold">
        {label}
      </span>
      <span className="text-[10px] text-fg-dim">{count}</span>
    </div>
  );
}

// ── Panel ────────────────────────────────────────────────────────────────────

export interface PhasePanelProps {
  exposureId: number | undefined;
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const indicesQ = useIndices(exposureId);
  const groupsQ  = useGroups(exposureId);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const active = (groupsQ.data && activeGroup(groupsQ.data)) ?? undefined;
  const addMember    = useAddIndexToGroup(exposureId ?? 0, active?.id ?? 0);
  const removeMember = useRemoveIndexFromGroup(exposureId ?? 0, active?.id ?? 0);

  if (exposureId === undefined) {
    return (
      <div className="p-4">
        <HintText>No exposure selected.</HintText>
      </div>
    );
  }
  if (indicesQ.isPending || groupsQ.isPending) {
    return (
      <div className="p-4">
        <HintText>Loading phase assignments…</HintText>
      </div>
    );
  }

  const indices = (indicesQ.data ?? []).slice().sort(
    (a, b) => (b.score ?? 0) - (a.score ?? 0),
  );
  const memberIds = new Set(active?.members ?? []);
  const activeMembers = indices.filter((ix) =>  memberIds.has(ix.id));
  const alternatives  = indices.filter((ix) => !memberIds.has(ix.id));

  return (
    <div className="flex flex-col h-full min-h-0">

      {/* ── Sticky header ── */}
      <div className="card-header">
        <div className="flex flex-col justify-center min-w-0">
          <div className="text-[13px] font-semibold tracking-tight">Index choices</div>
          <div className="text-[11px] text-fg-dim leading-tight">
            Hover a candidate to preview peaks
          </div>
        </div>
      </div>

      {/* ── Scrollable list ── */}
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
                  data-alternative-id={ix.id}
                  onAction={() => addMember.mutate(ix.id)}
                  onHover={() => setHoveredIndex(ix.id)}
                  onLeave={() => setHoveredIndex(undefined)}
                />
              ))}
            </ul>
          )}
        </div>

      </div>
    </div>
  );
}
