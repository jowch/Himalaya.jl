import { useIndices, useGroups, useAddIndexToGroup, useRemoveIndexFromGroup } from "../queries";
import { useAppState } from "../state";
import { phaseColor } from "../phases";
import type { GroupEntry } from "../api";

export interface PhasePanelProps {
  exposureId: number | undefined;
}

function activeGroup(groups: GroupEntry[]): GroupEntry | undefined {
  return groups.find((g) => g.active);
}

function formatLattice(d: number | null): string {
  return d != null ? d.toFixed(2) : "—";
}

function formatR2(r: number | null): string {
  return r != null ? r.toFixed(3) : "—";
}

export function PhasePanel({ exposureId }: PhasePanelProps): JSX.Element {
  const indicesQ = useIndices(exposureId);
  const groupsQ  = useGroups(exposureId);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const active = (groupsQ.data && activeGroup(groupsQ.data)) ?? undefined;
  const addMember    = useAddIndexToGroup(exposureId ?? 0, active?.id ?? 0);
  const removeMember = useRemoveIndexFromGroup(exposureId ?? 0, active?.id ?? 0);

  if (exposureId === undefined) {
    return <p className="text-fg-muted italic">No exposure selected.</p>;
  }
  if (indicesQ.isPending || groupsQ.isPending) {
    return <p className="text-fg-muted">Loading phase assignments…</p>;
  }

  const indices = (indicesQ.data ?? []).slice().sort(
    (a, b) => (b.score ?? 0) - (a.score ?? 0),
  );
  const memberIds = new Set(active?.members ?? []);
  const activeMembers = indices.filter((ix) => memberIds.has(ix.id));
  const alternatives  = indices.filter((ix) => !memberIds.has(ix.id));

  return (
    <div className="flex flex-col gap-3">
      <section>
        <h3 className="text-fg-muted text-[12px] uppercase tracking-wide mb-1">Active group</h3>
        {activeMembers.length === 0 ? (
          <p className="text-fg-muted italic text-[13px]">No indices in the active group.</p>
        ) : (
          <ul className="flex flex-col gap-1">
            {activeMembers.map((ix) => (
              <li
                key={ix.id}
                data-index-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md bg-bg-elevated"
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-error focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Remove index ${ix.id}`}
                  onClick={() => { if (active) removeMember.mutate(ix.id); }}
                >
                  −
                </button>
              </li>
            ))}
          </ul>
        )}
      </section>

      <section>
        <h3 className="text-fg-muted text-[12px] uppercase tracking-wide mb-1">
          Alternatives ({alternatives.length})
        </h3>
        {alternatives.length === 0 ? (
          <p className="text-fg-muted italic text-[13px]">No alternatives.</p>
        ) : (
          <ul className="flex flex-col gap-1">
            {alternatives.map((ix) => (
              <li
                key={ix.id}
                data-alternative-id={ix.id}
                className="flex items-center gap-2 px-2 py-1 rounded-md hover:bg-bg-hover cursor-default"
                onMouseEnter={() => setHoveredIndex(ix.id)}
                onMouseLeave={() => setHoveredIndex(undefined)}
              >
                <span
                  className="w-2 h-2 rounded-full"
                  style={{ background: phaseColor(ix.phase) }}
                  aria-hidden
                />
                <span className="font-medium">{ix.phase}</span>
                <span className="text-fg-muted text-[13px]">
                  a={formatLattice(ix.lattice_d)} nm · R²={formatR2(ix.r_squared)}
                </span>
                <button
                  className="ml-auto text-fg-muted hover:text-success focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-md px-1"
                  aria-label={`Add index ${ix.id}`}
                  onClick={() => { if (active) addMember.mutate(ix.id); }}
                  disabled={active === undefined}
                >
                  +
                </button>
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  );
}
