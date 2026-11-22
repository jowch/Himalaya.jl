import type { MutationCache, QueryClient } from "@tanstack/react-query";
import type { OpKind, Mutator } from "./types";
import type { PersistedOpForResolution } from "./mutatorRegistry";
import { queueMutationOptions } from "./useQueueMutation";

export const STORAGE_KEY = "himalaya-ui:queue";
export const SCHEMA_VERSION = 6;

interface PersistedOp {
  schemaVersion: number;
  kind: OpKind;
  clientOpId: string;
  payload: unknown;
}

export type MutatorResolver = (
  op: PersistedOpForResolution,
) => Mutator<any, any, any> | undefined;
export function attachPersistence(mc: MutationCache): () => void {
  const dump = () => mirrorToSessionStorage(mc);
  dump();
  return mc.subscribe(() => dump());
}
function mirrorToSessionStorage(mc: MutationCache): void {
  const pending = mc.getAll().filter((m) => m.state.status === "pending");
  const ops: PersistedOp[] = [];
  for (const m of pending) {
    const vars = m.state.variables as
      | (Record<string, unknown> & { kind?: OpKind; clientOpId?: string })
      | undefined;
    if (!vars?.kind || !vars.clientOpId) continue;
    ops.push({
      schemaVersion: SCHEMA_VERSION,
      kind: vars.kind,
      clientOpId: vars.clientOpId,
      payload: vars,
    });
  }
  try {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify(ops));
  } catch {
  }
}
export interface RehydrateResult {
  dropped: number;
  replayed: number;
  failed: number;
  lost: number;
}
export async function rehydrate(
  qc: QueryClient,
  mutators: Map<OpKind, Mutator<any, any, any>> | MutatorResolver,
): Promise<RehydrateResult> {
  const raw = sessionStorage.getItem(STORAGE_KEY);
  if (!raw) return { dropped: 0, replayed: 0, failed: 0, lost: 0 };
  sessionStorage.removeItem(STORAGE_KEY);
  let ops: PersistedOp[] = [];
  try {
    const parsed = JSON.parse(raw) as unknown;
    if (Array.isArray(parsed)) ops = parsed as PersistedOp[];
  } catch {
    return { dropped: 0, replayed: 0, failed: 0, lost: 0 };
  }
  const lookup = (op: PersistedOp): Mutator<any, any, any> | undefined =>
    typeof mutators === "function"
      ? mutators({ kind: op.kind, payload: op.payload })
      : mutators.get(op.kind);
  let dropped = 0;
  let replayed = 0;
  let failed = 0;
  let lost = 0;
  const fires: Promise<unknown>[] = [];
  const mc = qc.getMutationCache();
  const inFlight = new Set(
    mc.getAll()
      .filter((m) => m.state.status === "pending")
      .map((m) => (m.state.variables as { clientOpId?: string } | undefined)?.clientOpId),
  );
  for (const op of ops) {
    try {
      if (op.schemaVersion !== SCHEMA_VERSION) {
        dropped++;
        continue;
      }
      if (inFlight.has(op.clientOpId)) continue;
      const mutator = lookup(op);
      if (!mutator) {
        dropped++;
        continue;
      }
      const options = queueMutationOptions(mutator, qc);
      let unappliable = false;
      const onMutate = options.onMutate!;
      const mutation = mc.build(qc, {
        ...options,
        onMutate: (payload, ctx) => {
          try {
            return onMutate(payload, ctx);
          } catch (err) {
            unappliable = true;
            throw err;
          }
        },
      });
      fires.push(
        mutation
          .execute(op.payload as Record<string, unknown> & { kind: string; clientOpId: string })
          .then(() => { replayed++; })
          .catch(() => {
            if (unappliable) lost++;
            else failed++;
          }),
      );
    } catch {
      lost++;
    }
  }
  await Promise.all(fires);
  return { dropped, replayed, failed, lost };
}