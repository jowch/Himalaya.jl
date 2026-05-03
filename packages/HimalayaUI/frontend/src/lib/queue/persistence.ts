import type { MutationCache, QueryClient } from "@tanstack/react-query";
import type { OpKind, Mutator } from "./types";
import type { PersistedOpForResolution } from "./mutatorRegistry";

export const STORAGE_KEY = "himalaya-ui:queue";
export const SCHEMA_VERSION = 1;

interface PersistedOp {
  schemaVersion: number;
  kind: OpKind;
  clientOpId: string;
  payload: unknown;
}

/**
 * A function that, given a persisted op, returns the Mutator that should
 * replay it (or undefined if the op should be dropped). Used for kinds whose
 * dispatch depends on the payload shape (e.g. `add_tag` is sample-scoped or
 * exposure-scoped depending on whether the payload carries `experimentId`).
 */
export type MutatorResolver = (
  op: PersistedOpForResolution,
) => Mutator<any, any> | undefined;

/**
 * Subscribe to MutationCache events; mirror the current pending set to
 * sessionStorage on every add/update/remove. Returns an unsubscribe.
 *
 * Mirror happens unconditionally on every cache event — cheap because the
 * pending set is bounded (one in-flight op per user gesture, typically <10
 * concurrent). The full set is rewritten so transitions from pending →
 * success/error correctly drop the entry.
 */
export function attachPersistence(mc: MutationCache): () => void {
  const dump = () => mirrorToSessionStorage(mc);
  // Initial mirror so existing pending state at attach time lands in storage.
  dump();
  return mc.subscribe(() => dump());
}

function mirrorToSessionStorage(mc: MutationCache): void {
  const pending = mc.getAll().filter((m) => m.state.status === "pending");
  const ops: PersistedOp[] = [];
  for (const m of pending) {
    // useQueueMutation flat-spreads scope + input into `variables`; there is
    // no nested `payload` key. Persist the full variables object so rehydrate
    // can pass it back into mutator.onMutate / .request unchanged. Mutators
    // already destructure scope+input from this flat shape via their `flat()`
    // helper, so the rehydrate path matches the live-mutate path.
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
    // Quota exceeded or sessionStorage unavailable. Best-effort: skip.
  }
}

export interface RehydrateResult {
  dropped: number;
  replayed: number;
  failed: number;
}

/**
 * On tab reload, replay each persisted op through its registered mutator.
 * Server-side request-level idempotency (X-Client-Op-Id; M0.3) returns the
 * cached response on retry, so duplicate execution is impossible — this
 * call settles the local cache to match the server's view.
 *
 * Returns counts; the caller surfaces a toast for `dropped > 0` so the user
 * knows some edits from the prior session couldn't be replayed (schema-
 * version drift after a deploy, or unknown kind after a feature removal).
 *
 * The storage key is cleared after rehydrate so a subsequent reload starts
 * from a clean slate.
 */
export async function rehydrate(
  qc: QueryClient,
  mutators: Map<OpKind, Mutator<any, any>> | MutatorResolver,
): Promise<RehydrateResult> {
  const raw = sessionStorage.getItem(STORAGE_KEY);
  if (!raw) return { dropped: 0, replayed: 0, failed: 0 };

  let ops: PersistedOp[] = [];
  try {
    const parsed = JSON.parse(raw) as unknown;
    if (Array.isArray(parsed)) ops = parsed as PersistedOp[];
  } catch {
    sessionStorage.removeItem(STORAGE_KEY);
    return { dropped: 0, replayed: 0, failed: 0 };
  }

  const lookup = (op: PersistedOp): Mutator<any, any> | undefined =>
    typeof mutators === "function"
      ? mutators({ kind: op.kind, payload: op.payload })
      : mutators.get(op.kind);

  let dropped = 0;
  let replayed = 0;
  let failed = 0;
  const fires: Promise<unknown>[] = [];

  for (const op of ops) {
    if (op.schemaVersion !== SCHEMA_VERSION) {
      dropped++;
      continue;
    }
    const mutator = lookup(op);
    if (!mutator) {
      dropped++;
      continue;
    }
    // Re-run optimistic effect against fresh cache.
    mutator.onMutate(op.payload, qc);
    // Re-fire the request. AbortController is not wired through to rehydrate
    // because the original mutate() call's signal is gone; rehydrate-fired
    // requests run to completion and rely on idempotency for safety.
    const ctrl = new AbortController();
    fires.push(
      mutator
        .request(op.payload, ctrl.signal)
        .then((response) => {
          mutator.onSuccess(op.payload, response, qc);
          replayed++;
        })
        .catch(() => {
          // Retried request failed. The HTTP retry semantics in M1.5 will
          // also kick in for any in-flight session; here we count it so
          // the caller can surface a "couldn't replay" toast accurately.
          failed++;
        }),
    );
  }

  // Wait for all replays to settle so callers can act on the result.
  await Promise.all(fires);
  sessionStorage.removeItem(STORAGE_KEY);
  return { dropped, replayed, failed };
}
