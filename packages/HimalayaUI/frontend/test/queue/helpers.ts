import type { Mutation, MutationStatus } from "@tanstack/react-query";
import type { SseEvent } from "../../src/lib/queue/types";

interface FakeMutationOpts {
  status?: MutationStatus;
  context?: Record<string, unknown>;
  variables?: Record<string, unknown>;
  onMutate?: (vars: unknown) => void;
  mutationKey?: readonly unknown[];
}

/**
 * Build a Mutation-shaped object exposing the fields the replay coordinator
 * reads. Not a full TanStack Mutation — just enough to satisfy the type
 * surface that handleRemoteEvent touches.
 */
export function makeFakeMutation(opts: FakeMutationOpts = {}): Mutation {
  const state = {
    status: opts.status ?? "pending",
    context: opts.context ?? {},
    variables: opts.variables ?? {},
    data: undefined,
    error: null,
    failureCount: 0,
    failureReason: null,
    isPaused: false,
    submittedAt: Date.now(),
  };
  const options = {
    mutationKey: opts.mutationKey,
    onMutate: opts.onMutate,
  };
  return { state, options } as unknown as Mutation;
}

let foreignCounter = 0;
export function remoteForeignEvent(overrides: Partial<SseEvent> = {}): SseEvent {
  foreignCounter++;
  return {
    id: foreignCounter,
    kind: "peak_added",
    entity_type: "exposure",
    entity_id: 42,
    client_op_id: `foreign-${foreignCounter}`,
    payload: { q: 1.5 },
    ...overrides,
  };
}
