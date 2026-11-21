import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { act, renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { patchExperimentMutator } from "../../src/lib/queue/mutators/experiment";
import { useQueueMutation } from "../../src/lib/queue/useQueueMutation";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { setToastImpl } from "../../src/lib/toast";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { Experiment } from "../../src/api";

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "patch_experiment", payload: {}, experimentId: 3 } as const;

function experiment(over: Partial<Experiment> = {}): Experiment {
  return { id: 3, name: "BL-19", description: null, ...over } as Experiment;
}

function seeded(): QueryClient {
  const qc = new QueryClient({ defaultOptions: { mutations: { retry: false } } });
  qc.setQueryData(queryKeys.experiment(3), experiment());
  qc.setQueryData(queryKeys.experiments, [experiment(), experiment({ id: 4, name: "BL-20" })]);
  return qc;
}

describe("experiment settings mutations", () => {
  beforeEach(() => pendingDeferreds.clear());
  afterEach(() => vi.unstubAllGlobals());

  it("a remote rename updates the open experiment and the experiments list", () => {
    const qc = seeded();
    const frame = { kind: "experiment_updated", entity_id: 3, entity_type: "experiment", payload: { name: "BL-19b", description: null } } as unknown as SseEvent;
    applyRemoteToCache(frame, qc);
    expect(qc.getQueryData<Experiment>(queryKeys.experiment(3))?.name).toBe("BL-19b");
    expect(qc.getQueryData<Experiment[]>(queryKeys.experiments)?.map((e) => e.name)).toEqual(["BL-19b", "BL-20"]);
  });

  it("routes the outbound op and its echo to the experiment mutator", () => {
    expect(resolveMutator({ kind: "patch_experiment", payload: {} })).toBe(patchExperimentMutator);
    expect(resolveMutatorForEvent("experiment_updated")).toBe(patchExperimentMutator);
  });

  it("renames at once in both caches and rolls back cleanly", () => {
    const qc = seeded();
    const ctx = patchExperimentMutator.onMutate({ ...scope, name: "New" }, qc);
    expect(qc.getQueryData<Experiment>(queryKeys.experiment(3))?.name).toBe("New");
    expect(qc.getQueryData<Experiment[]>(queryKeys.experiments)?.[0]?.name).toBe("New");
    ctx.restore();
    expect(qc.getQueryData<Experiment>(queryKeys.experiment(3))?.name).toBe("BL-19");
    expect(qc.getQueryData<Experiment[]>(queryKeys.experiments)?.[0]?.name).toBe("BL-19");
  });

  it("takes only the fields the response carries, so an SSE-synthesized diff keeps the rest", () => {
    const qc = seeded();
    qc.setQueryData(queryKeys.experiment(3), experiment({ description: "kept" }));
    patchExperimentMutator.onSuccess({ ...scope, name: "Server" }, { name: "Server" } as Experiment, qc);
    expect(qc.getQueryData<Experiment>(queryKeys.experiment(3))).toMatchObject({ name: "Server", description: "kept" });
  });

  it("sends a PATCH with only the patched fields and the auth headers", async () => {
    const fetchMock = vi.fn().mockImplementation(async () => new Response("{}", { status: 200, headers: { "Content-Type": "application/json" } }));
    vi.stubGlobal("fetch", fetchMock);
    await patchExperimentMutator.request({ ...scope, description: null }, new AbortController().signal);
    const [url, init] = fetchMock.mock.calls[0]!;
    expect(String(url)).toMatch(/\/api\/experiments\/3$/);
    expect(init.method).toBe("PATCH");
    expect(JSON.parse(init.body)).toEqual({ description: null });
    expect(new Headers(init.headers).get("X-Client-Op-Id")).toBe("op");
  });

  it("a rejected save puts the old name back and says so", async () => {
    const status = 400;
    const toasts: [string, string][] = [];
    setToastImpl((msg, kind) => toasts.push([msg, kind]));
    vi.stubGlobal("fetch", vi.fn().mockImplementation(async () => new Response(JSON.stringify({ error: "nope" }), { status, headers: { "Content-Type": "application/json" } })));
    const qc = seeded();
    const wrapper = ({ children }: { children: ReactNode }) => <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
    const { result } = renderHook(() => useQueueMutation(patchExperimentMutator, { experimentId: 3, username: "jc", clientId: "c" }), { wrapper });
    act(() => result.current.mutate({ name: "Broken" }));
    await waitFor(() => expect(toasts).toHaveLength(1));
    expect(toasts[0]![0].startsWith("Couldn't save the experiment: ")).toBe(true);
    expect(toasts[0]![1]).toBe("error");
    expect(qc.getQueryData<Experiment>(queryKeys.experiment(3))?.name).toBe("BL-19");
  });
});
