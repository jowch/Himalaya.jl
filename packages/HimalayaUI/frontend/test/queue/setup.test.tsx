import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { patchSetupMutator } from "../../src/lib/queue/mutators/setup";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { Experiment, Setup } from "../../src/api";

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "patch_setup", payload: {}, experimentId: 3 } as const;

function setup(id: number, folder: string, over: Partial<Setup> = {}): Setup {
  return { id, experiment_id: 3, folder, display_name: null, name: folder, included: true, ...over } as Setup;
}

function seeded(): QueryClient {
  const qc = new QueryClient({ defaultOptions: { mutations: { retry: false } } });
  const exp = { id: 3, name: "BL", setups: [setup(7, "setup1"), setup(8, "setup2")] } as Experiment;
  qc.setQueryData(queryKeys.experiment(3), exp);
  qc.setQueryData(queryKeys.experiments, [exp, { id: 4, name: "Other", setups: [setup(9, "setup1")] } as Experiment]);
  return qc;
}

const setupsOf = (qc: QueryClient, id = 3): Setup[] => qc.getQueryData<Experiment>(queryKeys.experiment(id))?.setups ?? [];

describe("setup mutations", () => {
  beforeEach(() => pendingDeferreds.clear());
  afterEach(() => vi.unstubAllGlobals());

  it("routes the outbound op and its echo to the setup mutator", () => {
    expect(resolveMutator({ kind: "patch_setup", payload: {} })).toBe(patchSetupMutator);
    expect(resolveMutatorForEvent("setup_updated")).toBe(patchSetupMutator);
  });

  it("a remote hide or rename updates that experiment's setup in both caches", () => {
    const qc = seeded();
    applyRemoteToCache({ kind: "setup_updated", entity_id: 8, entity_type: "setup",
      payload: { experiment_id: 3, included: false } } as unknown as SseEvent, qc);
    applyRemoteToCache({ kind: "setup_updated", entity_id: 7, entity_type: "setup",
      payload: { experiment_id: 3, name: "Long pipe", display_name: "Long pipe" } } as unknown as SseEvent, qc);
    expect(setupsOf(qc).map((s) => [s.name, s.display_name, s.included])).toEqual([["Long pipe", "Long pipe", true], ["setup2", null, false]]);
    const list = qc.getQueryData<Experiment[]>(queryKeys.experiments)!;
    expect(list[0]!.setups!.map((s) => s.included)).toEqual([true, false]);
    expect(list[1]!.setups![0]!.name).toBe("setup1");
  });

  it("renames at once, falls back to the folder name when cleared, and rolls back", () => {
    const qc = seeded();
    const ctx = patchSetupMutator.onMutate({ ...scope, setupId: 7, name: " Long pipe " }, qc);
    expect(setupsOf(qc)[0]).toMatchObject({ name: "Long pipe", display_name: "Long pipe" });
    ctx.restore();
    expect(setupsOf(qc)[0]).toMatchObject({ name: "setup1", display_name: null });
    patchSetupMutator.onMutate({ ...scope, setupId: 7, name: null }, qc);
    expect(setupsOf(qc)[0]).toMatchObject({ name: "setup1", display_name: null });
  });

  it("sends a PATCH with only the patched fields and the auth headers", async () => {
    const fetchMock = vi.fn().mockImplementation(async () => new Response("{}", { status: 200, headers: { "Content-Type": "application/json" } }));
    vi.stubGlobal("fetch", fetchMock);
    await patchSetupMutator.request({ ...scope, setupId: 8, included: true }, new AbortController().signal);
    const [url, init] = fetchMock.mock.calls[0]!;
    expect(String(url)).toMatch(/\/api\/setups\/8$/);
    expect(init.method).toBe("PATCH");
    expect(JSON.parse(init.body)).toEqual({ included: true });
    expect(new Headers(init.headers).get("X-Client-Op-Id")).toBe("op");
  });
});
