import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { addLabelMutator, removeLabelMutator, restoreLabelSuggestionMutator } from "../../src/lib/queue/mutators/labels";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { SampleLabel } from "../../src/api";

function frame(kind: string, entity_id: number, payload: Record<string, unknown>): SseEvent {
  return { kind, entity_id, entity_type: "sample", payload } as unknown as SseEvent;
}

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "add_label", payload: {}, experimentId: 7 } as const;

function row(id: number, sample_id: number, label: string, state: SampleLabel["state"] = "confirmed"): SampleLabel {
  return { id, sample_id, label, state, source: "user" };
}

describe("label mutations", () => {
  it("a remote label event refetches that experiment's labels and nothing else", () => {
    for (const kind of ["label_added", "label_removed", "label_suggestion_restored"]) {
      const qc = new QueryClient();
      const spy = vi.spyOn(qc, "invalidateQueries");
      applyRemoteToCache(frame(kind, 42, { label: "37C", experiment_id: 7 }), qc);
      const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
      expect(keys).toEqual([JSON.stringify(queryKeys.labels(7))]);
    }
  });

  it("routes both the outbound op and the echo back to the label mutators", () => {
    expect(resolveMutator({ kind: "add_label", payload: {} })).toBe(addLabelMutator);
    expect(resolveMutator({ kind: "remove_label", payload: {} })).toBe(removeLabelMutator);
    expect(resolveMutatorForEvent("label_added")).toBe(addLabelMutator);
    expect(resolveMutatorForEvent("label_removed")).toBe(removeLabelMutator);
    expect(resolveMutator({ kind: "restore_label_suggestion", payload: {} })).toBe(restoreLabelSuggestionMutator);
    expect(resolveMutatorForEvent("label_suggestion_restored")).toBe(restoreLabelSuggestionMutator);
  });

  it("adding shows the tag at once and rolls back cleanly", () => {
    const qc = new QueryClient();
    const before = [row(1, 3, "buffer", "suggested")];
    qc.setQueryData(queryKeys.labels(7), before);
    const ctx = addLabelMutator.onMutate({ ...scope, sampleId: 3, label: "Buffer" }, qc);
    expect(qc.getQueryData<SampleLabel[]>(queryKeys.labels(7))).toEqual([row(1, 3, "buffer", "confirmed")]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.labels(7))).toEqual(before);
  });

  it("removing hides the tag at once and rolls back cleanly", () => {
    const qc = new QueryClient();
    const before = [row(1, 3, "retake"), row(2, 4, "retake")];
    qc.setQueryData(queryKeys.labels(7), before);
    const ctx = removeLabelMutator.onMutate({ ...scope, kind: "remove_label", sampleId: 3, label: "retake" }, qc);
    expect(qc.getQueryData<SampleLabel[]>(queryKeys.labels(7))).toEqual([row(2, 4, "retake")]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.labels(7))).toEqual(before);
  });

  it("restoring a suggestion shows it as suggested again, even after a dismissal hid it, and rolls back cleanly", () => {
    const qc = new QueryClient();
    const before = [row(1, 3, "buffer"), row(2, 4, "retake")];
    qc.setQueryData(queryKeys.labels(7), before);
    const ctx = restoreLabelSuggestionMutator.onMutate({ ...scope, kind: "restore_label_suggestion", sampleId: 3, label: "Buffer" }, qc);
    expect(qc.getQueryData<SampleLabel[]>(queryKeys.labels(7))).toEqual([{ ...row(1, 3, "buffer", "suggested"), source: "ingest" }, row(2, 4, "retake")]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.labels(7))).toEqual(before);
    restoreLabelSuggestionMutator.onMutate({ ...scope, kind: "restore_label_suggestion", sampleId: 4, label: "calibrant" }, qc);
    expect(qc.getQueryData<SampleLabel[]>(queryKeys.labels(7))!.at(-1)).toEqual({ id: -1, sample_id: 4, label: "calibrant", state: "suggested", source: "ingest" });
  });

  it("takes the server's list when the response is one, and refetches when it is a partial from the live stream", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    const server = [row(9, 3, "37C")];
    addLabelMutator.onSuccess({ ...scope, sampleId: 3, label: "37C" }, server, qc);
    expect(qc.getQueryData(queryKeys.labels(7))).toEqual(server);
    addLabelMutator.onSuccess({ ...scope, sampleId: 3, label: "37C" }, { id: 1 } as unknown as SampleLabel[], qc);
    expect(qc.getQueryData(queryKeys.labels(7))).toEqual(server);
    expect(spy).toHaveBeenCalledWith({ queryKey: queryKeys.labels(7) });
  });

  it("sends the label to the right route with the auth headers", async () => {
    const fetchMock = vi.fn().mockImplementation(async () => new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }));
    vi.stubGlobal("fetch", fetchMock);
    await addLabelMutator.request({ ...scope, sampleId: 3, label: "37C" }, new AbortController().signal);
    await removeLabelMutator.request({ ...scope, kind: "remove_label", sampleId: 3, label: "empty capillary" }, new AbortController().signal);
    await restoreLabelSuggestionMutator.request({ ...scope, kind: "restore_label_suggestion", sampleId: 3, label: "Buffer" }, new AbortController().signal);
    vi.unstubAllGlobals();
    const [addUrl, addInit] = fetchMock.mock.calls[0]!;
    expect(String(addUrl)).toMatch(/\/api\/samples\/3\/labels$/);
    expect(addInit.method).toBe("POST");
    expect(JSON.parse(addInit.body)).toEqual({ label: "37C" });
    expect(new Headers(addInit.headers).get("X-Username")).toBe("jc");
    expect(new Headers(addInit.headers).get("X-Client-Op-Id")).toBe("op");
    const [delUrl, delInit] = fetchMock.mock.calls[1]!;
    expect(String(delUrl)).toMatch(/\/api\/samples\/3\/labels\?label=empty%20capillary$/);
    expect(delInit.method).toBe("DELETE");
    const [restoreUrl, restoreInit] = fetchMock.mock.calls[2]!;
    expect(String(restoreUrl)).toMatch(/\/api\/samples\/3\/labels\/suggested$/);
    expect(restoreInit.method).toBe("POST");
    expect(JSON.parse(restoreInit.body)).toEqual({ label: "Buffer" });
    expect(new Headers(restoreInit.headers).get("X-Client-Op-Id")).toBe("op");
  });
});
