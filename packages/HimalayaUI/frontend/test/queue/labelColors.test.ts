import { afterEach, describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { setLabelColorMutator } from "../../src/lib/queue/mutators/labelColors";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { LabelColorEntry } from "../../src/api";

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "set_tag_color" } as const;

describe("label colour mutations", () => {
  afterEach(() => vi.unstubAllGlobals());

  it("a remote colour change refetches the shared colour map and nothing else", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    const frame = { kind: "tag_color_set", entity_id: 0, entity_type: "tag", payload: { label: "37C", color: "teal" } } as unknown as SseEvent;
    applyRemoteToCache(frame, qc);
    expect(spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey))).toEqual([JSON.stringify(queryKeys.labelColors())]);
  });

  it("routes the outbound op and its echo to the colour mutator", () => {
    expect(resolveMutator({ kind: "set_tag_color", payload: {} })).toBe(setLabelColorMutator);
    expect(resolveMutatorForEvent("tag_color_set")).toBe(setLabelColorMutator);
  });

  it("recolours at once, matching the name without case, and rolls back cleanly", () => {
    const qc = new QueryClient();
    const before: LabelColorEntry[] = [{ label: "37C", color: "teal" }, { label: "Buffer", color: "rose" }];
    qc.setQueryData(queryKeys.labelColors(), before);
    const ctx = setLabelColorMutator.onMutate({ ...scope, label: "37c", color: "stone" }, qc);
    expect(qc.getQueryData(queryKeys.labelColors())).toEqual([{ label: "Buffer", color: "rose" }, { label: "37c", color: "stone" }]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.labelColors())).toEqual(before);
    setLabelColorMutator.onMutate({ ...scope, label: "buffer", color: null }, qc);
    expect(qc.getQueryData(queryKeys.labelColors())).toEqual([{ label: "37C", color: "teal" }]);
  });

  it("takes the server's map when the response is one", () => {
    const qc = new QueryClient();
    const server: LabelColorEntry[] = [{ label: "37C", color: "teal" }];
    setLabelColorMutator.onSuccess({ ...scope, label: "37c", color: "teal" }, server, qc);
    expect(qc.getQueryData(queryKeys.labelColors())).toEqual(server);
  });

  it("sends a PUT with the label, the colour and the auth headers", async () => {
    const fetchMock = vi.fn().mockImplementation(async () => new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }));
    vi.stubGlobal("fetch", fetchMock);
    await setLabelColorMutator.request({ ...scope, label: "37C", color: null }, new AbortController().signal);
    const [url, init] = fetchMock.mock.calls[0]!;
    expect(String(url)).toMatch(/\/api\/label-colors$/);
    expect(init.method).toBe("PUT");
    expect(JSON.parse(init.body)).toEqual({ label: "37C", color: null });
    expect(new Headers(init.headers).get("X-Username")).toBe("jc");
    expect(new Headers(init.headers).get("X-Client-Op-Id")).toBe("op");
  });
});
