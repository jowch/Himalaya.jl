import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { usePostSampleMessage, queryKeys } from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { SampleMessage } from "../src/api";

const SAMPLE_ID = 10;

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — sample message mutations (queue-driven)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("appends an optimistic placeholder with a negative id", () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.messages(SAMPLE_ID), [] as SampleMessage[]);
    vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
    const { result } = renderHook(
      () => usePostSampleMessage(SAMPLE_ID), { wrapper },
    );
    act(() => { result.current.mutate("hello world"); });
    const list = client.getQueryData(queryKeys.messages(SAMPLE_ID)) as SampleMessage[];
    expect(list).toHaveLength(1);
    expect(list[0].id).toBeLessThan(0);
    expect(list[0].body).toBe("hello world");
    expect(list[0].sample_id).toBe(SAMPLE_ID);
  });

  it("replaces the placeholder with the server message on success", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.messages(SAMPLE_ID), [] as SampleMessage[]);
    mockOnce(201, {
      id: 42,
      sample_id: SAMPLE_ID,
      author_id: 1,
      author: "alice",
      body: "hello world",
      created_at: "2026-05-02T12:00:00Z",
    });
    const { result } = renderHook(
      () => usePostSampleMessage(SAMPLE_ID), { wrapper },
    );
    act(() => { result.current.mutate("hello world"); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.messages(SAMPLE_ID)) as SampleMessage[];
      expect(list).toHaveLength(1);
      expect(list[0].id).toBe(42);
      expect(list[0].author).toBe("alice");
    });
  });

  it("appends to existing messages, preserving prior history", async () => {
    const { client, wrapper } = withClient();
    const prior: SampleMessage = {
      id: 1, sample_id: SAMPLE_ID, author_id: 2, author: "bob",
      body: "earlier", created_at: "2026-05-02T10:00:00Z",
    };
    client.setQueryData(queryKeys.messages(SAMPLE_ID), [prior]);
    mockOnce(201, {
      id: 43, sample_id: SAMPLE_ID, author_id: 1, author: "alice",
      body: "later", created_at: "2026-05-02T12:00:00Z",
    });
    const { result } = renderHook(
      () => usePostSampleMessage(SAMPLE_ID), { wrapper },
    );
    act(() => { result.current.mutate("later"); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.messages(SAMPLE_ID)) as SampleMessage[];
      expect(list).toHaveLength(2);
      expect(list[0].id).toBe(1);
      expect(list[1].id).toBe(43);
    });
  });

  it("rolls back the placeholder on 4xx", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.messages(SAMPLE_ID), [] as SampleMessage[]);
    mockOnce(400, { error: "empty body" });
    const { result } = renderHook(
      () => usePostSampleMessage(SAMPLE_ID), { wrapper },
    );
    act(() => { result.current.mutate("nope"); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.messages(SAMPLE_ID)) as SampleMessage[];
      expect(list).toHaveLength(0);
    });
  });
});
