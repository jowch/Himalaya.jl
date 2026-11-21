import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useSamples,
} from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

describe("useSamples", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("asks for nothing until an experiment is picked", async () => {
    const fetchSpy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSamples(0), { wrapper });
    await waitFor(() => expect(result.current.fetchStatus).toBe("idle"));
    expect(fetchSpy).not.toHaveBeenCalled();
  });
});
