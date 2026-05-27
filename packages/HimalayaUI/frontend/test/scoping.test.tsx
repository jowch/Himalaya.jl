import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, render, screen, fireEvent, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useCorpusSampleTags, useCorpusPickerSamples } from "../src/queries";
import { ScopingConfirmModal } from "../src/components/ScopingConfirmModal";

function wrap(client = makeClient()) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
}

describe("scoping read hooks", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("useCorpusSampleTags fetches GET /api/sample-tags", async () => {
    const seen: string[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      seen.push(url);
      return Promise.resolve(
        new Response(JSON.stringify([{ key: "ratio", value: "1:1" }]), {
          status: 200, headers: { "Content-Type": "application/json" },
        }),
      );
    });
    const { result } = renderHook(() => useCorpusSampleTags(), { wrapper: wrap() });
    await waitFor(() => expect(result.current.data).toEqual([{ key: "ratio", value: "1:1" }]));
    expect(seen.some((u) => u.endsWith("/api/sample-tags"))).toBe(true);
  });

  it("useCorpusPickerSamples fetches GET /api/picker-samples", async () => {
    const seen: string[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      seen.push(url);
      return Promise.resolve(
        new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }),
      );
    });
    const { result } = renderHook(() => useCorpusPickerSamples(), { wrapper: wrap() });
    await waitFor(() => expect(result.current.data).toEqual([]));
    expect(seen.some((u) => u.endsWith("/api/picker-samples"))).toBe(true);
  });
});

describe("ScopingConfirmModal", () => {
  it("summarizes the write and fires onConfirm", () => {
    const onConfirm = vi.fn();
    const onClose = vi.fn();
    render(
      <ScopingConfirmModal open orderingKey="ratio" count={3}
        onConfirm={onConfirm} onClose={onClose} />,
    );
    const dialog = screen.getByTestId("scoping-confirm-modal");
    expect(dialog).toHaveAttribute("role", "dialog");
    expect(dialog).toHaveTextContent("ratio");
    expect(dialog).toHaveTextContent("3");
    fireEvent.click(screen.getByTestId("scoping-confirm-build"));
    expect(onConfirm).toHaveBeenCalledTimes(1);
  });

  it("does not render when closed", () => {
    render(<ScopingConfirmModal open={false} orderingKey="ratio" count={3}
      onConfirm={() => {}} onClose={() => {}} />);
    expect(screen.queryByTestId("scoping-confirm-modal")).toBeNull();
  });

  it("calls onClose on Escape", () => {
    const onClose = vi.fn();
    render(<ScopingConfirmModal open orderingKey="ratio" count={3}
      onConfirm={() => {}} onClose={onClose} />);
    fireEvent.keyDown(screen.getByTestId("scoping-confirm-modal"), { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);
  });
});
