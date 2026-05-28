import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, render, screen, fireEvent, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { makeClient } from "./test-utils";
import { useCorpusSampleTags, useCorpusPickerSamples } from "../src/queries";
import { ScopingConfirmModal } from "../src/components/ScopingConfirmModal";
import { SeriesScopingPage } from "../src/pages/SeriesScopingPage";

function jsonRes(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status, headers: { "Content-Type": "application/json" },
  });
}

function pickerRow(id: number, name: string, ratio?: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: ratio === undefined ? [] : [{ id, key: "ratio", value: ratio, source: "manual" }],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

function renderScoping() {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={["/series/new"]}>
        <Routes>
          <Route path="/series/new" element={<SeriesScopingPage />} />
          {/* D1: post-confirm lands on the folio (exists since I3.3). */}
          <Route path="/series" element={<div data-testid="folio-stub" />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

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

  it("the confirm-build button uses Print ink tokens, not the accent (S-B)", () => {
    render(<ScopingConfirmModal open orderingKey="ratio" count={2}
      onConfirm={() => {}} onClose={() => {}} />);
    const btn = screen.getByTestId("scoping-confirm-build");
    expect(btn.className).toContain("bg-ink");
    expect(btn.className).toContain("text-paper");
    expect(btn.className).not.toContain("bg-accent");
  });
});

describe("SeriesScopingPage", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("proposes an ordering variable from corpus tags and lists member rows", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags"))
        return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }, { key: "ratio", value: "2:1" }]));
      if (url.endsWith("/api/picker-samples"))
        return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B", "2:1")]));
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    await waitFor(() => expect(screen.getByTestId("scoping-page")).toBeInTheDocument());
    // Once the queries resolve the proposal surfaces the ordering key in the
    // "Ordered by" field (the bare "Ordering variable: —" header is retired).
    await waitFor(() =>
      expect(screen.getByTestId("scoping-order-field")).toHaveTextContent("ratio"));
    await waitFor(() =>
      expect(screen.getByTestId("scoping-row-10")).toHaveTextContent("A"));
    // Confirm-not-fill-out (S-E): the value renders as ink text, not an input.
    expect(screen.getByTestId("scoping-value-10")).toHaveTextContent("1:1");
  });

  it("confirm-and-build is gated until no rows are flagged, then writes + navigates", async () => {
    const posted: any[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input, init) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1")]));
      if (url.endsWith("/api/samples/tags/batch")) {
        posted.push(init?.body ? JSON.parse(String(init.body)) : null);
        return Promise.resolve(jsonRes([{ id: 1, sample_id: 10, key: "ratio", value: "1:1", source: "scoping" }], 201));
      }
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    // Build is gated until the queries resolve + the (single, non-flagged)
    // row seeds → button enabled.
    await waitFor(() =>
      expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled());
    fireEvent.click(screen.getByTestId("scoping-open-confirm"));
    fireEvent.click(await screen.findByTestId("scoping-confirm-build"));
    await waitFor(() =>
      expect(posted[0]).toMatchObject({ key: "ratio", source: "scoping",
        tags: [{ sample_id: 10, value: "1:1" }] }));
    // D1: post-confirm navigates to the folio (/series), not /series/:id.
    await waitFor(() => expect(screen.getByTestId("folio-stub")).toBeInTheDocument());
  });

  it("does NOT navigate when the batch write fails (no silent data loss)", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1")]));
      if (url.endsWith("/api/samples/tags/batch"))
        // 400 = validation error → not retried (settles to .error immediately),
        // unlike a 5xx which the queue retries. Mirrors a route-rejected batch.
        return Promise.resolve(jsonRes({ error: "bad batch" }, 400));
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    await waitFor(() =>
      expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled());
    fireEvent.click(screen.getByTestId("scoping-open-confirm"));
    fireEvent.click(await screen.findByTestId("scoping-confirm-build"));
    // The mutation settles to error; navigation must NOT fire (no data loss).
    await waitFor(() => expect(screen.getByTestId("scoping-error-banner")).toBeInTheDocument());
    expect(screen.queryByTestId("folio-stub")).toBeNull();
  });

  it("D5: a sample without the ordering value is a loose match, omitted from the batch", async () => {
    const posted: any[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input, init) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples"))
        // sample 11 has no ratio → loose match (excluded by default).
        return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B")]));
      if (url.endsWith("/api/samples/tags/batch")) {
        posted.push(init?.body ? JSON.parse(String(init.body)) : null);
        return Promise.resolve(jsonRes([], 201));
      }
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    await waitFor(() =>
      expect(screen.getByTestId("scoping-open-confirm")).not.toBeDisabled());
    // sample 11 is in Himalaya-also-found, not a member row.
    expect(screen.queryByTestId("scoping-row-11")).toBeNull();
    expect(screen.getByTestId("scoping-loose-11")).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("scoping-open-confirm"));
    fireEvent.click(await screen.findByTestId("scoping-confirm-build"));
    await waitFor(() => expect(posted[0]).toMatchObject({
      key: "ratio", source: "scoping", tags: [{ sample_id: 10, value: "1:1" }] }));
    expect((posted[0].tags as { sample_id: number }[]).some((t) => t.sample_id === 11)).toBe(false);
  });

  it("renders the worksheet plate with autogroup summary + serif title + member row", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples"))
        return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B", "2:1")]));
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    await waitFor(() => expect(screen.getByTestId("scoping-autogroup")).toBeInTheDocument());
    expect(screen.getByTestId("scoping-plate")).toBeInTheDocument();
    expect(screen.getByTestId("scoping-title")).toBeInTheDocument();
    await waitFor(() => expect(screen.getByTestId("scoping-row-10")).toBeInTheDocument());
  });

  it("adding a loose match folds it into the member list", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples"))
        return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1"), pickerRow(11, "B")]));
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    fireEvent.click(await screen.findByTestId("scoping-loose-add-11"));
    // 11 becomes a member row (still flagged — it has no value yet).
    await waitFor(() => expect(screen.getByTestId("scoping-row-11")).toBeInTheDocument());
    expect(screen.queryByTestId("scoping-loose-11")).toBeNull();
  });

  it("editing a value to empty flags the row; Undo restores it (in-session ⌘Z)", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/sample-tags")) return Promise.resolve(jsonRes([{ key: "ratio", value: "1:1" }]));
      if (url.endsWith("/api/picker-samples")) return Promise.resolve(jsonRes([pickerRow(10, "A", "1:1")]));
      return Promise.resolve(jsonRes([]));
    });
    renderScoping();
    // Re-open the confident value (confirm-not-fill-out) and clear it → flagged.
    fireEvent.click(await screen.findByTestId("scoping-value-10"));
    const input = screen.getByTestId("scoping-value-input-10");
    fireEvent.change(input, { target: { value: "" } });
    fireEvent.blur(input);
    await waitFor(() =>
      expect(screen.getByTestId("scoping-row-10")).toHaveAttribute("data-flagged", "true"));
    // Undo restores the confirmed value.
    fireEvent.click(screen.getByTestId("scoping-undo"));
    await waitFor(() =>
      expect(screen.getByTestId("scoping-row-10")).not.toHaveAttribute("data-flagged"));
    expect(screen.getByTestId("scoping-value-10")).toHaveTextContent("1:1");
  });

  it("discard navigates back to the folio", async () => {
    vi.spyOn(global, "fetch").mockImplementation(() => Promise.resolve(jsonRes([])));
    renderScoping();
    fireEvent.click(await screen.findByTestId("scoping-discard"));
    await waitFor(() => expect(screen.getByTestId("folio-stub")).toBeInTheDocument());
  });
});
