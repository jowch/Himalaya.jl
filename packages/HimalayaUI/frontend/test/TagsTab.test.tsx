import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { TagsTab } from "../src/components/TagsTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined });
});

function mockSamples(samples: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(samples), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<TagsTab>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<TagsTab />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("renders existing tags as chips with a remove button each", async () => {
    useAppState.setState({ activeSampleId: 10 });
    mockSamples([{
      id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null,
      tags: [
        { id: 1, key: "lipid",   value: "DOPC",     source: "manifest" },
        { id: 2, key: "peptide", value: "melittin", source: "manual"   },
      ],
    }]);
    renderWithProviders(<TagsTab />);
    await waitFor(() =>
      expect(screen.getByText(/lipid/i)).toBeInTheDocument(),
    );
    expect(screen.getByText(/DOPC/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /remove lipid/i })).toBeInTheDocument();
  });

  it("submitting the add form posts {key, value}", async () => {
    const user = userEvent.setup();
    useAppState.setState({ activeSampleId: 10 });
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/experiments/1/samples")) {
        return new Response(JSON.stringify([{
          id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [],
        }]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url.endsWith("/api/samples/10/tags")) {
        return new Response(JSON.stringify({
          id: 5, sample_id: 10, key: "lipid", value: "DOPC", source: "manual",
        }), { status: 201, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });

    renderWithProviders(<TagsTab />);
    await user.type(screen.getByPlaceholderText(/key/i), "lipid");
    await user.type(screen.getByPlaceholderText(/value/i), "DOPC");
    fireEvent.click(screen.getByRole("button", { name: /add tag/i }));

    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      const postUrls = urls.filter((u) => u === "/api/samples/10/tags");
      expect(postUrls.length).toBeGreaterThan(0);
    });
  });

  it("clicking remove on an existing tag calls DELETE", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      const method = typeof input === "string" ? "GET" : (input as Request).method;
      if (url.endsWith("/api/experiments/1/samples") && method === "GET") {
        return new Response(JSON.stringify([{
          id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null,
          tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }],
        }]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url.endsWith("/api/samples/10/tags/1")) {
        return new Response(null, { status: 204 });
      }
      return new Response("not found", { status: 404 });
    });

    renderWithProviders(<TagsTab />);
    const rm = await screen.findByRole("button", { name: /remove lipid/i });
    fireEvent.click(rm);
    await waitFor(() => {
      const calls = fetchSpy.mock.calls.filter((c) => {
        const url = typeof c[0] === "string" ? c[0] : (c[0] as Request).url;
        return url === "/api/samples/10/tags/1";
      });
      expect(calls.length).toBeGreaterThan(0);
    });
  });
});
