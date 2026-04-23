import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { NotesTab } from "../src/components/NotesTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined });
});

function mockSample(notes: string | null): ReturnType<typeof vi.spyOn> {
  return vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/experiments/1/samples")) {
      return new Response(JSON.stringify([{
        id: 10, experiment_id: 1, label: "A1", name: "s1", notes, tags: [],
      }]), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (url.endsWith("/api/samples/10")) {
      return new Response(JSON.stringify({
        id: 10, experiment_id: 1, label: "A1", name: "s1",
        notes: "new notes", tags: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<NotesTab>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<NotesTab />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("pre-fills the textarea with existing notes", async () => {
    useAppState.setState({ activeSampleId: 10 });
    mockSample("existing value");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox");
    await waitFor(() => expect((ta as HTMLTextAreaElement).value).toBe("existing value"));
  });

  it("PATCHes /api/samples/:id on blur when value changed", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const spy = mockSample("");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox") as HTMLTextAreaElement;
    fireEvent.change(ta, { target: { value: "new notes" } });
    fireEvent.blur(ta);
    await waitFor(() => {
      const patchCalls = spy.mock.calls.filter((c) => {
        const url = typeof c[0] === "string" ? c[0] : (c[0] as Request).url;
        const init = c[1] as RequestInit | undefined;
        return url.endsWith("/api/samples/10") && init?.method === "PATCH";
      });
      expect(patchCalls.length).toBeGreaterThan(0);
    });
  });

  it("does NOT PATCH on blur if the value is unchanged", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const spy = mockSample("existing");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox") as HTMLTextAreaElement;
    await waitFor(() => expect(ta.value).toBe("existing"));
    fireEvent.blur(ta);
    await new Promise((r) => setTimeout(r, 10));
    const patches = spy.mock.calls.filter((c) =>
      (c[1] as RequestInit | undefined)?.method === "PATCH");
    expect(patches.length).toBe(0);
  });
});
