import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor, within } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { FolderPicker, crumbsOf, describeHint } from "../src/components/FolderPicker";
import type { FolderPickerMode } from "../src/components/FolderPicker";
import * as api from "../src/api";

const TREE: Record<string, string[]> = {
  "/srv/beamtimes": ["2026-10", "2026-11", "calibrants"],
  "/srv/beamtimes/2026-11": ["setup1", "setup2", "notes"],
  "/srv/beamtimes/calibrants": ["agbe"],
  "/srv/beamtimes/2026-11/setup1": ["analysis", "data"],
  "/srv": ["beamtimes"],
  "/": ["srv"],
};
const KINDS: Record<string, api.FolderDescription["kind"]> = {
  "/srv/beamtimes/2026-10": "beamtime",
  "/srv/beamtimes/2026-11": "beamtime",
  "/srv/beamtimes/2026-11/setup1": "setup",
  "/srv/beamtimes/2026-11/setup2": "setup",
};

const parentOf = (p: string): string | null => (p === "/" ? null : p.slice(0, p.lastIndexOf("/")) || "/");

function fakeFs(): void {
  vi.spyOn(api, "listFolder").mockImplementation(async (p) => {
    const path = p === "" ? "/srv/beamtimes" : p.replace(/\/+$/, "") || "/";
    const folders = TREE[path];
    if (folders === undefined) throw new api.ApiError(400, "that folder does not exist", null);
    return { path, parent: parentOf(path), folders, inside_setup: KINDS[path] === "setup" };
  });
  vi.spyOn(api, "describeFolder").mockImplementation(async (p) => {
    const kind = KINDS[p] ?? null;
    return { path: p, kind, setups: kind === "beamtime" ? 2 : null, added: p.endsWith("2026-10") ? { id: 3, name: "October" } : null };
  });
}

function renderPicker(mode: FolderPickerMode = "beamtime", start = "", onCancel?: () => void) {
  const onChoose = vi.fn();
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <FolderPicker mode={mode} start={start} onChoose={onChoose} onCancel={onCancel} autoFocus />
    </QueryClientProvider>,
  );
  return { onChoose, list: screen.getByTestId("folder-browser-list") };
}

const names = (list: HTMLElement): string[] =>
  within(list).queryAllByRole("option").map((o) => o.querySelector("span.truncate")?.textContent ?? "");
const selected = (list: HTMLElement): string | null =>
  within(list).queryAllByRole("option").find((o) => o.getAttribute("aria-selected") === "true")?.getAttribute("data-path") ?? null;

describe("FolderPicker", () => {
  beforeEach(fakeFs);
  afterEach(() => vi.restoreAllMocks());

  it("opens where the server says, lists folders, and labels them one at a time", async () => {
    const { list } = renderPicker();
    await waitFor(() => expect(names(list)).toEqual(["2026-10", "2026-11", "calibrants"]));
    expect(screen.getByRole("navigation", { name: "Location" })).toHaveTextContent("/›srv›beamtimes");
    await within(list).findByRole("option", { name: "2026-10, Beamtime · 2 setups · added as October" });
    await within(list).findByRole("option", { name: "2026-11, Beamtime · 2 setups" });
    await within(list).findByRole("option", { name: "calibrants" });
    expect(vi.mocked(api.describeFolder).mock.calls.map((c) => c[0])).toEqual([
      "/srv/beamtimes/2026-10", "/srv/beamtimes/2026-11", "/srv/beamtimes/calibrants",
    ]);
  });

  it("Enter on a beamtime picks it, Enter on a plain folder opens it", async () => {
    const { list, onChoose } = renderPicker();
    await within(list).findByRole("option", { name: /^calibrants/ });
    fireEvent.keyDown(list, { key: "Tab", shiftKey: true });
    expect(selected(list)).toBe("/srv/beamtimes/calibrants");
    fireEvent.keyDown(list, { key: "Enter" });
    await waitFor(() => expect(names(list)).toEqual(["agbe"]));
    expect(onChoose).not.toHaveBeenCalled();
    fireEvent.keyDown(list, { key: "Escape" });
    await waitFor(() => expect(selected(list)).toBe("/srv/beamtimes/calibrants"));
    fireEvent.keyDown(list, { key: "ArrowUp" });
    fireEvent.keyDown(list, { key: "Enter" });
    await waitFor(() => expect(onChoose).toHaveBeenCalledWith("/srv/beamtimes/2026-11"));
  });

  it("Enter picks a beamtime even before its label has arrived", async () => {
    const { list, onChoose } = renderPicker();
    await waitFor(() => expect(names(list)).toContain("2026-10"));
    fireEvent.keyDown(list, { key: "Enter" });
    await waitFor(() => expect(onChoose).toHaveBeenCalledWith("/srv/beamtimes/2026-10"));
  });

  it("→ opens a beamtime to pick one setup, and ← goes back up to it", async () => {
    const { list, onChoose } = renderPicker();
    await waitFor(() => expect(names(list)).toContain("2026-11"));
    fireEvent.keyDown(list, { key: "ArrowDown" });
    fireEvent.keyDown(list, { key: "ArrowRight" });
    await waitFor(() => expect(names(list)).toEqual(["setup1", "setup2", "notes"]));
    await within(list).findByRole("option", { name: "setup1, Setup" });
    fireEvent.keyDown(list, { key: "Enter" });
    await waitFor(() => expect(onChoose).toHaveBeenCalledWith("/srv/beamtimes/2026-11/setup1"));
    fireEvent.keyDown(list, { key: "ArrowLeft" });
    await waitFor(() => expect(selected(list)).toBe("/srv/beamtimes/2026-11"));
  });

  it("typing narrows the list by name, Backspace and Esc undo it", async () => {
    const { list } = renderPicker();
    await waitFor(() => expect(names(list)).toHaveLength(3));
    fireEvent.keyDown(list, { key: "C" });
    fireEvent.keyDown(list, { key: "a" });
    expect(names(list)).toEqual(["calibrants"]);
    expect(screen.getByTestId("folder-browser-filter")).toHaveTextContent("Showing names with “Ca”");
    fireEvent.keyDown(list, { key: "l" });
    fireEvent.keyDown(list, { key: "x" });
    expect(screen.getByTestId("folder-browser-status")).toHaveTextContent("No folder here has “Calx” in its name.");
    fireEvent.keyDown(list, { key: "Backspace" });
    expect(names(list)).toEqual(["calibrants"]);
    fireEvent.keyDown(list, { key: "Escape" });
    expect(names(list)).toHaveLength(3);
  });

  it("/ opens a path field, Enter there goes to the folder or picks a beamtime, Esc closes it", async () => {
    const { list, onChoose } = renderPicker();
    await waitFor(() => expect(names(list)).toHaveLength(3));
    fireEvent.keyDown(list, { key: "/" });
    const field = screen.getByRole("textbox", { name: "Folder path" });
    expect(field).toHaveValue("/srv/beamtimes/");
    fireEvent.change(field, { target: { value: "/nope" } });
    fireEvent.keyDown(field, { key: "Enter" });
    expect(await screen.findByRole("alert")).toHaveTextContent("That folder does not exist.");
    fireEvent.change(field, { target: { value: "/srv" } });
    fireEvent.keyDown(field, { key: "Enter" });
    await waitFor(() => expect(names(screen.getByTestId("folder-browser-list"))).toEqual(["beamtimes"]));
    expect(onChoose).not.toHaveBeenCalled();
    fireEvent.keyDown(screen.getByTestId("folder-browser-list"), { key: "/" });
    fireEvent.change(screen.getByRole("textbox", { name: "Folder path" }), { target: { value: "/srv/beamtimes/2026-11/" } });
    fireEvent.keyDown(screen.getByRole("textbox", { name: "Folder path" }), { key: "Enter" });
    await waitFor(() => expect(onChoose).toHaveBeenCalledWith("/srv/beamtimes/2026-11"));
    fireEvent.keyDown(screen.getByTestId("folder-browser-list"), { key: "/" });
    fireEvent.keyDown(screen.getByRole("textbox", { name: "Folder path" }), { key: "Escape" });
    expect(screen.queryByRole("textbox", { name: "Folder path" })).toBeNull();
  });

  it("a click on a location segment goes there", async () => {
    const { list } = renderPicker();
    await waitFor(() => expect(names(list)).toHaveLength(3));
    fireEvent.click(screen.getByRole("button", { name: "srv" }));
    await waitFor(() => expect(names(list)).toEqual(["beamtimes"]));
  });

  it("stops labelling a folder's rows once you leave it", async () => {
    const pending: Array<(d: api.FolderDescription) => void> = [];
    vi.mocked(api.describeFolder).mockImplementation((p, opts) => new Promise((resolve, reject) => {
      opts?.signal?.addEventListener("abort", () => reject(new DOMException("aborted", "AbortError")));
      pending.push(() => resolve({ path: p, kind: null, setups: null, added: null }));
    }));
    const { list } = renderPicker();
    await waitFor(() => expect(pending).toHaveLength(1));
    expect(within(list).getByRole("option", { name: "2026-10, Checking…" })).toBeInTheDocument();
    fireEvent.keyDown(list, { key: "ArrowLeft" });
    await waitFor(() => expect(names(list)).toEqual(["beamtimes"]));
    await waitFor(() => expect(api.describeFolder).toHaveBeenLastCalledWith("/srv/beamtimes", expect.anything()));
    expect(vi.mocked(api.describeFolder).mock.calls.map((c) => c[0])).toEqual(["/srv/beamtimes/2026-10", "/srv/beamtimes"]);
  });

  it("in folder mode Enter uses the highlighted folder without checking what it holds", async () => {
    const { list, onChoose } = renderPicker("any", "/srv/beamtimes/calibrants");
    await waitFor(() => expect(names(list)).toEqual(["agbe"]));
    expect(within(list).getByRole("option", { name: "agbe" })).toBeInTheDocument();
    fireEvent.keyDown(list, { key: "Enter" });
    expect(onChoose).toHaveBeenCalledWith("/srv/beamtimes/calibrants/agbe");
    expect(api.describeFolder).not.toHaveBeenCalled();
  });

  it("opens the nearest folder above one that no longer exists, and ← works from there", async () => {
    const { list } = renderPicker("beamtime", "/srv/beamtimes/gone/deeper");
    await waitFor(() => expect(names(list)).toEqual(["2026-10", "2026-11", "calibrants"]));
    fireEvent.keyDown(list, { key: "ArrowLeft" });
    await waitFor(() => expect(names(list)).toEqual(["beamtimes"]));
  });

  it("inside a setup labels nothing and Enter opens the highlighted folder", async () => {
    const { list, onChoose } = renderPicker("beamtime", "/srv/beamtimes/2026-11/setup1");
    await waitFor(() => expect(names(list)).toEqual(["analysis", "data"]));
    expect(within(list).getByRole("option", { name: "analysis" })).toBeInTheDocument();
    fireEvent.keyDown(list, { key: "Enter" });
    await waitFor(() => expect(api.listFolder).toHaveBeenLastCalledWith("/srv/beamtimes/2026-11/setup1/analysis", expect.anything()));
    expect(onChoose).not.toHaveBeenCalled();
    expect(api.describeFolder).not.toHaveBeenCalled();
  });

  it("with a cancel handler, Esc clears the filter first and then cancels instead of going up", async () => {
    const onCancel = vi.fn();
    const { list } = renderPicker("any", "/srv/beamtimes", onCancel);
    await waitFor(() => expect(names(list)).toHaveLength(3));
    fireEvent.keyDown(list, { key: "c" });
    fireEvent.keyDown(list, { key: "Escape" });
    expect(onCancel).not.toHaveBeenCalled();
    expect(names(list)).toHaveLength(3);
    fireEvent.keyDown(list, { key: "Escape" });
    expect(onCancel).toHaveBeenCalledTimes(1);
    expect(names(list)).toHaveLength(3);
  });

  it("says when a folder has no folders", async () => {
    TREE["/srv/empty"] = [];
    const { list } = renderPicker("beamtime", "/srv/empty");
    await waitFor(() => expect(screen.getByTestId("folder-browser-status")).toHaveTextContent("No folders here."));
    expect(names(list)).toEqual([]);
    delete TREE["/srv/empty"];
    fireEvent.click(screen.getByRole("button", { name: /Type or paste a path/ }));
    expect(screen.getByRole("textbox", { name: "Folder path" })).toHaveValue("/srv/empty/");
  });
});

describe("crumbsOf and describeHint", () => {
  it("split a path into clickable segments", () => {
    expect(crumbsOf("/srv/beamtimes")).toEqual([
      { label: "/", path: "/" }, { label: "srv", path: "/srv" }, { label: "beamtimes", path: "/srv/beamtimes" },
    ]);
    expect(crumbsOf("/")).toEqual([{ label: "/", path: "/" }]);
  });

  it("name beamtimes, setups and where they were added", () => {
    expect(describeHint(undefined)).toBe("Checking…");
    expect(describeHint({ path: "/a", kind: "beamtime", setups: 1, added: null })).toBe("Beamtime · 1 setup");
    expect(describeHint({ path: "/a", kind: "setup", setups: null, added: { id: 1, name: "Nov" } })).toBe("Setup · added as Nov");
    expect(describeHint({ path: "/a", kind: null, setups: null, added: null })).toBe("");
  });
});
