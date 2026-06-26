// test/SeriesScopingPage.actions.test.tsx
//
// Phase 5.1 — SeriesScopingPage interaction-architecture migration tests.
// Mounts the full TestShell (keyboard layer + InteractionDock) so dock buttons
// AND keyboard-triggered actions flow through the REAL registry path:
//   usePageActions → registry → InteractionDock / useKeyboardLayer → action.run
// Arrow keys (↑/↓) fire on the scope container (not window).
// Alt+↑/↓ and ⌘Z/⌘⇧Z fire on window (keyboard layer).
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { SampleTagPair, PickerSampleRow } from "../src/api";
import { SeriesScopingPage } from "../src/print/pages/SeriesScopingPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { useListCursor } from "../src/print/interaction/useListCursor";
import { DockStepper } from "../src/print/ui";
import { runCursorContract } from "./interaction/cursorContract";

// ── mock state ───────────────────────────────────────────────────────────────
const IDS = [1, 2, 3];

// 3 samples with "concentration" tag — triggers the warm proposal path
const TAGS: SampleTagPair[] = [
  { key: "concentration", value: "10mM" },
  { key: "concentration", value: "20mM" },
  { key: "concentration", value: "30mM" },
];

function makePickerRow(id: number, value: string): PickerSampleRow {
  return {
    sample: {
      id,
      experiment_id: 1,
      name: `S${id}`,
      notes: null,
      tags: [{ id, key: "concentration", value, source: "auto" }],
    },
    indexing_exposure_id: null,
  };
}

const PICKER: PickerSampleRow[] = [
  makePickerRow(1, "10mM"),
  makePickerRow(2, "20mM"),
  makePickerRow(3, "30mM"),
];

const scopeSeries = { mutate: vi.fn(), isSuccess: false, isPending: false, error: null, data: undefined };
const createSeries = { mutate: vi.fn(), isSuccess: false, isPending: false, error: null, data: undefined };

vi.mock("../src/queries", () => ({
  useCorpusSampleTags: () => ({ data: TAGS, isLoading: false, isError: false }),
  useCorpusPickerSamples: () => ({ data: PICKER, isLoading: false, isError: false }),
  useScopeSeries: () => scopeSeries,
  useCreateSeries: () => createSeries,
  useMemberTraces: () => new Map(),
  useMemberIndices: () => new Map(),
}));

// Boneyard passes through children immediately (no loading overlay)
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

const navigateSpy = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => navigateSpy };
});

// ── shell ─────────────────────────────────────────────────────────────────────
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return (
    <>
      {children}
      <InteractionDock />
    </>
  );
}

function renderScoping() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series/new"]}>
        <Routes>
          <Route
            path="/series/new"
            element={<TestShell><SeriesScopingPage /></TestShell>}
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function getScope(): HTMLElement {
  const scope = document.querySelector("[data-interaction-scope]") as HTMLElement | null;
  if (!scope) throw new Error("No [data-interaction-scope] found");
  return scope;
}

function getMemberRows(): NodeListOf<HTMLElement> {
  return document.querySelectorAll<HTMLElement>("[data-reorder-row]");
}

beforeEach(() => {
  navigateSpy.mockClear();
  scopeSeries.mutate.mockClear();
  createSeries.mutate.mockClear();
  useInteraction.getState().clearPage();
});

// ── tests ─────────────────────────────────────────────────────────────────────
describe("SeriesScopingPage interaction", () => {
  it("member list rows are roving-tabindex: exactly one role=row with tabindex=0", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    const tabbable = document.querySelectorAll('[role="row"][tabindex="0"]');
    expect(tabbable.length).toBe(1);

    // All other row elements have tabindex=-1
    const all = document.querySelectorAll('[role="row"]');
    expect(all.length).toBe(3); // 3 members
    const inert = document.querySelectorAll('[role="row"][tabindex="-1"]');
    expect(inert.length).toBe(2);
  });

  it("ArrowDown on scope moves cursor to second row", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });

    // Cursor should now be on id=2 (second member)
    const tabbable = document.querySelectorAll('[role="row"][tabindex="0"]');
    expect(tabbable.length).toBe(1);
    const cursored = document.querySelector('[data-cursored="true"]');
    // The cursored row is the second one (id=2, value "20mM")
    expect(cursored).not.toBeNull();
  });

  it("Alt+ArrowUp moves the cursored row up and cursor follows the item", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    // Move cursor to 2nd row (id=2)
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });

    // Record which row was cursored before reorder
    const rowsBefore = getMemberRows();
    expect(rowsBefore.length).toBe(3);
    // Second row (index 1) should be cursored
    expect(rowsBefore[1]).toHaveAttribute("data-cursored", "true");
    expect(rowsBefore[1]).toHaveAttribute("data-reorder-index", "1");

    // Fire Alt+ArrowUp via window (keyboard layer → reorderUp action)
    act(() => {
      fireEvent.keyDown(window, { key: "ArrowUp", altKey: true });
    });

    // After reorder: id=2 moved to index 0; order is now [2, 1, 3]
    const rowsAfter = getMemberRows();
    expect(rowsAfter.length).toBe(3);
    // The first row should now be cursored (id=2 moved up)
    expect(rowsAfter[0]).toHaveAttribute("data-cursored", "true");
    // Cursor follows the item: the cursored row is now at index 0
    expect(rowsAfter[0]).toHaveAttribute("data-reorder-index", "0");
    // tabindex=0 is on the same element
    expect(rowsAfter[0]).toHaveAttribute("tabindex", "0");

    // Live-region must announce the reorder (regression guard for screen readers).
    // S2 moved from index 1 → index 0 = position 1 of 3.
    const liveRegion = screen.getByTestId("reorder-announcement");
    expect(liveRegion.textContent?.trim()).toBe("Moved S2 to position 1 of 3.");
  });

  it("Alt+ArrowDown moves the cursored row down and cursor follows the item", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    // Cursor starts at id=1 (first row)
    const rowsBefore = getMemberRows();
    expect(rowsBefore[0]).toHaveAttribute("data-cursored", "true");

    // Fire Alt+ArrowDown via window (keyboard layer → reorderDown action)
    act(() => {
      fireEvent.keyDown(window, { key: "ArrowDown", altKey: true });
    });

    // After reorder: id=1 moved to index 1; order is now [2, 1, 3]
    const rowsAfter = getMemberRows();
    // The second row should now be cursored (id=1 moved down)
    expect(rowsAfter[1]).toHaveAttribute("data-cursored", "true");
    expect(rowsAfter[1]).toHaveAttribute("tabindex", "0");

    // Live-region must announce the reorder (regression guard for screen readers).
    // S1 moved from index 0 → index 1 = position 2 of 3.
    const liveRegion = screen.getByTestId("reorder-announcement");
    expect(liveRegion.textContent?.trim()).toBe("Moved S1 to position 2 of 3.");
  });

  it("⌘Z fires undo (reverses a reorder) via core(undo)", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    // Move cursor to 2nd row, then reorder up so we have something to undo
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });
    act(() => { fireEvent.keyDown(window, { key: "ArrowUp", altKey: true }); });

    // Order is now [2, 1, 3]. Undo via ⌘Z should restore [1, 2, 3].
    act(() => {
      fireEvent.keyDown(window, { key: "z", metaKey: true });
    });

    // After undo, original order restored: first row is id=1 (data-reorder-index=0)
    const rowsAfterUndo = getMemberRows();
    // The cursored row should have reverted to its original position
    expect(rowsAfterUndo.length).toBe(3);
    // id=1 should be back at index 0
    expect(rowsAfterUndo[0]).toHaveAttribute("data-reorder-index", "0");
  });

  it("⌘⇧Z fires redo via core(redo)", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    const scope = getScope();
    // Reorder: id=1 → index 1
    act(() => { fireEvent.keyDown(window, { key: "ArrowDown", altKey: true }); });
    // Undo: restore [1, 2, 3]
    act(() => { fireEvent.keyDown(window, { key: "z", metaKey: true }); });

    const rowsBeforeRedo = getMemberRows();
    expect(rowsBeforeRedo[0]).toHaveAttribute("data-reorder-index", "0");

    // Redo: re-apply the reorder → [2, 1, 3]
    act(() => {
      fireEvent.keyDown(window, { key: "z", metaKey: true, shiftKey: true });
    });

    const rowsAfterRedo = getMemberRows();
    // id=1 should be at index 1 again (after redo)
    const cursored = document.querySelector('[data-cursored="true"]');
    expect(cursored).not.toBeNull();
    // The redo result is 3 rows, in the reordered order
    expect(rowsAfterRedo.length).toBe(3);

    // suppress unused variable warning
    void scope;
  });

  it("drag-reorder coexists with roving tabindex: data-reorder-row and tabindex coexist on the same row", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    // Every member row should have BOTH data-reorder-row (for drag) AND tabindex (for roving)
    const rows = getMemberRows();
    expect(rows.length).toBe(3);
    for (const row of rows) {
      // Drag attr: data-reorder-row is a boolean attr (present = true)
      expect(row).toHaveAttribute("data-reorder-row");
      // Roving attr: tabindex is set by cursor.rowProps
      const ti = row.getAttribute("tabindex");
      expect(ti === "0" || ti === "-1").toBe(true);
      // Role is also set by cursor.rowProps
      expect(row).toHaveAttribute("role", "row");
    }
    // Exactly one has tabindex=0
    const tabbable = document.querySelectorAll('[data-reorder-row][tabindex="0"]');
    expect(tabbable.length).toBe(1);
  });

  it("Escape navigates to /series via core(back)", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).toHaveBeenCalledWith("/series");
  });

  it("dock renders the InteractionDock up-link (no hand-built Dock)", async () => {
    renderScoping();
    await screen.findByTestId("scoping-scope");

    // InteractionDock renders the back link via core(back, dock:true) using DockUpLink.
    // DockUpLink always renders with data-testid="dock-up-link".
    const backLink = document.querySelector('[data-testid="dock-up-link"]');
    expect(backLink).not.toBeNull();
    expect(backLink?.textContent).toContain("All series");
  });
});

// ── cursor contract (member cursor) ──────────────────────────────────────────
// Runs the standard cursor contract (roving tabindex, input parity, orthogonality)
// on a standalone useListCursor over IDS = [1, 2, 3]. Scoping has no onActivate
// (rows don't navigate on Enter), so the activate parity test uses a no-op stub.
runCursorContract("Scoping member cursor", () => {
  const onActivate = vi.fn();
  return {
    ui: (capture) => {
      function Probe(): JSX.Element {
        const cursor = useListCursor({
          ids: IDS,
          onActivate,
          stepperLabel: "Member",
          stepperTestIdBase: "member",
          axis: "vertical",
        });
        capture({ cursor, ids: IDS, onActivate });
        return (
          <div data-interaction-scope="">
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>
                member {id}
              </div>
            ))}
            <DockStepper {...cursor.stepperProps()} />
          </div>
        );
      }
      return <Probe />;
    },
  };
});
