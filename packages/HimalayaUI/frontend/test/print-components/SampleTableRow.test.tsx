import { render, screen, fireEvent, within } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import type { ReactNode } from "react";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";
import type { GalleryExposure } from "../../src/print/components/ThumbnailGallery";
import type { Tag } from "../../src/print/ui/tag";
import {
  useRovingGrid,
  RovingGridProvider,
} from "../../src/lib/grid/useRovingGrid";

/** Wrap a row in a LIVE roving-grid context (the grid is "on") so the row's
 *  gridcell roles + roving tabindex wiring are exercised. Mirrors how SheetTable
 *  supplies the provider; multi-widget cols Exposures(2) + Tags(4). */
function RovingHarness({ children, rows = 3 }: { children: ReactNode; rows?: number }) {
  const grid = useRovingGrid({ rows, cols: 6, interactionCols: new Set([2, 4]) });
  return <RovingGridProvider value={grid}>{children}</RovingGridProvider>;
}

const EXPOSURES: GalleryExposure[] = [
  { id: 37, src: null, frameNo: "37", representative: true },
  { id: 65, src: null, frameNo: "65" },
  { id: 66, src: null, frameNo: "66", rejected: true },
];

const TAGS: Tag[] = [
  { key: "LL37" },
  { key: "temperature", value: "37C" },
];

function renderRow(overrides: Partial<Parameters<typeof SampleTableRow>[0]> = {}) {
  return render(
    <SampleTableRow
      name="Sample A"
      sampleId="s-001"
      exposures={EXPOSURES}
      kept={2}
      total={3}
      tags={TAGS}
      {...overrides}
    />
  );
}

describe("<SampleTableRow> structure", () => {
  it("(a) renders the root data-testid=sample-table-row", () => {
    renderRow();
    expect(screen.getByTestId("sample-table-row")).toBeInTheDocument();
  });

  it("(b) contains the SpecCell name and sampleId text", () => {
    renderRow();
    expect(screen.getByText("Sample A")).toBeInTheDocument();
    expect(screen.getByText("s-001")).toBeInTheDocument();
  });

  it("(c) contains a thumbnail-gallery with one thumbnail per exposure", () => {
    renderRow();
    const gallery = screen.getByTestId("thumbnail-gallery");
    expect(gallery).toBeInTheDocument();
    const thumbs = screen.getAllByTestId("thumbnail");
    expect(thumbs).toHaveLength(EXPOSURES.length);
  });

  it("renders exactly 5 grid-child cell wrappers", () => {
    renderRow();
    const row = screen.getByTestId("sample-table-row");
    const grid = row.firstElementChild!;
    expect(grid.children).toHaveLength(5);
  });

  it("exposes role=row with 5 cells (no checkbox)", () => {
    renderRow();
    const row = screen.getByRole("row");
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("role", "row");
    expect(within(row).getAllByRole("cell")).toHaveLength(5);
  });

  it("exposes 6 cells when a checkbox column is rendered (onCheck set)", () => {
    renderRow({ onCheck: () => {} });
    const row = screen.getByRole("row");
    expect(within(row).getAllByRole("cell")).toHaveLength(6);
  });

  it("sticky identity cells: Sample sticks at left 0 when there is no checkbox column", () => {
    renderRow();
    const cells = within(screen.getByRole("row")).getAllByRole("cell");
    expect(cells[0]).toHaveAttribute("data-sticky", "true");
    expect(cells[0]).toHaveStyle({ left: "0px" });
    cells.slice(1).forEach((c) => expect(c).not.toHaveAttribute("data-sticky"));
  });

  it("sticky identity cells: checkbox at 0 + Sample offset by the checkbox track (shared const)", () => {
    renderRow({ onCheck: () => {} });
    const cells = within(screen.getByRole("row")).getAllByRole("cell");
    expect(cells[0]).toHaveAttribute("data-sticky", "true"); // checkbox cell
    expect(cells[1]).toHaveAttribute("data-sticky", "true"); // Sample cell
    expect(cells[1]).toHaveStyle({ left: "36px" });
    cells.slice(2).forEach((c) => expect(c).not.toHaveAttribute("data-sticky"));
  });

  it("(d) contains the KeptCell counts (kept and total)", () => {
    const { container } = renderRow();
    expect(container.querySelector("[data-role='kept-count']")).toHaveTextContent("2");
    expect(container.querySelector("[data-role='kept-total']")).toHaveTextContent("/ 3");
  });

  it("(h) renders the tags via TagList", () => {
    renderRow();
    expect(screen.getByText("LL37")).toBeInTheDocument();
    expect(screen.getByText("temperature")).toBeInTheDocument();
  });
});

describe("<SampleTableRow> status column", () => {
  it("(e) phase set → PhaseChip present and status-unset absent", () => {
    const { container } = renderRow({ phase: "Pn3m" });
    expect(screen.getByTestId("phase-chip")).toBeInTheDocument();
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("(e) phase null → status-unset present", () => {
    const { container } = renderRow({ phase: null });
    expect(container.querySelector("[data-role='status-unset']")).toBeInTheDocument();
  });
});

describe("<SampleTableRow> screened reflection", () => {
  it("(f) data-screened is 'true' when screened", () => {
    renderRow({ screened: true });
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "true");
  });

  it("(f) data-screened is 'false' when not screened", () => {
    renderRow({ screened: false });
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "false");
  });

  it("(f) data-screened is 'false' when screened is omitted", () => {
    renderRow();
    expect(screen.getByTestId("sample-table-row")).toHaveAttribute("data-screened", "false");
  });
});

describe("<SampleTableRow> exposure selection", () => {
  it("(g) selecting a thumb bubbles onSelectExposure with that exposure id", () => {
    const onSelectExposure = vi.fn();
    renderRow({ onSelectExposure });
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[1]);
    expect(onSelectExposure).toHaveBeenCalledTimes(1);
    expect(onSelectExposure).toHaveBeenCalledWith(65);
  });

  it("selectedExposureIds highlights EVERY matching thumb (multi-select cull model)", () => {
    // EXPOSURES order: [37, 65, 66] — select the first and last.
    renderRow({ selectedExposureIds: new Set([37, 66]) });
    const thumbs = screen.getAllByTestId("thumbnail");
    expect((thumbs[0].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
    expect((thumbs[1].getAttribute("data-state") ?? "").split(" ")).not.toContain("selected");
    expect((thumbs[2].getAttribute("data-state") ?? "").split(" ")).toContain("selected");
  });
});

// ---------------------------------------------------------------------------
// Nav seams (Task 3)
// ---------------------------------------------------------------------------

const baseProps = {
  name: "POPC",
  sampleId: "#42",
  exposures: [{ id: 1, src: null, frameNo: 1 }],
  kept: 1,
  total: 1,
  tags: [] as import("../../src/print/ui/tag").Tag[],
};

describe("<SampleTableRow> nav seams", () => {
  it("forwards onOpenLoupe to the SpecCell name button", () => {
    const onOpenLoupe = vi.fn();
    render(<SampleTableRow {...baseProps} onOpenLoupe={onOpenLoupe} />);
    fireEvent.click(screen.getByRole("button", { name: /POPC/ }));
    expect(onOpenLoupe).toHaveBeenCalled();
  });

  it("forwards onActivateExposure to the gallery on thumb double-click", () => {
    const onActivateExposure = vi.fn();
    render(<SampleTableRow {...baseProps} onActivateExposure={onActivateExposure} />);
    fireEvent.doubleClick(screen.getAllByTestId("thumbnail")[0]);
    expect(onActivateExposure).toHaveBeenCalledWith(1);
  });

  it("makes the status cell a Focus door (button) when onOpenFocus is set", () => {
    const onOpenFocus = vi.fn();
    render(<SampleTableRow {...baseProps} phase={null} onOpenFocus={onOpenFocus} />);
    fireEvent.click(screen.getByRole("button", { name: /index/i }));
    expect(onOpenFocus).toHaveBeenCalled();
  });

  it("status cell is NOT a button when onOpenFocus is absent", () => {
    render(<SampleTableRow {...baseProps} phase={null} />);
    expect(screen.queryByRole("button", { name: /index/i })).toBeNull();
  });
});

// ---------------------------------------------------------------------------
// Roving data grid (SA-ROVING)
// ---------------------------------------------------------------------------

const rovingProps = {
  ...baseProps,
  rowIndex: 1,
  onCheck: () => {},
  onOpenLoupe: () => {},
  onOpenFocus: () => {},
};

describe("<SampleTableRow> roving grid wiring", () => {
  it("renders gridcell roles (not cell) when inside a roving context with a rowIndex", () => {
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} phase={null} />
      </RovingHarness>,
    );
    const row = screen.getByRole("row");
    expect(within(row).getAllByRole("gridcell")).toHaveLength(6);
    expect(within(row).queryAllByRole("cell")).toHaveLength(0);
  });

  it("renders role=cell (unchanged) when OUTSIDE a roving context, even with a rowIndex", () => {
    render(<SampleTableRow {...rovingProps} phase={null} />);
    const row = screen.getByRole("row");
    // Inert context default → tabIndexFor undefined → static path.
    expect(within(row).getAllByRole("cell")).toHaveLength(6);
    expect(within(row).queryAllByRole("gridcell")).toHaveLength(0);
  });

  it("the active (first) row exposes exactly one tabindex=0 widget; the rest are -1", () => {
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} phase={null} />
      </RovingHarness>,
    );
    const row = screen.getByRole("row");
    const tabbable = row.querySelectorAll('[tabindex="0"]');
    expect(tabbable).toHaveLength(1);
    // Default active cell is the Sample widget (col 1).
    expect(tabbable[0]).toHaveAttribute("data-role", "spec-name");
  });

  it("a NON-active row (rowIndex 2) has no tabindex=0 widget", () => {
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} rowIndex={2} phase={null} />
      </RovingHarness>,
    );
    const row = screen.getByRole("row");
    expect(row.querySelectorAll('[tabindex="0"]')).toHaveLength(0);
  });

  it("clicking a thumbnail STILL bubbles onSelectExposure under roving (pointer parity)", () => {
    const onSelectExposure = vi.fn();
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} phase={null} onSelectExposure={onSelectExposure} />
      </RovingHarness>,
    );
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]);
    expect(onSelectExposure).toHaveBeenCalledWith(1);
  });

  it("the Focus-door button still fires onOpenFocus under roving (pointer parity)", () => {
    const onOpenFocus = vi.fn();
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} phase={null} onOpenFocus={onOpenFocus} />
      </RovingHarness>,
    );
    fireEvent.click(screen.getByRole("button", { name: /index/i }));
    expect(onOpenFocus).toHaveBeenCalled();
  });

  it("mousedown on a cell makes it the active roving cell (pointer parity, BUG D model)", () => {
    // Pointer parity is now wired on mousedown (not onFocus / onClick): pressing
    // on the Kept cell (col 3) moves the single tabindex=0 to that gridcell. We
    // assert the OBSERVABLE roving STATE shift (tabindex), not live focus — the
    // no-focus-yank / Tab-out behavior is the render-verify / e2e gate.
    render(
      <RovingHarness>
        <SampleTableRow {...rovingProps} phase={null} />
      </RovingHarness>,
    );
    const row = screen.getByRole("row");
    const cells = within(row).getAllByRole("gridcell");
    // 6 cells: checkbox(0) Sample(1) Exposures(2) Kept(3) Tags(4) Status(5).
    const keptCell = cells[3];
    fireEvent.mouseDown(keptCell);
    expect(keptCell).toHaveAttribute("tabindex", "0");
    // The Sample widget (was active) is no longer the tab stop.
    expect(within(row).getByText("POPC").closest("button")).toHaveAttribute("tabindex", "-1");
  });

  it("multi-widget cell descendants are made inert (tabindex=-1) in nav mode", () => {
    render(
      <RovingHarness>
        <SampleTableRow
          {...rovingProps}
          rowIndex={2}
          phase={null}
          exposures={[
            { id: 1, src: null, frameNo: "1" },
            { id: 2, src: null, frameNo: "2" },
          ]}
        />
      </RovingHarness>,
    );
    // The two thumbnails (gallery cell focusables) must not be extra tab stops.
    screen.getAllByTestId("thumbnail").forEach((t) => {
      expect(t).toHaveAttribute("tabindex", "-1");
    });
  });
});
