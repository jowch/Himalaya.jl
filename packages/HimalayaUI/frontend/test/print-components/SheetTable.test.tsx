import { describe, it, expect, vi } from "vitest";
import { render, screen, within, fireEvent } from "@testing-library/react";
import { SheetTable } from "../../src/print/components/SheetTable";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";

const row = (id: string) => (
  <SampleTableRow key={id} name={`Sample ${id}`} sampleId={id} exposures={[]} kept={1} total={1} tags={[]} />
);

/** A row wired the way SamplesPage wires it under roving: checkbox column +
 *  rowIndex + the nav-seam handlers (so each cell has a real focusable widget). */
const gridRow = (id: string, rowIndex: number) => (
  <SampleTableRow
    key={id}
    rowIndex={rowIndex}
    name={`Sample ${id}`}
    sampleId={id}
    exposures={[{ id: 1, src: null, frameNo: "1" }]}
    kept={1}
    total={1}
    tags={[{ key: "t" }]}
    onCheck={() => {}}
    onOpenLoupe={() => {}}
    onOpenFocus={() => {}}
  />
);

const checkRow = (id: string) => (
  <SampleTableRow
    key={id}
    name={`Sample ${id}`}
    sampleId={id}
    exposures={[]}
    kept={1}
    total={1}
    tags={[]}
    onCheck={() => {}}
  />
);

describe("SheetTable", () => {
  it("renders the five aligned column headers", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const head = screen.getByTestId("sheet-head");
    ["Sample", "Exposures", "Frames kept", "Tags", "Status"].forEach((label) =>
      expect(head.textContent).toContain(label),
    );
  });

  it("slots the provided rows", () => {
    render(<SheetTable>{[row("a"), row("b"), row("c")]}</SheetTable>);
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(3);
  });

  it("renders the empty node when there are no rows", () => {
    render(<SheetTable empty={<span>No samples</span>}>{[]}</SheetTable>);
    expect(screen.getByTestId("sheet-empty")).toBeInTheDocument();
    expect(screen.queryByTestId("sample-table-row")).toBeNull();
  });
});

describe("SheetTable a11y semantics (WCAG 1.3.1 / 4.1.2)", () => {
  it("exposes a labelled table region", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    expect(screen.getByRole("table", { name: "Samples" })).toBeInTheDocument();
  });

  it("exposes the five column headers by name (no checkbox column)", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const headers = screen.getAllByRole("columnheader").map((h) => h.textContent);
    expect(headers).toEqual(["Sample", "Exposures", "Frames kept", "Tags", "Status"]);
  });

  it("column-header kickers use the soft tone (WCAG 1.4.3 on the sunk header band)", () => {
    // ink-faint is 2.92:1 on paper-sunk — fails AA even for large text. The
    // header band is paper-sunk, so its informational labels must be ink-soft.
    render(<SheetTable checkboxColumn>{[checkRow("a")]}</SheetTable>);
    const head = screen.getByTestId("sheet-head");
    const kickers = head.querySelectorAll("[data-tone]");
    expect(kickers).toHaveLength(5);
    kickers.forEach((k) => expect(k).toHaveAttribute("data-tone", "soft"));
  });

  it("adds a 'Select' columnheader as column 0 when checkboxColumn is set", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a")]}</SheetTable>);
    const headers = screen.getAllByRole("columnheader");
    expect(headers).toHaveLength(6);
    expect(headers[0]).toHaveAccessibleName("Select");
  });

  it("rows count = header row + N body rows", () => {
    render(<SheetTable>{[row("a"), row("b"), row("c")]}</SheetTable>);
    // 1 header row + 3 body rows
    expect(screen.getAllByRole("row")).toHaveLength(4);
  });

  it("each body row has one cell per column header (5, no checkbox)", () => {
    render(<SheetTable>{[row("a"), row("b")]}</SheetTable>);
    const bodyRows = screen
      .getAllByRole("row")
      .filter((r) => within(r).queryAllByRole("columnheader").length === 0);
    expect(bodyRows).toHaveLength(2);
    bodyRows.forEach((r) => {
      expect(within(r).getAllByRole("cell")).toHaveLength(5);
    });
  });

  it("each body row has 6 cells when the checkbox column is present", () => {
    render(
      <SheetTable checkboxColumn>
        {[checkRow("a"), checkRow("b")]}
      </SheetTable>,
    );
    const bodyRows = screen
      .getAllByRole("row")
      .filter((r) => within(r).queryAllByRole("columnheader").length === 0);
    expect(bodyRows).toHaveLength(2);
    bodyRows.forEach((r) => {
      expect(within(r).getAllByRole("cell")).toHaveLength(6);
    });
  });

  it("table role tree is contiguous: rowgroups are DIRECT children of the table, rows of a rowgroup", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a"), checkRow("b")]}</SheetTable>);
    const table = screen.getByRole("table", { name: "Samples" });
    const rowgroups = screen.getAllByRole("rowgroup");
    expect(rowgroups).toHaveLength(2); // header + rows region
    rowgroups.forEach((rg) => {
      expect(rg.closest('[role="table"]')).toBe(table);
      // No role-bearing element between table and rowgroup: direct DOM child.
      expect(rg.parentElement).toBe(table);
    });
    screen.getAllByRole("row").forEach((r) => {
      expect(r.parentElement?.getAttribute("role")).toBe("rowgroup");
    });
  });
});

describe("SheetTable sortable columns (SA-SORT, WAI-ARIA sortable table)", () => {
  function getHeader(name: string): HTMLElement {
    return screen
      .getAllByRole("columnheader")
      .find((h) => within(h).queryByRole("button", { name })) ?? screen
      .getAllByRole("columnheader")
      .find((h) => h.textContent?.includes(name))!;
  }

  it("renders the four data headers as sortable buttons when onSort is given", () => {
    render(
      <SheetTable checkboxColumn onSort={() => {}} sort={{ key: null, dir: "asc" }}>
        {[checkRow("a")]}
      </SheetTable>,
    );
    for (const label of ["Sample", "Exposures", "Frames kept", "Status"]) {
      expect(screen.getByRole("button", { name: label })).toBeInTheDocument();
    }
    // Tags stays a non-button label (multi-valued → not sortable).
    expect(screen.queryByRole("button", { name: "Tags" })).toBeNull();
  });

  it("clicking a sortable header calls onSort with that column key", () => {
    const onSort = vi.fn();
    render(
      <SheetTable checkboxColumn onSort={onSort} sort={{ key: null, dir: "asc" }}>
        {[checkRow("a")]}
      </SheetTable>,
    );
    fireEvent.click(screen.getByRole("button", { name: "Exposures" }));
    expect(onSort).toHaveBeenCalledWith("exposures");
    fireEvent.click(screen.getByRole("button", { name: "Frames kept" }));
    expect(onSort).toHaveBeenCalledWith("kept");
  });

  it("aria-sort reflects the active column's direction; other sortable headers are 'none'", () => {
    render(
      <SheetTable checkboxColumn onSort={() => {}} sort={{ key: "sample", dir: "desc" }}>
        {[checkRow("a")]}
      </SheetTable>,
    );
    expect(getHeader("Sample")).toHaveAttribute("aria-sort", "descending");
    for (const label of ["Exposures", "Frames kept", "Status"]) {
      expect(getHeader(label)).toHaveAttribute("aria-sort", "none");
    }
  });

  it("ascending renders aria-sort='ascending' on the active column only", () => {
    render(
      <SheetTable checkboxColumn onSort={() => {}} sort={{ key: "kept", dir: "asc" }}>
        {[checkRow("a")]}
      </SheetTable>,
    );
    expect(getHeader("Frames kept")).toHaveAttribute("aria-sort", "ascending");
    expect(getHeader("Sample")).toHaveAttribute("aria-sort", "none");
  });

  it("non-sortable columns (checkbox, Tags) carry NO aria-sort", () => {
    render(
      <SheetTable checkboxColumn onSort={() => {}} sort={{ key: "sample", dir: "asc" }}>
        {[checkRow("a")]}
      </SheetTable>,
    );
    expect(screen.getAllByRole("columnheader")[0]).not.toHaveAttribute("aria-sort"); // checkbox
    expect(getHeader("Tags")).not.toHaveAttribute("aria-sort");
  });

  it("the static-header fallback still renders when onSort is absent (other consumers unchanged)", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    // No sort buttons; the labels remain plain kicker text.
    expect(screen.queryByRole("button", { name: "Sample" })).toBeNull();
    screen
      .getAllByRole("columnheader")
      .forEach((h) => expect(h).not.toHaveAttribute("aria-sort"));
  });
});

describe("SheetTable horizontal scroll + sticky identity (WCAG 1.4.10)", () => {
  it("header rowgroup and rows region share ONE sheet-scroll container", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const scroll = screen.getByTestId("sheet-scroll");
    expect(within(scroll).getByTestId("sheet-head")).toBeInTheDocument();
    expect(within(scroll).getByTestId("sheet-rows")).toBeInTheDocument();
  });

  it("the sheet-scroll horizontal scroller carries tabindex=-1 (must NOT be an auto Tab-stop)", () => {
    // Chrome auto-adds an overflow-x-auto container with overflowing content to
    // the Tab order; tabindex=-1 opts the table's own scroller out so it never
    // traps focus before the roving grid is exited (SA-ROVING BUG D). Horizontal
    // scroll stays reachable via the grid's arrow-key cell navigation.
    render(<SheetTable>{[row("a")]}</SheetTable>);
    expect(screen.getByTestId("sheet-scroll")).toHaveAttribute("tabindex", "-1");
  });

  it("the scroller sits OUTSIDE the table role element (Card > scroller > table)", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const table = screen.getByRole("table", { name: "Samples" });
    const scroll = screen.getByTestId("sheet-scroll");
    expect(scroll.contains(table)).toBe(true);
    expect(table.contains(scroll)).toBe(false);
    // The Card stays the visual plate around the scroller.
    expect(screen.getByTestId("sheet-table").contains(scroll)).toBe(true);
  });

  it("sticky headers: checkbox + Sample carry data-sticky; the rest do not", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a")]}</SheetTable>);
    const headers = screen.getAllByRole("columnheader");
    expect(headers[0]).toHaveAccessibleName("Select");
    expect(headers[0]).toHaveAttribute("data-sticky", "true");
    expect(headers[1]).toHaveTextContent("Sample");
    expect(headers[1]).toHaveAttribute("data-sticky", "true");
    headers.slice(2).forEach((h) => expect(h).not.toHaveAttribute("data-sticky"));
  });

  it("sticky headers without checkbox column: only Sample carries data-sticky", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const headers = screen.getAllByRole("columnheader");
    expect(headers[0]).toHaveTextContent("Sample");
    expect(headers[0]).toHaveAttribute("data-sticky", "true");
    headers.slice(1).forEach((h) => expect(h).not.toHaveAttribute("data-sticky"));
  });

  it("sticky row cells: checkbox + Sample cells carry data-sticky; the rest do not", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a")]}</SheetTable>);
    const bodyRow = screen
      .getAllByRole("row")
      .find((r) => within(r).queryAllByRole("columnheader").length === 0)!;
    const cells = within(bodyRow).getAllByRole("cell");
    expect(cells[0]).toHaveAttribute("data-sticky", "true"); // checkbox
    expect(cells[1]).toHaveAttribute("data-sticky", "true"); // Sample
    cells.slice(2).forEach((c) => expect(c).not.toHaveAttribute("data-sticky"));
  });

  it("still renders the empty slot (outside the table role tree) when no children", () => {
    render(<SheetTable empty={<span>No samples</span>}>{[]}</SheetTable>);
    const empty = screen.getByTestId("sheet-empty");
    expect(empty).toBeInTheDocument();
    expect(empty.closest('[role="table"]')).toBeNull();
  });
});

describe("SheetTable roving data grid (SA-ROVING)", () => {
  function renderGrid() {
    return render(
      <SheetTable
        checkboxColumn
        roving
        dataRowCount={2}
        onSort={() => {}}
        sort={{ key: null, dir: "asc" }}
      >
        {[gridRow("a", 1), gridRow("b", 2)]}
      </SheetTable>,
    );
  }

  it("becomes a role=grid (not role=table) when roving is on", () => {
    renderGrid();
    expect(screen.getByRole("grid", { name: "Samples" })).toBeInTheDocument();
    expect(screen.queryByRole("table")).toBeNull();
  });

  it("body cells become role=gridcell (was role=cell)", () => {
    renderGrid();
    const bodyRows = screen
      .getAllByRole("row")
      .filter((r) => within(r).queryAllByRole("columnheader").length === 0);
    expect(bodyRows).toHaveLength(2);
    bodyRows.forEach((r) => {
      expect(within(r).getAllByRole("gridcell")).toHaveLength(6);
      expect(within(r).queryAllByRole("cell")).toHaveLength(0);
    });
  });

  it("focusOnMountRow lands focus on that data row's Sample cell on mount (LO-FOCUSRET)", () => {
    render(
      <SheetTable
        checkboxColumn
        roving
        dataRowCount={3}
        focusOnMountRow={2}
        onSort={() => {}}
        sort={{ key: null, dir: "asc" }}
      >
        {[gridRow("a", 1), gridRow("b", 2), gridRow("c", 3)]}
      </SheetTable>,
    );
    const bodyRows = screen
      .getAllByRole("row")
      .filter((r) => within(r).queryAllByRole("columnheader").length === 0);
    expect(bodyRows).toHaveLength(3);
    // Focus restored to the SECOND data row (row "b"), not <body> or row 1.
    expect(document.activeElement).not.toBe(document.body);
    expect(bodyRows[1]!.contains(document.activeElement)).toBe(true);
    expect(bodyRows[0]!.contains(document.activeElement)).toBe(false);
    // And it is the active (tabindex=0) Sample gridcell of that row.
    const active = document.activeElement as HTMLElement;
    expect(active.getAttribute("tabindex")).toBe("0");
  });

  it("without focusOnMountRow an ordinary mount steals no focus", () => {
    render(
      <SheetTable
        checkboxColumn
        roving
        dataRowCount={2}
        onSort={() => {}}
        sort={{ key: null, dir: "asc" }}
      >
        {[gridRow("a", 1), gridRow("b", 2)]}
      </SheetTable>,
    );
    expect(document.activeElement).toBe(document.body);
  });

  it("has EXACTLY one element with tabindex=0 in the grid at rest", () => {
    const { container } = renderGrid();
    const grid = screen.getByRole("grid", { name: "Samples" });
    const tabbable = grid.querySelectorAll('[tabindex="0"]');
    expect(tabbable).toHaveLength(1);
    // Everything else in the grid is tabindex=-1 (one tab stop).
    expect(grid.querySelectorAll('[tabindex="-1"]').length).toBeGreaterThan(0);
    expect(container).toBeTruthy();
  });

  it("the single tab stop defaults to the first data Sample cell's widget", () => {
    renderGrid();
    const grid = screen.getByRole("grid", { name: "Samples" });
    const active = grid.querySelector('[tabindex="0"]')!;
    // Sample col widget is the spec-name button.
    expect(active).toHaveAttribute("data-role", "spec-name");
  });

  it("keeps aria-sort headers under roving (SA-SORT preserved)", () => {
    renderGrid();
    const sample = screen
      .getAllByRole("columnheader")
      .find((h) => h.textContent?.includes("Sample"))!;
    expect(sample).toHaveAttribute("aria-sort", "none");
    // The sort buttons still exist and are clickable.
    expect(screen.getByRole("button", { name: "Exposures" })).toBeInTheDocument();
  });

  it("clicking a sortable header still sorts (pointer parity) under roving", () => {
    const onSort = vi.fn();
    render(
      <SheetTable checkboxColumn roving dataRowCount={1} onSort={onSort} sort={{ key: null, dir: "asc" }}>
        {[gridRow("a", 1)]}
      </SheetTable>,
    );
    fireEvent.click(screen.getByRole("button", { name: "Frames kept" }));
    expect(onSort).toHaveBeenCalledWith("kept");
  });

  it("the static (roving-off) path keeps role=table and role=cell", () => {
    render(<SheetTable checkboxColumn>{[row("a")]}</SheetTable>);
    expect(screen.getByRole("table", { name: "Samples" })).toBeInTheDocument();
    expect(screen.queryByRole("grid")).toBeNull();
  });
});
