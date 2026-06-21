import { describe, it, expect, vi } from "vitest";
import { render, screen, within, fireEvent } from "@testing-library/react";
import { SheetTable } from "../../src/print/components/SheetTable";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";

const row = (id: string) => (
  <SampleTableRow key={id} name={`Sample ${id}`} sampleId={id} exposures={[]} kept={1} total={1} tags={[]} />
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
    ["Sample", "Exposures", "Frames kept", "Tags", "Phase"].forEach((label) =>
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
    expect(headers).toEqual(["Sample", "Exposures", "Frames kept", "Tags", "Phase"]);
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
    for (const label of ["Sample", "Exposures", "Frames kept", "Phase"]) {
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
    for (const label of ["Exposures", "Frames kept", "Phase"]) {
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
    // becomes a stray Tab stop. Horizontal scroll stays reachable via the
    // table's keyboard navigation.
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

describe("SheetTable plain table semantics", () => {
  it("renders a plain table (role=table, not role=grid)", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    expect(screen.queryByRole("grid")).toBeNull();
    expect(screen.getByRole("table")).toBeInTheDocument();
  });

  it("always renders role=table (never role=grid)", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a"), checkRow("b")]}</SheetTable>);
    expect(screen.getByRole("table", { name: "Samples" })).toBeInTheDocument();
    expect(screen.queryByRole("grid")).toBeNull();
  });

  it("body cells are always role=cell (never gridcell)", () => {
    render(<SheetTable checkboxColumn>{[checkRow("a"), checkRow("b")]}</SheetTable>);
    const bodyRows = screen
      .getAllByRole("row")
      .filter((r) => within(r).queryAllByRole("columnheader").length === 0);
    expect(bodyRows).toHaveLength(2);
    bodyRows.forEach((r) => {
      expect(within(r).getAllByRole("cell")).toHaveLength(6);
      expect(within(r).queryAllByRole("gridcell")).toHaveLength(0);
    });
  });
});
