import { describe, it, expect } from "vitest";
import { render, screen, within } from "@testing-library/react";
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
    ["Sample", "Exposures", "Kept", "Tags", "Status"].forEach((label) =>
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
    expect(headers).toEqual(["Sample", "Exposures", "Kept", "Tags", "Status"]);
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

describe("SheetTable horizontal scroll + sticky identity (WCAG 1.4.10)", () => {
  it("header rowgroup and rows region share ONE sheet-scroll container", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const scroll = screen.getByTestId("sheet-scroll");
    expect(within(scroll).getByTestId("sheet-head")).toBeInTheDocument();
    expect(within(scroll).getByTestId("sheet-rows")).toBeInTheDocument();
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
