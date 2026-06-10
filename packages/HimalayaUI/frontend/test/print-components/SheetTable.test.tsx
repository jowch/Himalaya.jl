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
});
