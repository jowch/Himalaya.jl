import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { SheetTable } from "../../src/print/components/SheetTable";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";

const row = (id: string) => (
  <SampleTableRow key={id} name={`Sample ${id}`} sampleId={id} exposures={[]} kept={1} total={1} tags={[]} />
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
