import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { MetaList } from "../../src/print/ui/MetaList";
import type { MetaEntry } from "../../src/print/ui/MetaList";

const entries: MetaEntry[] = [
  { key: "frame", value: "037" },
  { key: "integration", value: "1.2 s" },
  { key: "collected", value: "2024-03-14" },
];

describe("MetaList", () => {
  it("renders each key and value text", () => {
    render(<MetaList entries={entries} />);
    for (const e of entries) {
      expect(screen.getByText(e.key)).toBeTruthy();
      expect(screen.getByText(e.value as string)).toBeTruthy();
    }
  });

  it("uses semantic <dl>/<dt>/<dd>: term + definition per entry", () => {
    render(<MetaList entries={entries} />);
    // In the a11y tree, <dt> exposes role "term" and <dd> role "definition".
    expect(screen.getAllByRole("term")).toHaveLength(entries.length);
    expect(screen.getAllByRole("definition")).toHaveLength(entries.length);
  });

  it('has data-testid="meta-list"', () => {
    render(<MetaList entries={entries} />);
    expect(screen.getByTestId("meta-list")).toBeTruthy();
  });

  it("forwards a placement className", () => {
    render(<MetaList entries={entries} className="mt-4" />);
    expect(screen.getByTestId("meta-list").className).toContain("mt-4");
  });
});
