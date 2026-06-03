import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ReadingPanel, type ReadingDatum } from "../../src/print/components/ReadingPanel";

const readings: ReadingDatum[] = [
  { phase: "Pn3m", span: "1:0 → 1:0.5", lattice: "a 205 → 195 Å" },
  { phase: "Lamellar", span: "1:0.5 → 1:1", lattice: "d 60 → 55 Å" },
];

describe("ReadingPanel", () => {
  it("renders one ReadingRow per reading", () => {
    const { getAllByTestId } = render(<ReadingPanel readings={readings} />);
    expect(getAllByTestId("reading-row")).toHaveLength(2);
  });
  it("renders coexistence and form-factor notes when provided", () => {
    const { getByText } = render(
      <ReadingPanel readings={readings} coexistenceNote="coexistence at 1:0.5" formFactorNote="form factor only at 1:1.5" />,
    );
    expect(getByText("coexistence at 1:0.5")).toBeTruthy();
    expect(getByText("form factor only at 1:1.5")).toBeTruthy();
  });
  it("omits notes when not provided", () => {
    const { queryByTestId } = render(<ReadingPanel readings={readings} />);
    expect(queryByTestId("reading-coex")).toBeNull();
  });
});
