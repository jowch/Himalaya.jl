import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SeriesMemberRow } from "../../src/print/components/SeriesMemberRow";

describe("SeriesMemberRow", () => {
  it("renders a single-phase member: solid swatch, one colored name, variable + data line", () => {
    const { getByTestId, getByText, container } = render(
      <SeriesMemberRow phases={["Pn3m"]} variableValue="1:0" dataLine="a = 205 Å · q₁ 0.061 Å⁻¹" />,
    );
    expect(getByTestId("series-member-row")).toBeTruthy();
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("1:0")).toBeTruthy();
    expect(getByText("a = 205 Å · q₁ 0.061 Å⁻¹")).toBeTruthy();
    const sw = container.querySelector("[data-swatch]")!;
    expect(sw.getAttribute("data-coexist")).toBeNull();
    expect(sw.getAttribute("data-empty")).toBeNull();
  });

  it("renders a coexistence member: gradient swatch + both names", () => {
    const { getByText, container } = render(
      <SeriesMemberRow phases={["Pn3m", "Lamellar"]} variableValue="1:0.5" dataLine="a 195 · d 60 Å · q₁ 0.057 Å⁻¹" />,
    );
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("Lamellar")).toBeTruthy();
    expect(container.querySelector("[data-swatch]")!.getAttribute("data-coexist")).toBe("Lamellar");
  });

  it("renders a form-factor member: empty swatch + 'Form factor'", () => {
    const { getByText, container } = render(
      <SeriesMemberRow phases={[]} variableValue="1:1.5" dataLine="no Bragg peaks · q₁ —" />,
    );
    expect(getByText("Form factor")).toBeTruthy();
    expect(container.querySelector("[data-swatch]")!.getAttribute("data-empty")).toBe("true");
  });

  it("marks hot state and fires hover handlers", () => {
    const onHover = vi.fn(), onLeave = vi.fn();
    const { getByTestId } = render(
      <SeriesMemberRow phases={["Pn3m"]} variableValue="1:0" dataLine="x" hot onHover={onHover} onLeave={onLeave} />,
    );
    const row = getByTestId("series-member-row");
    expect(row.getAttribute("data-hot")).toBe("true");
    fireEvent.mouseEnter(row); fireEvent.mouseLeave(row);
    expect(onHover).toHaveBeenCalledTimes(1);
    expect(onLeave).toHaveBeenCalledTimes(1);
  });
});
