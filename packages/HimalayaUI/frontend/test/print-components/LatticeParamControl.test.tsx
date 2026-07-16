import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { LatticeParamControl } from "../../src/print/components/LatticeParamControl";

describe("LatticeParamControl", () => {
  it("renders a slider + a mono number field + unit, both reflecting value", () => {
    const { getByTestId, getByText } = render(
      <LatticeParamControl value="252" min={120} max={360} step={1} onValueChange={() => {}} unit="Å" />,
    );
    const slider = getByTestId("slider") as HTMLInputElement;
    expect(slider.value).toBe("252");
    const num = getByTestId("lattice-num").querySelector("input") as HTMLInputElement;
    expect(num.value).toBe("252");
    expect(num.className).toContain("font-mono");
    expect(getByText("Å")).toBeTruthy();
  });
  it("fires onValueChange from both controls", () => {
    const onValueChange = vi.fn();
    const { getByTestId } = render(
      <LatticeParamControl value="252" min={120} max={360} step={1} onValueChange={onValueChange} unit="Å" />,
    );
    fireEvent.change(getByTestId("slider"), { target: { value: "260" } });
    fireEvent.change(getByTestId("lattice-num").querySelector("input")!, { target: { value: "300" } });
    expect(onValueChange).toHaveBeenNthCalledWith(1, "260");
    expect(onValueChange).toHaveBeenNthCalledWith(2, "300");
  });
});
