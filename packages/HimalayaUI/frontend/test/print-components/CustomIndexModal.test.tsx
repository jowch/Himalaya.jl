import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CustomIndexModal } from "../../src/print/components/CustomIndexModal";
import { PN3M } from "../../src/print/comb/comb.fixtures";

const base = {
  open: true,
  symmetries: ["Pn3m", "Im3m", "Ia3d", "Lamellar", "Hexagonal"],
  symmetry: "Im3m",
  paramName: "a",
  paramValue: "252",
  paramMin: 120, paramMax: 360, paramStep: 1,
  previewSeries: PN3M,
  observed: [0.03, 0.045, 0.06],
  fit: { landed: 3, total: 5, snapped: true },
};

describe("CustomIndexModal", () => {
  it("renders head, symmetry seg, lattice control, preview, fit line and footer when open", () => {
    const { getByTestId, getByText, getAllByRole } = render(
      <CustomIndexModal {...base} onClose={() => {}} onCancel={() => {}} onAdd={() => {}}
        onSymmetryChange={() => {}} onParamChange={() => {}} />,
    );
    expect(getByText("Custom index")).toBeTruthy();
    expect(getByText("Speculative")).toBeTruthy();
    expect(getAllByRole("button").some((b) => b.textContent === "Im3m")).toBe(true);
    expect(getByTestId("lattice-param")).toBeTruthy();
    expect(getByTestId("custom-preview")).toBeTruthy();
    expect(getByTestId("fit-metadata").textContent).toContain("Lands on");
    expect(getByText("Add to assignment")).toBeTruthy();
    expect(getByText("Cancel")).toBeTruthy();
  });
  it("renders nothing when open=false", () => {
    const { queryByText } = render(
      <CustomIndexModal {...base} open={false} onClose={() => {}} onCancel={() => {}} onAdd={() => {}}
        onSymmetryChange={() => {}} onParamChange={() => {}} />,
    );
    expect(queryByText("Custom index")).toBeNull();
  });
  it("wires the actions", () => {
    const onAdd = vi.fn(), onCancel = vi.fn(), onSymmetryChange = vi.fn();
    const { getByText } = render(
      <CustomIndexModal {...base} onClose={() => {}} onCancel={onCancel} onAdd={onAdd}
        onSymmetryChange={onSymmetryChange} onParamChange={() => {}} />,
    );
    fireEvent.click(getByText("Add to assignment"));
    fireEvent.click(getByText("Cancel"));
    fireEvent.click(getByText("Pn3m"));
    expect(onAdd).toHaveBeenCalledTimes(1);
    expect(onCancel).toHaveBeenCalledTimes(1);
    expect(onSymmetryChange).toHaveBeenCalledWith("Pn3m");
  });
});
