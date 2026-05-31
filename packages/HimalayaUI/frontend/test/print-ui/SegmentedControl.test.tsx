import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SegmentedControl } from "../../src/print/ui/SegmentedControl";

const opts = [
  { value: "comb", label: "Comb" },
  { value: "resid", label: "Resid" },
] as const;

describe("<SegmentedControl> xs size", () => {
  it("reflects size=xs on data-size", () => {
    render(
      <SegmentedControl aria-label="view" options={opts} value="comb" onChange={() => {}} size="xs" />,
    );
    expect(screen.getByRole("group").getAttribute("data-size")).toBe("xs");
  });
  it("still toggles via onChange at xs", () => {
    const onChange = vi.fn();
    render(
      <SegmentedControl aria-label="view" options={opts} value="comb" onChange={onChange} size="xs" />,
    );
    const btns = screen.getAllByRole("button");
    btns[1].click();
    expect(onChange).toHaveBeenCalledWith("resid");
  });
});
