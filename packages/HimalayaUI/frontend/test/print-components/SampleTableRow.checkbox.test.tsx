import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";

const BASE_PROPS = {
  name: "DOPC-1:1",
  sampleId: "#10",
  exposures: [],
  kept: 2,
  total: 3,
  tags: [],
};

describe("SampleTableRow — checkbox column", () => {
  it("renders no checkbox cell when onCheck is absent (backward compat)", () => {
    render(<SampleTableRow {...BASE_PROPS} />);
    expect(screen.queryByRole("checkbox")).toBeNull();
  });

  it("renders a checkbox cell when onCheck is provided", () => {
    render(<SampleTableRow {...BASE_PROPS} onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toBeInTheDocument();
  });

  it("checkbox has data-checked=false when unchecked", () => {
    render(<SampleTableRow {...BASE_PROPS} checked={false} onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toHaveAttribute("data-checked", "false");
  });

  it("checkbox has data-checked=true when checked", () => {
    render(<SampleTableRow {...BASE_PROPS} checked onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toHaveAttribute("data-checked", "true");
  });

  it("calls onCheck when checkbox is clicked", async () => {
    const user = userEvent.setup();
    const onCheck = vi.fn();
    render(<SampleTableRow {...BASE_PROPS} onCheck={onCheck} />);
    await user.click(screen.getByRole("checkbox"));
    expect(onCheck).toHaveBeenCalledTimes(1);
  });
});
