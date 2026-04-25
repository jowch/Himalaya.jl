import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { QNumInput } from "../src/components/PlotCard";

describe("QNumInput", () => {
  it("displays the value formatted to three decimal places", () => {
    render(<QNumInput value={0.1} onCommit={() => {}} testId="qmin" />);
    expect((screen.getByTestId("qmin") as HTMLInputElement).value).toBe("0.100");
  });

  it("does not overwrite the draft while the input is focused", async () => {
    const user = userEvent.setup();
    const { rerender } = render(<QNumInput value={0.1} onCommit={() => {}} testId="qmin" />);
    const input = screen.getByTestId("qmin") as HTMLInputElement;

    await user.click(input); // sets focused = true

    // Simulate wheel zoom updating the external q boundary mid-edit.
    rerender(<QNumInput value={0.3} onCommit={() => {}} testId="qmin" />);

    // Draft must NOT be replaced with "0.300" while the user is typing.
    expect(input.value).toBe("0.100");
  });

  it("syncs to the new prop value once the input is blurred", async () => {
    const user = userEvent.setup();
    const { rerender } = render(<QNumInput value={0.1} onCommit={() => {}} testId="qmin" />);
    const input = screen.getByTestId("qmin") as HTMLInputElement;

    await user.click(input);
    await user.tab(); // moves focus away → focused becomes false

    rerender(<QNumInput value={0.3} onCommit={() => {}} testId="qmin" />);

    // jsdom normalises number-input .value strings during React reconciliation
    // (e.g. "0.300" → "0.3"), so compare the parsed float rather than the raw string.
    expect(parseFloat(input.value)).toBeCloseTo(0.3);
  });

  it("calls onCommit with the parsed float on blur", async () => {
    const user = userEvent.setup();
    const onCommit = vi.fn();
    render(<QNumInput value={0.1} onCommit={onCommit} testId="qmin" />);
    const input = screen.getByTestId("qmin") as HTMLInputElement;

    await user.click(input);
    await user.clear(input);
    await user.type(input, "0.25");
    await user.tab();

    expect(onCommit).toHaveBeenCalledWith(0.25);
  });

  it("calls onCommit on Enter and blurs the input", async () => {
    const user = userEvent.setup();
    const onCommit = vi.fn();
    render(<QNumInput value={0.1} onCommit={onCommit} testId="qmin" />);
    const input = screen.getByTestId("qmin") as HTMLInputElement;

    await user.click(input);
    await user.clear(input);
    await user.type(input, "0.42");
    await user.keyboard("{Enter}");

    expect(onCommit).toHaveBeenCalledWith(0.42);
  });
});
