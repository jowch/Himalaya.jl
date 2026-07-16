import { test, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { KeptCell } from "../src/print/components/KeptCell";

test("KeptCell renders a plain kept / total count", () => {
  render(<KeptCell kept={2} total={3} />);
  expect(screen.getByTestId("kept-cell")).toHaveTextContent("2 / 3");
});

test("KeptCell has no Restore button (binary model: keep un-drops the frame)", () => {
  render(<KeptCell kept={2} total={3} />);
  expect(screen.queryByRole("button", { name: /restore/i })).toBeNull();
});
