import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PrintApp } from "../src/print/App";

describe("PrintApp shell", () => {
  it("renders the greenfield shell landmark", () => {
    render(<PrintApp />);
    expect(screen.getByTestId("print-shell")).toBeInTheDocument();
  });
});
