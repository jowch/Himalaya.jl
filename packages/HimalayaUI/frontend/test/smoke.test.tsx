import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { App } from "../src/App";

describe("App smoke", () => {
  it("renders the loading text", () => {
    render(<App />);
    expect(screen.getByText(/loading/i)).toBeInTheDocument();
  });
});
