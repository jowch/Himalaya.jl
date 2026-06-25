// test/TopNav.test.tsx
import { describe, test, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { TopNav } from "../src/print/shell/TopNav";

function renderWithRouter(ui: React.ReactElement, { route = "/" }: { route?: string } = {}) {
  render(<MemoryRouter initialEntries={[route]}>{ui}</MemoryRouter>);
}

describe("TopNav", () => {
  test("TopNav: wordmark→/experiments, two tabs, no Samples/beamtime/gear", () => {
    renderWithRouter(<TopNav />, { route: "/experiments" });
    expect(screen.getByRole("link", { name: /himalaya/i })).toHaveAttribute("href", "/experiments");
    expect(screen.getByRole("link", { name: "Experiments" })).toBeInTheDocument();
    expect(screen.getByRole("link", { name: "Series" })).toBeInTheDocument();
    expect(screen.queryByText("Samples")).toBeNull();
    expect(screen.queryByRole("button", { name: /beamtime/i })).toBeNull();
  });
});
