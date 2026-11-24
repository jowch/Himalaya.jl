import { describe, expect, it } from "vitest";
import { render, screen, within } from "@testing-library/react";
import { MemoryRouter } from "react-router";
import { WorkspaceHeader } from "../../src/ui/WorkspaceHeader";

const header = (reviewed: number, total: number) => (
  <MemoryRouter>
    <WorkspaceHeader experimentName="Run 1" reviewed={reviewed} total={total} onHelp={() => {}} />
  </MemoryRouter>
);

describe("WorkspaceHeader", () => {
  it("shows the review count", () => {
    render(header(2, 4));
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("2 of 4 reviewed");
  });

  it("shows no progress when no sample can be reviewed, but keeps its place hidden from sight and from assistive tech", () => {
    const { container } = render(header(0, 0));
    expect(screen.queryByTestId("workspace-progress")).toBeNull();
    const reserved = container.querySelector("[data-reserved]");
    expect(reserved).not.toBeNull();
    expect(reserved).toHaveAttribute("aria-hidden", "true");
  });

  it("carries the experiments page's wordmark and a breadcrumb back to all experiments", () => {
    render(header(2, 4));
    expect(screen.getByRole("link", { name: "Himalaya SAXS, go to experiments" })).toHaveAttribute("href", "/experiments");
    expect(screen.getByTestId("wordmark")).toHaveTextContent("Himalaya · SAXS");
    const crumbs = screen.getByRole("navigation", { name: "Breadcrumb" });
    expect(within(crumbs).getByRole("link", { name: "Experiments" })).toHaveAttribute("href", "/experiments");
    expect(within(crumbs).getByRole("heading", { level: 1, name: "Run 1" })).toHaveAttribute("aria-current", "page");
  });
});
