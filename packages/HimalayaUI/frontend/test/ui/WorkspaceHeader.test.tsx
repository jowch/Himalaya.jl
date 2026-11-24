import { describe, expect, it, vi } from "vitest";
import userEvent from "@testing-library/user-event";
import { render, screen, within } from "@testing-library/react";
import { MemoryRouter } from "react-router";
import { WorkspaceHeader } from "../../src/ui/WorkspaceHeader";

const onActions = vi.fn();
const header = () => (
  <MemoryRouter>
    <WorkspaceHeader experimentName="Run 1" onActions={onActions} onHelp={() => {}} />
  </MemoryRouter>
);

describe("WorkspaceHeader", () => {
  it("opens the command palette from an Actions button that shows its shortcut", async () => {
    render(header());
    const actions = screen.getByRole("button", { name: /^Actions/ });
    expect(within(actions).getByText("K")).toBeInTheDocument();
    await userEvent.click(actions);
    expect(onActions).toHaveBeenCalledTimes(1);
  });

  it("carries the experiments page's wordmark and a breadcrumb back to all experiments", () => {
    render(header());
    expect(screen.getByRole("link", { name: "Himalaya SAXS, go to experiments" })).toHaveAttribute("href", "/experiments");
    expect(screen.getByTestId("wordmark")).toHaveTextContent("Himalaya · SAXS");
    const crumbs = screen.getByRole("navigation", { name: "Breadcrumb" });
    expect(within(crumbs).getByRole("link", { name: "Experiments" })).toHaveAttribute("href", "/experiments");
    expect(within(crumbs).getByRole("heading", { level: 1, name: "Run 1" })).toHaveAttribute("aria-current", "page");
  });
});
