import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter } from "react-router-dom";
import { TabRocker } from "../src/components/TabRocker";
import { AppHeader } from "../src/components/AppHeader";
import { useAppState } from "../src/state";

// TabRocker uses `useNavigate()` (Plan §Phase 4: hybrid Zustand/router model
// — Compare nav goes through URL, Index/Inspect through Zustand). Wrap in
// MemoryRouter so tests have a routing context.
const renderInRouter = (ui: React.ReactElement) =>
  render(<MemoryRouter>{ui}</MemoryRouter>);

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activePage: "compare",
    username: "tester",
    theme: "dark",
  });
});

describe("<TabRocker>", () => {
  it("renders only the Compare tab (Index retired #181, Inspect retired #163)", () => {
    renderInRouter(<TabRocker />);
    const tabs = screen.getAllByRole("tab");
    expect(tabs).toHaveLength(1);
    expect(tabs[0]).toHaveTextContent("Compare");
    expect(screen.queryByTestId("tab-index")).toBeNull();
    expect(screen.queryByTestId("tab-inspect")).toBeNull();
  });

  it("marks the Compare tab active", () => {
    renderInRouter(<TabRocker />);
    const cmp = screen.getByTestId("tab-compare");
    expect(cmp).toHaveAttribute("aria-selected", "true");
    expect(cmp).toHaveAttribute("data-active", "true");
  });

  it("clicking the Compare tab keeps activePage 'compare'", async () => {
    const user = userEvent.setup();
    renderInRouter(<TabRocker />);
    await user.click(screen.getByTestId("tab-compare"));
    expect(useAppState.getState().activePage).toBe("compare");
  });
});

describe("<AppHeader>", () => {
  it("renders the utility cluster", () => {
    // The TabRocker now lives in AppShell (its own row below AppHeader),
    // not inside AppHeader itself.
    render(<AppHeader />);
    expect(screen.getByTestId("utility-cluster")).toBeInTheDocument();
  });

  it("theme toggle flips dark ↔ light", async () => {
    const user = userEvent.setup();
    render(<AppHeader />);
    expect(useAppState.getState().theme).toBe("dark");
    await user.click(screen.getByTestId("theme-toggle"));
    expect(useAppState.getState().theme).toBe("light");
    await user.click(screen.getByTestId("theme-toggle"));
    expect(useAppState.getState().theme).toBe("dark");
  });

  it("avatar click clears username (re-triggers onboarding)", async () => {
    const user = userEvent.setup();
    render(<AppHeader />);
    expect(useAppState.getState().username).toBe("tester");
    await user.click(screen.getByTestId("user-avatar"));
    expect(useAppState.getState().username).toBeUndefined();
  });

  it("avatar shows initials of username", () => {
    useAppState.setState({ username: "alice_smith" });
    render(<AppHeader />);
    expect(screen.getByTestId("user-avatar")).toHaveTextContent("AS");
  });
});
