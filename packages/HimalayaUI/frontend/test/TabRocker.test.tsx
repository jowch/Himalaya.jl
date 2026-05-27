import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter } from "react-router-dom";
import { TabRocker } from "../src/components/TabRocker";
import { AppHeader } from "../src/components/AppHeader";
import { useAppState } from "../src/state";

// Wrap in MemoryRouter so the component has a routing context.
const renderInRouter = (ui: React.ReactElement) =>
  render(<MemoryRouter>{ui}</MemoryRouter>);

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activePage: "none",
    username: "tester",
    theme: "dark",
  });
});

describe("<TabRocker>", () => {
  // I3.6 (#177): Compare retired (Index #181, Inspect #163 before it). No legacy
  // page tabs remain, so the rocker renders an empty tablist. The component +
  // the whole dual-nav model are deleted in I5.1.
  it("renders no tabs (all legacy surfaces retired)", () => {
    renderInRouter(<TabRocker />);
    expect(screen.getByTestId("tab-rocker")).toBeInTheDocument();
    expect(screen.queryAllByRole("tab")).toHaveLength(0);
    expect(screen.queryByTestId("tab-compare")).toBeNull();
    expect(screen.queryByTestId("tab-index")).toBeNull();
    expect(screen.queryByTestId("tab-inspect")).toBeNull();
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
