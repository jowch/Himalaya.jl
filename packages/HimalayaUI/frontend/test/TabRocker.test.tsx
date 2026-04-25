import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { TabRocker } from "../src/components/TabRocker";
import { AppHeader } from "../src/components/AppHeader";
import { useAppState } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activePage: "index",
    username: "tester",
    theme: "dark",
  });
});

describe("<TabRocker>", () => {
  it("renders both tabs with the active page marked", () => {
    render(<TabRocker />);
    const idx = screen.getByTestId("tab-index");
    const cmp = screen.getByTestId("tab-compare");
    expect(idx).toHaveAttribute("aria-selected", "true");
    expect(idx).toHaveAttribute("data-active", "true");
    expect(cmp).toHaveAttribute("aria-selected", "false");
  });

  it("clicking a tab updates activePage in the store", async () => {
    const user = userEvent.setup();
    render(<TabRocker />);
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
