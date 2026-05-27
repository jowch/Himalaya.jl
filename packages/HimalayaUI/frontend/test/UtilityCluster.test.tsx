import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { CorpusTopbar } from "../src/components/CorpusTopbar";
import { useAppState } from "../src/state";

// I5.1 (#182): the UtilityCluster (theme toggle + user-avatar) moved from the
// retired AppHeader/AppShell into CorpusTopbar — the corpus app's only surviving
// topbar and its only theme-toggle + user-switch affordance. These tests
// (salvaged from the deleted TabRocker.test.tsx <AppHeader> describe) render it
// through CorpusTopbar so they exercise the real production mount site.

function makeQc() {
  return new QueryClient({
    defaultOptions: { queries: { retry: false }, mutations: { retry: false } },
  });
}

const renderTopbar = () =>
  render(
    <QueryClientProvider client={makeQc()}>
      <MemoryRouter>
        <CorpusTopbar />
      </MemoryRouter>
    </QueryClientProvider>,
  );

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ username: "tester", theme: "dark" });
});

describe("UtilityCluster (in CorpusTopbar)", () => {
  it("renders the utility cluster in the corpus topbar", () => {
    renderTopbar();
    expect(screen.getByTestId("utility-cluster")).toBeInTheDocument();
  });

  it("theme toggle flips dark ↔ light", async () => {
    const user = userEvent.setup();
    renderTopbar();
    expect(useAppState.getState().theme).toBe("dark");
    await user.click(screen.getByTestId("theme-toggle"));
    expect(useAppState.getState().theme).toBe("light");
    await user.click(screen.getByTestId("theme-toggle"));
    expect(useAppState.getState().theme).toBe("dark");
  });

  it("avatar click clears username (re-triggers onboarding)", async () => {
    const user = userEvent.setup();
    renderTopbar();
    expect(useAppState.getState().username).toBe("tester");
    await user.click(screen.getByTestId("user-avatar"));
    expect(useAppState.getState().username).toBeUndefined();
  });

  it("avatar shows initials of username", () => {
    useAppState.setState({ username: "alice_smith" });
    renderTopbar();
    expect(screen.getByTestId("user-avatar")).toHaveTextContent("AS");
  });
});
