import { describe, it, expect, afterEach } from "vitest";
import { act, render } from "@testing-library/react";
import { useDocumentTheme, useWorkspaceTheme } from "../../src/lib/theme";

const KEY = "himalaya-ui:workspace-theme";
let toggle: () => void = () => {};

function App({ withWorkspace }: { withWorkspace: boolean }): JSX.Element {
  useDocumentTheme();
  return withWorkspace ? <Workspace /> : <p>home</p>;
}

function Workspace(): JSX.Element {
  const theme = useWorkspaceTheme();
  toggle = theme.toggle;
  return <div data-testid="ws" data-workspace-theme={theme.theme} />;
}

function pickInAnotherTab(value: "light" | "dark"): void {
  window.localStorage.setItem(KEY, value);
  window.dispatchEvent(new StorageEvent("storage", { key: KEY, newValue: value }));
}

describe("the theme picked in the workspace", () => {
  afterEach(() => {
    window.localStorage.removeItem(KEY);
    document.documentElement.removeAttribute("data-workspace-theme");
    document.documentElement.classList.remove("dark");
  });

  it("applies to the whole app, so pages outside the workspace follow it", () => {
    window.localStorage.setItem(KEY, "dark");
    const { rerender } = render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "dark");
    expect(document.documentElement).toHaveClass("dark");
    rerender(<App withWorkspace />);
    act(() => toggle());
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "light");
    expect(document.documentElement).not.toHaveClass("dark");
  });

  it("follows the system setting when nothing was picked", () => {
    render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "system");
  });

  it("keeps the page and the workspace toggle in step when another tab picks a theme", () => {
    window.localStorage.setItem(KEY, "light");
    const { getByTestId } = render(<App withWorkspace />);
    act(() => pickInAnotherTab("dark"));
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "dark");
    expect(getByTestId("ws")).toHaveAttribute("data-workspace-theme", "dark");
    act(() => toggle());
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "light");
    expect(getByTestId("ws")).toHaveAttribute("data-workspace-theme", "light");
  });
});
