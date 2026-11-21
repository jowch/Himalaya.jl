import { describe, it, expect, afterEach } from "vitest";
import { act, render } from "@testing-library/react";
import { useDocumentTheme, useAppTheme } from "../src/lib/theme";

const KEY = "himalaya-ui:theme";
let toggle: () => void = () => {};

function App({ withWorkspace }: { withWorkspace: boolean }): JSX.Element {
  useDocumentTheme();
  return withWorkspace ? <Workspace /> : <p>home</p>;
}

function Workspace(): JSX.Element {
  const theme = useAppTheme();
  toggle = theme.toggle;
  return <div data-testid="ws" data-theme={theme.theme} />;
}

function pickInAnotherTab(value: "light" | "dark"): void {
  window.localStorage.setItem(KEY, value);
  window.dispatchEvent(new StorageEvent("storage", { key: KEY, newValue: value }));
}

describe("the theme picked in the workspace", () => {
  afterEach(() => {
    window.localStorage.removeItem(KEY);
    window.localStorage.removeItem("himalaya-ui:workspace-theme");
    document.documentElement.removeAttribute("data-theme");
    document.documentElement.classList.remove("dark");
  });

  it("applies to the whole app, so pages outside the workspace follow it", () => {
    window.localStorage.setItem(KEY, "dark");
    const { rerender } = render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-theme", "dark");
    expect(document.documentElement).toHaveClass("dark");
    rerender(<App withWorkspace />);
    act(() => toggle());
    expect(document.documentElement).toHaveAttribute("data-theme", "light");
    expect(document.documentElement).not.toHaveClass("dark");
  });

  it("follows the system setting when nothing was picked", () => {
    render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-theme", "system");
  });

  it("keeps the page and the workspace toggle in step when another tab picks a theme", () => {
    window.localStorage.setItem(KEY, "light");
    const { getByTestId } = render(<App withWorkspace />);
    act(() => pickInAnotherTab("dark"));
    expect(document.documentElement).toHaveAttribute("data-theme", "dark");
    expect(getByTestId("ws")).toHaveAttribute("data-theme", "dark");
    act(() => toggle());
    expect(document.documentElement).toHaveAttribute("data-theme", "light");
    expect(getByTestId("ws")).toHaveAttribute("data-theme", "light");
  });

  it("keeps a pick saved under the old workspace-only key, and moves it to the new key", () => {
    window.localStorage.setItem("himalaya-ui:workspace-theme", "dark");
    render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-theme", "dark");
    expect(window.localStorage.getItem(KEY)).toBe("dark");
    expect(window.localStorage.getItem("himalaya-ui:workspace-theme")).toBeNull();
  });
});
