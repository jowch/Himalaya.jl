import { describe, it, expect, afterEach } from "vitest";
import { act, render } from "@testing-library/react";
import { useDocumentTheme, useWorkspaceTheme } from "../../src/workspace/theme";

let toggle: () => void = () => {};

function App({ withWorkspace }: { withWorkspace: boolean }): JSX.Element {
  useDocumentTheme();
  return withWorkspace ? <Workspace /> : <p>home</p>;
}

function Workspace(): null {
  toggle = useWorkspaceTheme().toggle;
  return null;
}

describe("the theme picked in the workspace", () => {
  afterEach(() => {
    window.localStorage.removeItem("himalaya-ui:workspace-theme");
    document.documentElement.removeAttribute("data-workspace-theme");
  });

  it("applies to the whole app, so pages outside the workspace follow it", () => {
    window.localStorage.setItem("himalaya-ui:workspace-theme", "dark");
    const { rerender } = render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "dark");
    rerender(<App withWorkspace />);
    act(() => toggle());
    rerender(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "light");
  });

  it("follows the system setting when nothing was picked", () => {
    render(<App withWorkspace={false} />);
    expect(document.documentElement).toHaveAttribute("data-workspace-theme", "system");
  });
});
