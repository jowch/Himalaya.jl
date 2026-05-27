import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ConflictModalShell, type ConflictPanelData } from "../src/components/ConflictModalShell";

const server: ConflictPanelData = {
  label: "Server (current)", testId: "conflict-panel-server",
  title: "Server", memberCount: 2, description: "srv", updatedAt: "2026-05-02",
};
const local: ConflictPanelData = {
  label: "Your draft", testId: "conflict-panel-local",
  title: "Mine", memberCount: 1, description: null, updatedAt: null,
};

function shell(over: Partial<React.ComponentProps<typeof ConflictModalShell>> = {}) {
  const props: React.ComponentProps<typeof ConflictModalShell> = {
    open: true,
    heading: "Heading", subtitle: "Subtitle",
    serverPanel: server, localPanel: local,
    onClose: vi.fn(), onDiscard: vi.fn(), discardLabel: "Discard",
    onOverwrite: vi.fn(), overwriteBusy: false,
    ...over,
  };
  return { props, ...render(<ConflictModalShell {...props} />) };
}

describe("ConflictModalShell", () => {
  it("renders nothing when closed", () => {
    const { container } = shell({ open: false });
    expect(container.firstChild).toBeNull();
  });

  it("renders both panels + heading/subtitle", () => {
    shell();
    expect(screen.getByText("Heading")).toBeInTheDocument();
    expect(screen.getByText("Subtitle")).toBeInTheDocument();
    expect(screen.getByTestId("conflict-panel-server-title")).toHaveTextContent("Server");
    expect(screen.getByTestId("conflict-panel-local-title")).toHaveTextContent("Mine");
  });

  it("Esc calls onClose without committing", () => {
    const onClose = vi.fn();
    shell({ onClose });
    fireEvent.keyDown(document, { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("Discard + Overwrite forward to their handlers; busy disables overwrite", () => {
    const onDiscard = vi.fn();
    const onOverwrite = vi.fn();
    const { rerender, props } = shell({ onDiscard, onOverwrite });
    fireEvent.click(screen.getByTestId("conflict-discard"));
    expect(onDiscard).toHaveBeenCalled();
    fireEvent.click(screen.getByTestId("conflict-overwrite"));
    expect(onOverwrite).toHaveBeenCalled();
    rerender(<ConflictModalShell {...props} overwriteBusy />);
    expect(screen.getByTestId("conflict-overwrite")).toBeDisabled();
    expect(screen.getByTestId("conflict-overwrite")).toHaveTextContent("Saving…");
  });

  it("renders an extraAction slot when provided", () => {
    shell({ extraAction: <button data-testid="conflict-fork">Fork</button> });
    expect(screen.getByTestId("conflict-fork")).toBeInTheDocument();
  });
});
