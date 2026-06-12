import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { LoupeSidePanel } from "../../src/print/components/LoupeSidePanel";

const meta = [
  { key: "frame", value: "65" },
  { key: "exposure", value: "0.40 s" },
];
const tags = [{ key: "LL37" }, { key: "temp", value: "37C" }];

function setup(over = {}) {
  const props = { meta, dropped: false, isRepresentative: false, tags, ...over };
  render(<LoupeSidePanel {...props} />);
  return props;
}

describe("LoupeSidePanel", () => {
  it("renders the section headers + composed leaves", () => {
    setup();
    expect(screen.getByTestId("loupe-side-panel")).toBeInTheDocument();
    // LO-TERM: the loupe speaks ONE word for a single detector image — "frame"
    // (matching the BigFrame caption, toasts, and meta keys); not "exposure".
    expect(screen.getByText("This frame")).toBeInTheDocument();
    expect(screen.queryByText("This exposure")).toBeNull();
    expect(screen.getByText("Sample tags")).toBeInTheDocument();
    expect(screen.getByText("Keys")).toBeInTheDocument();
    expect(screen.getByTestId("meta-list")).toBeInTheDocument();
    expect(screen.getByTestId("tag-list")).toBeInTheDocument();
    expect(screen.getByTestId("kb-legend")).toBeInTheDocument();
  });

  it("reflects kept verdict and wires the drop toggle", async () => {
    const onToggleDrop = vi.fn();
    setup({ onToggleDrop });
    await userEvent.click(screen.getByRole("button", { name: /drop/i }));
    expect(onToggleDrop).toHaveBeenCalledOnce();
  });

  it("threads the keep state and wires the keep toggle (SA-SCREENED)", async () => {
    const onToggleKeep = vi.fn();
    setup({ kept: true, onToggleKeep });
    // kept=true → the Verdict reads "Kept" and the keep toggle reads Restore.
    expect(screen.getByText("Kept")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: "Restore" }));
    expect(onToggleKeep).toHaveBeenCalledOnce();
  });

  it("offers Keep on an unscreened exposure", async () => {
    const onToggleKeep = vi.fn();
    setup({ onToggleKeep });
    expect(screen.getByText("Unscreened")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: "Keep" }));
    expect(onToggleKeep).toHaveBeenCalledOnce();
  });

  it("documents K alongside X in the default key legend", () => {
    setup();
    expect(screen.getByText("drop / restore")).toBeInTheDocument();
    expect(screen.getByText("keep / restore")).toBeInTheDocument();
  });

  it("wires set-representative", async () => {
    const onSetRepresentative = vi.fn();
    setup({ onSetRepresentative });
    await userEvent.click(screen.getByRole("button", { name: /set as representative/i }));
    expect(onSetRepresentative).toHaveBeenCalledOnce();
  });

  it("threads representativeDropped through to the RepresentativeBox warning", () => {
    setup({ representativeDropped: true });
    expect(screen.getByTestId("rep-dropped-warning")).toBeInTheDocument();
  });

  it("no dropped-representative warning by default", () => {
    setup();
    expect(screen.queryByTestId("rep-dropped-warning")).not.toBeInTheDocument();
  });

  it("shows the default loupe keys", () => {
    setup();
    expect(screen.getByText("flip frames")).toBeInTheDocument();
    expect(screen.getByText("set representative")).toBeInTheDocument();
  });

  it("shows the add-tag affordance persistently (the loupe rule) and commits a new tag", async () => {
    const onAddTag = vi.fn();
    setup({ onAddTag });
    // Loupe tags keep the dashed add chip always visible, even with tags present.
    expect(screen.getByTestId("tag-list")).toHaveAttribute(
      "data-persistent-add",
      "true",
    );
    await userEvent.click(screen.getByRole("button", { name: "Add" }));
    const editor = screen.getByTestId("tag-editor");
    const keyInput = editor.querySelector(
      'input[placeholder="key"]',
    ) as HTMLInputElement;
    await userEvent.type(keyInput, "pH");
    await userEvent.keyboard("{Enter}");
    expect(onAddTag).toHaveBeenCalledWith({ key: "pH" });
  });
});
