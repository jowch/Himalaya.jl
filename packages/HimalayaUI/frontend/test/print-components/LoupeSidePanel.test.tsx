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
    expect(screen.getByTestId("meta-list")).toBeInTheDocument();
    expect(screen.getByTestId("tag-list")).toBeInTheDocument();
    // The "Keys" legend was removed — redundant with the dock.
    expect(screen.queryByText("Keys")).toBeNull();
    expect(screen.queryByTestId("kb-legend")).toBeNull();
  });

  it("suppresses the 'This frame' block when the only row is the redundant position (LO-EXPSPARSE)", () => {
    // The position ("N of M") is already shown in the subtitle + BigFrame
    // caption, so a lone-row block reads as an unfinished section.
    setup({ meta: [{ key: "frame", value: "1 of 1" }] });
    expect(screen.queryByText("This frame")).toBeNull();
    expect(screen.queryByTestId("meta-list")).toBeNull();
    // The rest of the panel is unaffected.
    expect(screen.getByText("Sample tags")).toBeInTheDocument();
  });

  it("shows the 'This frame' block once a second, real row exists (e.g. a rejection reason)", () => {
    setup({
      meta: [
        { key: "frame", value: "2 of 4" },
        { key: "reason", value: "beam dropout" },
      ],
    });
    expect(screen.getByText("This frame")).toBeInTheDocument();
    expect(screen.getByTestId("meta-list")).toBeInTheDocument();
    expect(screen.getByText("beam dropout")).toBeInTheDocument();
  });

  it("reflects a kept verdict and wires the single Drop toggle", async () => {
    const onToggleDrop = vi.fn();
    setup({ onToggleDrop });
    // dropped=false → the Verdict reads "Kept".
    expect(screen.getByText("Kept")).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: /keep/i })).toBeNull();
    await userEvent.click(screen.getByRole("button", { name: "Drop" }));
    expect(onToggleDrop).toHaveBeenCalledOnce();
  });

  it("a dropped verdict still offers the Drop toggle (press again to bring back)", async () => {
    const onToggleDrop = vi.fn();
    setup({ dropped: true, onToggleDrop });
    expect(screen.getByText("Dropped")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: "Drop" }));
    expect(onToggleDrop).toHaveBeenCalledOnce();
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
