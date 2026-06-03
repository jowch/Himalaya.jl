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
    expect(screen.getByText("This exposure")).toBeInTheDocument();
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

  it("wires set-representative", async () => {
    const onSetRepresentative = vi.fn();
    setup({ onSetRepresentative });
    await userEvent.click(screen.getByRole("button", { name: /set as representative/i }));
    expect(onSetRepresentative).toHaveBeenCalledOnce();
  });

  it("shows the default loupe keys", () => {
    setup();
    expect(screen.getByText("flip frames")).toBeInTheDocument();
    expect(screen.getByText("set representative")).toBeInTheDocument();
  });
});
