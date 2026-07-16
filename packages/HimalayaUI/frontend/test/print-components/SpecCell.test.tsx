import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SpecCell } from "../../src/print/components/SpecCell";

describe("<SpecCell> onOpenLoupe", () => {
  it("renders the name as a button that fires onOpenLoupe when provided", () => {
    const onOpenLoupe = vi.fn();
    render(<SpecCell name="POPC" sampleId="#42" onOpenLoupe={onOpenLoupe} />);
    fireEvent.click(screen.getByRole("button", { name: /POPC/ }));
    expect(onOpenLoupe).toHaveBeenCalled();
  });

  it("renders the name as static text (no button) when onOpenLoupe is absent", () => {
    render(<SpecCell name="POPC" sampleId="#42" />);
    expect(screen.queryByRole("button", { name: /POPC/ })).toBeNull();
  });
});
