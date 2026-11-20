import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SpecCell } from "../../src/print/components/SpecCell";

describe("<SpecCell> onOpenSample", () => {
  it("renders the name as a button that fires onOpenSample when provided", () => {
    const onOpenSample = vi.fn();
    render(<SpecCell name="POPC" sampleId="#42" onOpenSample={onOpenSample} />);
    fireEvent.click(screen.getByRole("button", { name: /POPC/ }));
    expect(onOpenSample).toHaveBeenCalled();
  });

  it("renders the name as static text (no button) when onOpenSample is absent", () => {
    render(<SpecCell name="POPC" sampleId="#42" />);
    expect(screen.queryByRole("button", { name: /POPC/ })).toBeNull();
  });
});
