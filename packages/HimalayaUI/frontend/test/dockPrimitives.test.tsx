// Net-new branches from ledger 4b: KbKey `frost` variant.
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { KbKey } from "../src/print/ui/KbKey";

describe("KbKey frost variant", () => {
  it("plate is the default; frost swaps to translucent currentColor fill", () => {
    const { rerender } = render(<KbKey>X</KbKey>);
    expect(screen.getByTestId("kbkey").className).toContain("bg-plate");
    rerender(<KbKey variant="frost">↵</KbKey>);
    const k = screen.getByTestId("kbkey");
    expect(k.className).toContain("bg-current/15");
    expect(k.className).toContain("text-current");
    expect(k.className).not.toContain("bg-plate");
  });
});
