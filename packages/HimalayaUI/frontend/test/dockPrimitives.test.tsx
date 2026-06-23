// Net-new branches from ledger 4b/4c: KbKey `frost` variant + CullBar `actions`
// override. Trivial conditional rendering, one assertion each.
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { KbKey } from "../src/print/ui/KbKey";
import { CullBar } from "../src/print/components/CullBar";

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

describe("CullBar actions override", () => {
  it("renders the supplied actions in place of Drop/Keep/Restore/Clear", () => {
    const onClick = vi.fn();
    render(
      <CullBar
        count={2}
        show
        actions={[{ label: "Merge", onClick, variant: "accent" }]}
      />,
    );
    expect(screen.getByRole("button", { name: "Merge" })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: /drop/i })).not.toBeInTheDocument();
    screen.getByRole("button", { name: "Merge" }).click();
    expect(onClick).toHaveBeenCalledOnce();
  });
});
