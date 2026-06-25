import { describe, it, expect } from "vitest";
import { render, cleanup } from "@testing-library/react";
import { useInteraction } from "../../src/print/interaction/registry";
import { usePageActions } from "../../src/print/interaction/usePageActions";
import { core, page } from "../../src/print/interaction/core";

function Page({ label }: { label: string }): null {
  usePageActions({ actions: [core("openFocus", { run: () => {} }), page("cull", { label, keys: ["x"], group: "Act", run: () => {} })] });
  return null;
}

describe("useInteraction + usePageActions", () => {
  it("registers a page's actions on mount", () => {
    render(<Page label="Cull" />);
    const ids = useInteraction.getState().actions.map((a) => a.id);
    expect(ids).toContain("openFocus");
    expect(ids).toContain("cull");
  });

  it("clears the store on unmount", () => {
    const { unmount } = render(<Page label="Cull" />);
    expect(useInteraction.getState().actions.length).toBeGreaterThan(0);
    unmount();
    expect(useInteraction.getState().actions).toEqual([]);
    expect(useInteraction.getState().cursor).toBeNull();
  });

  it("throws when a page verb collides with a core key", () => {
    function Bad(): null {
      usePageActions({ actions: [page("cull", { label: "x", keys: ["Enter"], group: "Act", run: () => {} })] });
      return null;
    }
    expect(() => render(<Bad />)).toThrow(/Enter/);
    cleanup();
  });
});
