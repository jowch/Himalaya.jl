import { render, cleanup } from "@testing-library/react";
import { describe, it, expect, afterEach } from "vitest";
import { useDocumentTitle } from "../../src/hooks/useDocumentTitle";

function Titled({ lead }: { lead: string | null | undefined }): null {
  useDocumentTitle(lead);
  return null;
}

afterEach(() => {
  cleanup();
  document.title = "";
});

describe("useDocumentTitle", () => {
  it("suffixes the lead with the app name", () => {
    render(<Titled lead="LL37 ratio series" />);
    expect(document.title).toBe("LL37 ratio series · Himalaya");
  });

  it("falls back to the bare app name for a nullish or blank lead", () => {
    render(<Titled lead={null} />);
    expect(document.title).toBe("Himalaya");
    cleanup();
    render(<Titled lead="   " />);
    expect(document.title).toBe("Himalaya");
  });

  it("restores the previous title on unmount", () => {
    document.title = "Previous";
    const { unmount } = render(<Titled lead="Focus sample" />);
    expect(document.title).toBe("Focus sample · Himalaya");
    unmount();
    expect(document.title).toBe("Previous");
  });

  it("updates when the lead changes", () => {
    const { rerender } = render(<Titled lead="A" />);
    expect(document.title).toBe("A · Himalaya");
    rerender(<Titled lead="B" />);
    expect(document.title).toBe("B · Himalaya");
  });
});
