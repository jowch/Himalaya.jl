import { describe, it, expect } from "vitest";
import { prettifyUnits } from "../src/lib/units";

describe("prettifyUnits", () => {
  it("renders A-1 as Å⁻¹", () => {
    expect(prettifyUnits("A-1")).toBe("Å⁻¹");
  });

  it("renders nm-1 as nm⁻¹ without touching the lowercase letters", () => {
    expect(prettifyUnits("nm-1")).toBe("nm⁻¹");
  });

  it("renders 1/A as Å⁻¹", () => {
    expect(prettifyUnits("1/A")).toBe("Å⁻¹");
  });

  it("renders 1/nm as nm⁻¹", () => {
    expect(prettifyUnits("1/nm")).toBe("nm⁻¹");
  });

  it("is idempotent on already-pretty strings", () => {
    expect(prettifyUnits("Å⁻¹")).toBe("Å⁻¹");
    expect(prettifyUnits("nm⁻¹")).toBe("nm⁻¹");
  });

  it("supports higher inverse powers", () => {
    expect(prettifyUnits("nm-2")).toBe("nm⁻²");
    expect(prettifyUnits("A-3")).toBe("Å⁻³");
  });

  it("returns empty string unchanged", () => {
    expect(prettifyUnits("")).toBe("");
  });

  it("does not mangle a standalone unit symbol with no power", () => {
    // Without a superscript present, the rule for A → Å doesn't fire — that's
    // intentional, since "A" alone could mean Ampere; we only confidently
    // convert when it's clearly being used as an inverse-Angstrom unit.
    expect(prettifyUnits("A")).toBe("A");
    expect(prettifyUnits("nm")).toBe("nm");
  });
});
