import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { MetaList } from "../../src/ui/MetaList";
import { IMAGE_INFO_TERMS } from "../../src/ui/explanations";

describe("MetaList", () => {
  it("explains a label that has an explanation, and leaves the others plain", () => {
    render(<MetaList rows={[["I₂", "1,234"], ["File", "a.tif"]]} explain={IMAGE_INFO_TERMS} />);
    expect(screen.getByText("I₂").closest("[tabindex]")).toHaveAccessibleDescription(IMAGE_INFO_TERMS["I₂"]!);
    expect(screen.getByText("File").closest("[tabindex]")).toBeNull();
  });

  it("explains every physics row of the image info", () => {
    expect(Object.keys(IMAGE_INFO_TERMS).sort()).toEqual(["Energy", "Exposure", "Flight tube", "I₀", "I₂", "Kind"].sort());
  });
});
