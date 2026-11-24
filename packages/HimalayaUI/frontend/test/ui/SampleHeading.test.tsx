import { describe, expect, it } from "vitest";
import { render, screen } from "@testing-library/react";
import { SampleHeading } from "../../src/ui/Section";

describe("SampleHeading", () => {
  it("puts the setup and the sample's place above its name, as the experiments page does", () => {
    render(<SampleHeading name="MO_PEG5" eyebrow="1.7 m · sample 3 of 13" />);
    const heading = screen.getByRole("heading", { level: 2, name: "MO_PEG5" });
    const eyebrow = screen.getByText("1.7 m · sample 3 of 13");
    expect(eyebrow.compareDocumentPosition(heading) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    expect(heading.className).toMatch(/(^| )font-semibold( |$)/);
  });
});
