import { describe, expect, it } from "vitest";
import { render, screen } from "@testing-library/react";
import { SampleProgress } from "../../src/ui/SampleProgress";

describe("SampleProgress", () => {
  it("counts the reviewed samples and draws one ribbon segment per sample in list order", () => {
    const { container } = render(<SampleProgress marks={["reviewed", "open", "none", "reviewed"]} />);
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("2 of 3 reviewed");
    expect(screen.getByRole("img", { name: "2 of 3 samples reviewed" })).toBeInTheDocument();
    const segments = [...container.querySelectorAll("[data-mark]")].map((n) => n.getAttribute("data-mark"));
    expect(segments).toEqual(["reviewed", "open", "none", "reviewed"]);
  });

  it("shows no count when no sample can be reviewed, but keeps its place hidden from sight and from assistive tech", () => {
    const { container } = render(<SampleProgress marks={["none", "none"]} />);
    expect(screen.queryByTestId("workspace-progress")).toBeNull();
    const reserved = container.querySelector("[data-reserved]");
    expect(reserved).not.toBeNull();
    expect(reserved).toHaveAttribute("aria-hidden", "true");
  });
});
