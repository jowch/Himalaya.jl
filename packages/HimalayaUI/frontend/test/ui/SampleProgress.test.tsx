import { describe, expect, it } from "vitest";
import { render, screen } from "@testing-library/react";
import { MAX_SEPARATE_SEGMENTS, SampleProgress, ribbonRuns, type ReviewMark } from "../../src/ui/SampleProgress";

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

  it("folds a long run of samples into one segment per stretch of the same state, so hundreds of samples never smear", () => {
    const marks: ReviewMark[] = Array.from({ length: MAX_SEPARATE_SEGMENTS + 40 }, (_, i) => (i < 30 ? "reviewed" : i < 90 ? "open" : "none"));
    expect(ribbonRuns(marks)).toEqual([{ mark: "reviewed", count: 30 }, { mark: "open", count: 60 }, { mark: "none", count: 10 }]);
    render(<SampleProgress marks={marks} />);
    const ribbon = screen.getByTestId("sample-ribbon");
    expect(ribbon).not.toHaveAttribute("data-separate");
    expect([...ribbon.querySelectorAll("[data-mark]")].map((n) => n.getAttribute("data-count"))).toEqual(["30", "60", "10"]);
  });
});
