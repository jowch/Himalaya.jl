// test/ExperimentTopNav.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { ExperimentTopNav } from "../src/print/shell/ExperimentTopNav";

function at(path: string) {
  render(
    <MemoryRouter initialEntries={[path]}>
      <ExperimentTopNav />
    </MemoryRouter>,
  );
}

describe("ExperimentTopNav (Phase E1)", () => {
  it("links to /experiments and /series", () => {
    at("/experiments");
    expect(screen.getByTestId("topnav-experiments")).toHaveAttribute("href", "/experiments");
    expect(screen.getByTestId("topnav-series")).toHaveAttribute("href", "/series");
  });

  it("marks Experiments active on an experiment route", () => {
    at("/experiments/7/corpus");
    expect(screen.getByTestId("topnav-experiments")).toHaveAttribute("aria-current", "page");
    expect(screen.getByTestId("topnav-series")).not.toHaveAttribute("aria-current");
  });

  it("marks Series active on /series", () => {
    at("/series");
    expect(screen.getByTestId("topnav-series")).toHaveAttribute("aria-current", "page");
  });
});
