// test/PageFrame.ingestion.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PageFrame, PAGE_WIDTHS } from "../src/print/components/PageFrame";

describe("PageFrame ingestion widths (Phase E1)", () => {
  it("exposes home + experiment width keys", () => {
    expect(PAGE_WIDTHS.home).toBeDefined();
    expect(PAGE_WIDTHS.experiment).toBeDefined();
  });
  it("renders with the new keys", () => {
    render(<PageFrame width="experiment">x</PageFrame>);
    expect(screen.getByTestId("page-frame")).toBeInTheDocument();
  });
});
