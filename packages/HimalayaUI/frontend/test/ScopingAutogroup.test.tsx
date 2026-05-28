import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { ScopingAutogroupCard } from "../src/components/ScopingAutogroupCard";
import { ScopingOrderField } from "../src/components/ScopingOrderField";

describe("ScopingAutogroupCard", () => {
  it("summarises the grouping: member count, key label, flag count", () => {
    render(<ScopingAutogroupCard memberCount={6} keyLabel="ratio" flagCount={1} />);
    const card = screen.getByTestId("scoping-autogroup");
    expect(card).toHaveTextContent("6 samples");
    expect(card).toHaveTextContent("ratio");
    expect(card).toHaveTextContent(/one needs a look/i);
  });
  it("reads cleanly when nothing is flagged", () => {
    render(<ScopingAutogroupCard memberCount={5} keyLabel="ratio" flagCount={0} />);
    expect(screen.getByTestId("scoping-autogroup")).toHaveTextContent(/parsed cleanly/i);
  });
});

describe("ScopingOrderField", () => {
  it("shows the humanised ordering key as the field value", () => {
    render(<ScopingOrderField keyLabel="ll37 lipid ratio" />);
    expect(screen.getByTestId("scoping-order-field")).toHaveTextContent("ll37 lipid ratio");
  });
});
