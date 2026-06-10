import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { FitMetadata } from "../../src/print/components/FitMetadata";

describe("<FitMetadata>", () => {
  it("renders with data-testid=fit-metadata", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toBeInTheDocument();
  });

  it("names reflections as the subject with a plural verb when landed=0", () => {
    render(
      <FitMetadata landed={0} total={6} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent(
      "0 of 6 reflections land on observed peaks"
    );
  });

  it("uses the singular verb 'lands' when landed=1", () => {
    render(
      <FitMetadata landed={1} total={6} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent(
      "1 of 6 reflections lands on observed peaks"
    );
  });

  it("uses the plural verb 'land' when landed=4", () => {
    render(
      <FitMetadata landed={4} total={6} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent(
      "4 of 6 reflections land on observed peaks"
    );
  });

  it("contains the landed count", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("7");
  });

  it("contains the total count", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("9");
  });

  it("contains the paramName", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("d");
  });

  it("contains the paramValue", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("81.2");
  });

  it("contains the unit", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("Å");
  });

  it("renders fit-snapped element when snapped=true", () => {
    render(
      <FitMetadata
        landed={7}
        total={9}
        paramName="d"
        paramValue="81.2"
        unit="Å"
        snapped
      />
    );
    expect(screen.getByTestId("fit-snapped")).toBeInTheDocument();
  });

  it("does NOT render fit-snapped element when snapped=false", () => {
    render(
      <FitMetadata
        landed={7}
        total={9}
        paramName="d"
        paramValue="81.2"
        unit="Å"
        snapped={false}
      />
    );
    expect(screen.queryByTestId("fit-snapped")).not.toBeInTheDocument();
  });

  it("does NOT render fit-snapped element when snapped is omitted", () => {
    render(
      <FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" unit="Å" />
    );
    expect(screen.queryByTestId("fit-snapped")).not.toBeInTheDocument();
  });

  it("uses default unit Å when unit is omitted", () => {
    render(<FitMetadata landed={7} total={9} paramName="d" paramValue="81.2" />);
    expect(screen.getByTestId("fit-metadata")).toHaveTextContent("Å");
  });
});
