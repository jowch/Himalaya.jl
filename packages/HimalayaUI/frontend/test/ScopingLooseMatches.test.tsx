import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScopingLooseMatches } from "../src/components/ScopingLooseMatches";

const loose = { sampleId: 21, sampleName: "Lipid 1-1 + LL37 1:1", value: "", flagged: true, include: false };

describe("ScopingLooseMatches", () => {
  it("lists each loose match with an Add affordance", () => {
    render(<ScopingLooseMatches rows={[loose]} traces={new Map()} phases={new Map()} onAdd={() => {}} />);
    expect(screen.getByTestId("scoping-loose-21")).toHaveTextContent("Lipid 1-1 + LL37 1:1");
    expect(screen.getByTestId("scoping-loose-add-21")).toBeInTheDocument();
  });
  it("clicking Add threads the sample id up", () => {
    const onAdd = vi.fn();
    render(<ScopingLooseMatches rows={[loose]} traces={new Map()} phases={new Map()} onAdd={onAdd} />);
    fireEvent.click(screen.getByTestId("scoping-loose-add-21"));
    expect(onAdd).toHaveBeenCalledWith(21);
  });
  it("shows the empty note when nothing else matches", () => {
    render(<ScopingLooseMatches rows={[]} traces={new Map()} phases={new Map()} onAdd={() => {}} />);
    expect(screen.getByTestId("scoping-loose-empty")).toBeInTheDocument();
  });
  it("collapses a long list behind 'show N more' and expands on click", () => {
    const many = Array.from({ length: 7 }, (_, i) => ({
      sampleId: 100 + i, sampleName: `S${i}`, value: "", flagged: true, include: false,
    }));
    render(<ScopingLooseMatches rows={many} traces={new Map()} phases={new Map()} onAdd={() => {}} />);
    // 4 shown, the 5th hidden until expanded.
    expect(screen.getByTestId("scoping-loose-103")).toBeInTheDocument();
    expect(screen.queryByTestId("scoping-loose-104")).toBeNull();
    fireEvent.click(screen.getByTestId("scoping-loose-more"));
    expect(screen.getByTestId("scoping-loose-106")).toBeInTheDocument();
  });
});
