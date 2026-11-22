import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { ExperimentNotice } from "../../src/workspace/ui/ExperimentNotice";
import { ProblemList } from "../../src/workspace/ui/ProblemList";

describe("setup folders a rescan skips", () => {
  it("are counted in the notice, which offers the details", () => {
    const onDetails = vi.fn();
    render(<ExperimentNotice progress={undefined} ingestFailed={false} waiting={0} failed={0} unrecognized={2} onDetails={onDetails} />);
    expect(screen.getByTestId("workspace-ingest-status")).toHaveTextContent("2 setup folders not recognised");
    fireEvent.click(screen.getByRole("button", { name: "Details" }));
    expect(onDetails).toHaveBeenCalled();
  });

  it("are named in the problem list, with a way to add them by hand", () => {
    const onAddSetups = vi.fn();
    render(
      <ProblemList
        problems={{ ingest_error: null, exposures: [], setups: [{ root: "/b/setup3", folder: "setup3" }] }}
        onPick={vi.fn()} onAddSetups={onAddSetups} onClose={vi.fn()}
      />,
    );
    expect(screen.getByTestId("problems-setups")).toHaveTextContent("A setup folder was not added: setup3.");
    expect(screen.queryByText("No problems.")).toBeNull();
    fireEvent.click(screen.getByTestId("problems-add-setups"));
    expect(onAddSetups).toHaveBeenCalled();
  });
});
