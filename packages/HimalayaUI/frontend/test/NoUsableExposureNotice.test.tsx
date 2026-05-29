import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { NoUsableExposureNotice } from "../src/components/NoUsableExposureNotice";

function renderNotice(props: { sampleId: number; allRejected: boolean }) {
  return render(
    <MemoryRouter>
      <NoUsableExposureNotice {...props} />
    </MemoryRouter>,
  );
}

describe("NoUsableExposureNotice", () => {
  it("links back to the loupe for the sample (the reversible-rejection surface)", () => {
    renderNotice({ sampleId: 42, allRejected: true });
    const link = screen.getByTestId("no-usable-exposure-loupe-link");
    expect(link).toHaveAttribute("href", "/samples/loupe/42");
  });

  it("explains the all-rejected case and offers to restore", () => {
    renderNotice({ sampleId: 7, allRejected: true });
    expect(screen.getByTestId("no-usable-exposure")).toBeInTheDocument();
    // The message must distinguish this from a generic load failure.
    expect(screen.getByText(/rejected/i)).toBeInTheDocument();
    expect(screen.getByTestId("no-usable-exposure-loupe-link")).toHaveTextContent(/restore/i);
  });

  it("explains the no-exposures case without the restore framing", () => {
    renderNotice({ sampleId: 7, allRejected: false });
    expect(screen.getByText("No exposures yet")).toBeInTheDocument();
    // "Restore" only makes sense when something was rejected.
    expect(screen.getByTestId("no-usable-exposure-loupe-link")).not.toHaveTextContent(/restore/i);
  });
});
