import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { FigureExportControls } from "../../src/components/FigureExportControls";
import * as renderer from "../../src/lib/figure-export/renderer";
import * as clipboard from "../../src/lib/figure-export/clipboard";
import * as toastModule from "../../src/lib/toast";
import type { ExportSpec } from "../../src/lib/figure-export/types";

const fakeSpec: ExportSpec = {
  title: { primary: "P" }, width: 100, height: 80,
  plot: { marks: [], width: 80, height: 60 },
};

let toastSpy: ReturnType<typeof vi.fn>;

beforeEach(() => {
  toastSpy = vi.fn();
  toastModule.setToastImpl(toastSpy);
});
afterEach(() => {
  toastModule.setToastImpl(null);
  vi.restoreAllMocks();
});

describe("FigureExportControls", () => {
  it("renders Copy and Download buttons with aria labels using ariaContext", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="himalaya-trace-x"
        ariaContext="trace plot"
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot to clipboard/i }))
      .toBeInTheDocument();
    expect(screen.getByRole("button", { name: /download trace plot as png/i }))
      .toBeInTheDocument();
  });

  it("Copy click invokes renderer + clipboard and emits success toast", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
    const writeSpy = vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: /copy trace plot/i }));

    await waitFor(() => expect(writeSpy).toHaveBeenCalledWith(blob));
    expect(toastSpy).toHaveBeenCalledWith(
      expect.stringContaining("Copied"),
      "success",
    );
  });

  it("Copy failure surfaces error toast (renderer rejects)", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: /copy trace plot/i }));

    await waitFor(() =>
      expect(toastSpy).toHaveBeenCalledWith(
        expect.stringContaining("Couldn't copy"),
        "error",
      ),
    );
  });

  it("buttons re-enable after a thrown render (try/finally pending flag)", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const copyBtn = screen.getByRole("button", { name: /copy trace plot/i });
    fireEvent.click(copyBtn);

    await waitFor(() => expect(copyBtn).not.toBeDisabled());
  });

  it("Copy is disabled when canCopyPngToClipboard() is false", () => {
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(false);
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot/i })).toBeDisabled();
  });

  it("disabled prop disables both buttons (parent gate)", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
        disabled
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot/i })).toBeDisabled();
    expect(screen.getByRole("button", { name: /download trace plot as png/i })).toBeDisabled();
  });

  it("Chevron toggles a popover with PNG/SVG menu rows", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const chevron = screen.getByRole("button", { name: /other download formats/i });
    expect(chevron).toHaveAttribute("aria-haspopup", "true");
    expect(chevron).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(chevron);
    expect(chevron).toHaveAttribute("aria-expanded", "true");
    expect(screen.getByText(/download as png/i)).toBeInTheDocument();
    expect(screen.getByText(/download as svg/i)).toBeInTheDocument();
  });

  it("Esc closes the popover", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const chevron = screen.getByRole("button", { name: /other download formats/i });
    fireEvent.click(chevron);
    fireEvent.keyDown(document, { key: "Escape" });
    expect(chevron).toHaveAttribute("aria-expanded", "false");
  });

  it("spec is invoked on each click (thunk captures fresh state)", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const specThunk = vi.fn().mockReturnValue(fakeSpec);
    render(
      <FigureExportControls
        spec={specThunk}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const copyBtn = screen.getByRole("button", { name: /copy trace plot/i });
    fireEvent.click(copyBtn);
    await waitFor(() => expect(specThunk).toHaveBeenCalledTimes(1));
    fireEvent.click(copyBtn);
    await waitFor(() => expect(specThunk).toHaveBeenCalledTimes(2));
  });
});
