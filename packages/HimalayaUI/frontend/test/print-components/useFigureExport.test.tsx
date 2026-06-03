import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { useFigureExport } from "../../src/print/components/useFigureExport";
import * as renderer from "../../src/lib/figure-export/renderer";
import * as clipboard from "../../src/lib/figure-export/clipboard";
import * as download from "../../src/lib/figure-export/download";
import * as filename from "../../src/lib/figure-export/filename";
import * as toastModule from "../../src/lib/toast";
import type { ExportSpec } from "../../src/lib/figure-export/types";

const fakeSpec: ExportSpec = {
  title: { primary: "P" }, width: 100, height: 80,
  plot: { marks: [], width: 80, height: 60 },
};
const specThunk = (): ExportSpec => fakeSpec;

let toastSpy: ReturnType<typeof vi.fn>;

beforeEach(() => {
  toastSpy = vi.fn();
  toastModule.setToastImpl(toastSpy);
  vi.spyOn(renderer, "canExportPng").mockReturnValue(true);
  vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
});
afterEach(() => {
  toastModule.setToastImpl(null);
  vi.restoreAllMocks();
});

describe("useFigureExport", () => {
  it("onCopy renders PNG, writes clipboard, emits success toast", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    const write = vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onCopy(); });

    await waitFor(() => expect(write).toHaveBeenCalledWith(blob));
    expect(toastSpy).toHaveBeenCalledWith(
      expect.stringContaining("Copied the trace to clipboard"), "success",
    );
  });

  it("onCopy failure emits an error toast", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onCopy(); });

    await waitFor(() =>
      expect(toastSpy).toHaveBeenCalledWith(
        expect.stringContaining("Couldn't copy figure"), "error",
      ),
    );
  });

  it("onDownloadPng renders PNG and downloads with the built filename", async () => {
    const blob = new Blob([new Uint8Array([2])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(filename, "buildFilename").mockReturnValue("stem-2099-01-01.png");
    const dl = vi.spyOn(download, "downloadBlob").mockImplementation(() => {});

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onDownloadPng(); });

    await waitFor(() => expect(dl).toHaveBeenCalledWith(blob, "stem-2099-01-01.png"));
  });

  it("onDownloadSvg serializes the SVG and downloads it", async () => {
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    vi.spyOn(renderer, "buildExportSvg").mockReturnValue(svg);
    vi.spyOn(filename, "buildFilename").mockReturnValue("stem-2099-01-01.svg");
    const dl = vi.spyOn(download, "downloadBlob").mockImplementation(() => {});

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "the trace"));
    act(() => { result.current.onDownloadSvg(); });

    await waitFor(() => {
      expect(dl).toHaveBeenCalledTimes(1);
      const [blobArg, nameArg] = dl.mock.calls[0];
      expect(blobArg).toBeInstanceOf(Blob);
      expect((blobArg as Blob).type).toBe("image/svg+xml");
      expect(nameArg).toBe("stem-2099-01-01.svg");
    });
  });

  it("copyDisabled is true when clipboard is unsupported", () => {
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(false);
    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    expect(result.current.copyDisabled).toBe(true);
  });

  it("copyDisabled and pngDisabled are true when PNG export is unsupported", () => {
    vi.spyOn(renderer, "canExportPng").mockReturnValue(false);
    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    expect(result.current.copyDisabled).toBe(true);
    expect(result.current.pngDisabled).toBe(true);
  });

  it("guards against a double-invocation while a render is in flight", async () => {
    const build = vi.spyOn(renderer, "buildExportPng").mockReturnValue(new Promise<Blob>(() => {}));
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const { result } = renderHook(() => useFigureExport(specThunk, "stem", "x"));
    act(() => { result.current.onCopy(); result.current.onCopy(); });

    expect(build).toHaveBeenCalledTimes(1);
  });
});
