// download.ts — trigger a browser download for an in-memory blob.
export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.rel = "noopener";
  a.style.display = "none";
  document.body.appendChild(a);
  a.click();
  // Defer cleanup. `a.click()` starts the download ASYNCHRONOUSLY, so revoking
  // the object URL (or tearing down the anchor) synchronously can race the
  // browser before it has read the blob AND the `download` filename. When it
  // loses that race the browser falls back to saving the blob URL's UUID as the
  // name — and may not place it in the downloads folder at all. Letting the
  // current task finish (and giving the download a moment to start) keeps the
  // anchor's `download` attribute authoritative and the figure correctly named.
  setTimeout(() => {
    a.remove();
    URL.revokeObjectURL(url);
  }, 10_000);
}
