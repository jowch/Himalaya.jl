/**
 * ComparePage — placeholder. A future multi-sample comparison workflow.
 *
 * Kept minimal so the tab-rocker can switch to it without breaking; real
 * content will be designed in a follow-up iteration.
 */
export function ComparePage(): JSX.Element {
  return (
    <div
      data-testid="compare-page"
      className="flex-1 min-h-0 flex items-center justify-center px-4 pb-4"
    >
      <div className="card max-w-lg w-full p-8 flex flex-col items-center gap-2 text-center">
        <div className="text-xs uppercase tracking-widest text-fg-dim">
          Compare
        </div>
        <h2 className="text-title">
          Coming soon
        </h2>
        <p className="text-base text-fg-muted max-w-sm">
          Multi-sample comparison view. For now, use <kbd className="border border-border rounded px-1 text-xs">/</kbd>
          {" "}to switch between samples one at a time on the Index page.
        </p>
      </div>
    </div>
  );
}
