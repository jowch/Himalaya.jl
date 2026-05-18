/**
 * SamplesPage — placeholder body for the corpus contact sheet.
 *
 * Issue #155 (I1.1) only establishes the route slot so the corpus shell has
 * a body to mount. Issue #160 (I1.4) builds the real contact-sheet table
 * into this file.
 */
export function SamplesPage(): JSX.Element {
  return (
    <div data-testid="samples-page" className="p-10 text-sm text-ink-faint">
      Contact sheet — coming soon.
    </div>
  );
}
