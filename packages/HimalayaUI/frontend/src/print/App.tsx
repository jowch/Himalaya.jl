/**
 * PrintApp — root of the greenfield ("The Print") UI tree.
 *
 * Foundation placeholder. Routing, SSE wiring, and the mutation-queue
 * persistence/rehydrate effects (mirrored from the old src/App.tsx) are added
 * in the pages-assembly plan. For now this proves the src/print/ entry builds,
 * renders, and is reachable via print.html.
 */
export function PrintApp(): JSX.Element {
  return (
    <main data-testid="print-shell" className="min-h-screen bg-paper text-ink">
      <h1 className="text-display">The Print · greenfield shell</h1>
    </main>
  );
}
