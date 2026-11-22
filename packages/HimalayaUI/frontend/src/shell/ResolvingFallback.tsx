export function ResolvingFallback(): JSX.Element {
  return (
    <div
      data-testid="resolving"
      className="flex-1 min-h-0 flex flex-col"
      aria-busy="true"
      aria-label="Loading page"
    />
  );
}
