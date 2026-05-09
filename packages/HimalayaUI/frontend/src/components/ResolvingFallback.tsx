// Spec §4.2 — minimal sub-paint flicker mask. Renders inside <main>
// while a resolve is in flight. No <Skeleton> from boneyard — that
// would itself flash on every navigation.

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
