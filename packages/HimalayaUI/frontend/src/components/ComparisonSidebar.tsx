/**
 * ComparisonSidebar — placeholder shell (Plan §Phase 4, Task 4.2).
 *
 * Real implementation lands in Task 4.2; this placeholder exists so the
 * ComparePage / ComparePageEdit page shells (Task 4.1) can mount without
 * pulling in the listing+search+toggle UI prematurely.
 */
interface Props {
  experimentId: number | undefined;
  scope: "experiment" | "all";
  activeComparisonId: number | undefined;
}

export function ComparisonSidebar(_props: Props): JSX.Element {
  return (
    <div data-testid="comparison-sidebar" className="flex-1 min-h-0 p-3" />
  );
}
