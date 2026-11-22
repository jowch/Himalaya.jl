import { showToast } from "../toast";
export function buildAssignmentDroppedMessage(dropped: string[]): string {
  const counts = new Map<string, number>();
  for (const phase of dropped) counts.set(phase, (counts.get(phase) ?? 0) + 1);
  const groups = [...counts.entries()].map(
    ([phase, n]) => (n > 1 ? `${n} ${phase}` : phase),
  );
  const head = groups.length === 1
    ? groups[0]!
    : `${groups.slice(0, -1).join(", ")} and ${groups[groups.length - 1]!}`;
  const total = dropped.length;
  const noun = total === 1 ? "index" : "indices";
  const tail = total === 1
    ? "Its index no longer fits the peaks."
    : "They no longer fit the peaks.";
  return `${head} ${noun} dropped from the call. ${tail}`;
}
export function announceAssignmentDropped(dropped: string[] | undefined): void {
  if (!dropped || dropped.length === 0) return;
  showToast(buildAssignmentDroppedMessage(dropped), "warning");
}