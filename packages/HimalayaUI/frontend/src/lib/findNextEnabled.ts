/**
 * Roving-focus traversal: from index `from`, step by `delta` (wrapping) to the
 * next item whose `disabled` is falsy. Returns that index, or null if none are
 * enabled. Callers decide what to do with it (focus, select, both).
 */
export function findNextEnabled(
  options: ReadonlyArray<{ disabled?: boolean }>,
  from: number,
  delta: number,
): number | null {
  const n = options.length;
  if (n === 0) return null;
  let i = from;
  for (let step = 0; step < n; step++) {
    i = (i + delta + n) % n;
    if (!options[i].disabled) return i;
  }
  return null;
}
