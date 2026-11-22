export function replacePlaceholder<T extends { id: number }>(
  list: T[],
  serverItem: T,
  matches: (item: T) => boolean,
  options?: { isDuplicate?: (item: T) => boolean },
): T[] {
  const isDup = options?.isDuplicate ?? ((item: T) => item.id === serverItem.id);
  let replaced = false;
  const out: T[] = [];
  for (const item of list) {
    if (isDup(item)) continue;
    if (!replaced && item.id < 0 && matches(item)) {
      out.push(serverItem);
      replaced = true;
      continue;
    }
    out.push(item);
  }
  if (!replaced) out.push(serverItem);
  return out;
}