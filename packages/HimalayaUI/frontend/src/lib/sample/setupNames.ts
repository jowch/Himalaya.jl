import type { Setup } from "../../api";

export function setupNames(setups: readonly Pick<Setup, "id" | "name">[]): ReadonlyMap<number, string> | undefined {
  return setups.length > 1 ? new Map(setups.map((st) => [st.id, st.name])) : undefined;
}

export function joinHints(...parts: readonly (string | undefined)[]): string | undefined {
  const kept = parts.filter((p): p is string => Boolean(p));
  return kept.length ? kept.join(" · ") : undefined;
}
