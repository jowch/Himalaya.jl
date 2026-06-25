/** Normalized combo: [Mod][+Alt][+Shift]+<key>. `Mod` folds meta/ctrl so the
 *  same declaration works on macOS and the rest. Letters lowercase; `?` and the
 *  Arrow/Enter/Escape names pass through; " " becomes "Space". */
export function comboOf(e: KeyboardEvent): string {
  const parts: string[] = [];
  if (e.metaKey || e.ctrlKey) parts.push("Mod");
  if (e.altKey) parts.push("Alt");
  if (e.shiftKey) parts.push("Shift");
  let k = e.key;
  if (k === " ") k = "Space";
  else if (k.length === 1) k = k.toLowerCase();
  parts.push(k);
  return parts.join("+");
}

export function matchesKeys(e: KeyboardEvent, keys: string[]): boolean {
  const got = comboOf(e);
  // Normalize each declared key through the same lowering so "K" === "k".
  return keys.some((want) => {
    const parts = want.split("+");
    const k = parts.pop()!;
    const norm = [...parts, k.length === 1 ? k.toLowerCase() : k].join("+");
    return norm === got;
  });
}

export function isBareKey(e: KeyboardEvent): boolean {
  return !e.metaKey && !e.ctrlKey && !e.altKey && e.key.length === 1;
}

export function isTyping(t: EventTarget | null): boolean {
  if (!(t instanceof HTMLElement)) return false;
  const tag = t.tagName;
  if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
  return t.isContentEditable || t.getAttribute("contenteditable") === "true";
}
