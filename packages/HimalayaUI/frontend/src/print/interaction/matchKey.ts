/** Normalized combo: [Mod][+Alt][+Shift]+<key>. `Mod` folds meta/ctrl so the
 *  same declaration works on macOS and the rest. Letters lowercase; `?` and the
 *  Arrow/Enter/Escape names pass through; " " becomes "Space". */
export function comboOf(e: KeyboardEvent): string {
  const parts: string[] = [];
  if (e.metaKey || e.ctrlKey) parts.push("Mod");
  if (e.altKey) parts.push("Alt");
  // Shift is a real modifier for named keys (Shift+ArrowUp) and alphanumerics
  // (Mod+Shift+z). But a shifted PUNCTUATION key already bakes Shift into the
  // produced character — a real `?` is Shift+/ and `e.key` is `"?"`, declared as
  // `"?"`. Adding "Shift" there yields "Shift+?" and the bare-"?" help shortcut
  // never matches. So skip Shift only for a single non-alphanumeric character.
  const shiftEncodedInChar = e.key.length === 1 && !/[A-Za-z0-9]/.test(e.key);
  if (e.shiftKey && !shiftEncodedInChar) parts.push("Shift");
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
    const kn = k === " " ? "Space" : k.length === 1 ? k.toLowerCase() : k;
    const norm = [...parts, kn].join("+");
    return norm === got;
  });
}

export function isTyping(t: EventTarget | null): boolean {
  if (!(t instanceof HTMLElement)) return false;
  const tag = t.tagName;
  if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
  if (t.isContentEditable) return true;
  const ce = t.getAttribute("contenteditable");
  return ce === "true" || ce === "";
}
