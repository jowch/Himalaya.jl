/**
 * suppressGlobalKeys — the shared guard for window-level page shortcuts
 * (loupe X/R/arrows, contact-sheet X/Escape, global `/`-`,`-`.`, scoping ⌘Z).
 *
 * Returns `true` (the page shortcut must NOT run) when either:
 *
 * 1. The event target is a typing/selection context — INPUT, TEXTAREA,
 *    SELECT, or inside a contenteditable region. Single letters and chords
 *    like ⌘Z belong to the editing surface there, never to the page.
 *
 * 2. A modal dialog is open ANYWHERE in the document:
 *    `[role="dialog"][aria-modal="true"]` (what ModalShell emits).
 *    This is deliberately a document-level check, not `target.closest(…)`:
 *    `aria-modal="true"` declares everything outside the dialog inert, so
 *    page shortcuts must be off no matter where the event lands — and focus
 *    can transiently sit on <body> (target = body/window) in the instant
 *    before the modal's focus trap engages.
 *
 * Modifier-key chords (⌘/Ctrl/Alt) are deliberately NOT handled here: some
 * callers bind chords themselves (useGlobalShortcuts owns ⌘K, the scoping
 * page owns ⌘Z), so each mutating call site keeps its own
 * `if (e.metaKey || e.ctrlKey || e.altKey) return;` line where appropriate.
 *
 * The contenteditable check pairs the spec'd `isContentEditable` property
 * (the browser source of truth, covers inheritance + designMode) with an
 * attribute-based `closest` fallback: JSDOM does not implement
 * `isContentEditable`, and the fallback also resolves inherited editability
 * while respecting the nearest explicit `contenteditable="false"` island.
 */
export function suppressGlobalKeys(e: KeyboardEvent): boolean {
  const t = e.target;
  if (t instanceof HTMLElement) {
    const tag = t.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
    if (t.isContentEditable) return true;
    const editableHost = t.closest("[contenteditable]");
    if (editableHost !== null && editableHost.getAttribute("contenteditable") !== "false") {
      return true;
    }
  }
  return document.querySelector('[role="dialog"][aria-modal="true"]') !== null;
}
