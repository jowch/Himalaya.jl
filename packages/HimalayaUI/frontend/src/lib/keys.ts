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
 * 2. The event target is inside an open NON-modal popover — a
 *    `[role="dialog"]` panel (what the Popover primitive emits) or a
 *    `[data-keys-trap]` inline editor (the tag editor). An open popover owns
 *    the keyboard while focus is inside it, so a single-letter page shortcut
 *    must not fire through a non-input control (a button/chip) sitting in the
 *    popover (LO-POPKEY). This is `target.closest(…)`, NOT document-level: a
 *    non-modal popover does not make the rest of the page inert, so shortcuts
 *    stay live everywhere OUTSIDE the popover's own subtree.
 *
 * 3. A modal dialog is open ANYWHERE in the document:
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
  // A handler closer to the source already consumed this key (e.g. ModalShell's
  // Escape preventDefault()s before closing). The open-dialog check below is NOT
  // enough on its own: in a real browser a microtask checkpoint runs between the
  // dialog's `document` listener and a page's `window` listener, so React can
  // flush the close and UNMOUNT the dialog mid-dispatch — by the time we query
  // for `[aria-modal]` it is already gone and the page shortcut (loupe Esc →
  // back to the sheet) wrongly fires. `defaultPrevented` survives that race.
  if (e.defaultPrevented) return true;
  const t = e.target;
  if (t instanceof HTMLElement) {
    const tag = t.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
    if (t.isContentEditable) return true;
    const editableHost = t.closest("[contenteditable]");
    if (editableHost !== null && editableHost.getAttribute("contenteditable") !== "false") {
      return true;
    }
    // An open non-modal popover (Popover primitive `[role="dialog"]`) or inline
    // editor (`[data-keys-trap]`, e.g. the tag editor) owns the keyboard while
    // focus is inside it — page shortcuts must not fire through a button/chip
    // sitting in it (LO-POPKEY). Target-scoped: the page stays live outside.
    if (t.closest('[role="dialog"], [data-keys-trap]') !== null) return true;
  }
  return document.querySelector('[role="dialog"][aria-modal="true"]') !== null;
}

/**
 * isNativeInteractiveTarget — guard for `Enter` actions registered in the
 * keyboard layer (e.g. `openFocus` / "Apply"). When Enter lands on a native
 * interactive control, the control owns the event: a button activates, a link
 * follows, an input submits/edits, a select opens, a contenteditable inserts a
 * break. The page-level Enter action must yield rather than hijacking.
 *
 * Covers:
 *   • Typing contexts  — INPUT, TEXTAREA, SELECT, contenteditable
 *     (Enter inserts or submits; page navigation must not fire)
 *   • Activation targets — button, a, [role=button], sortable column header
 *     (Enter natively clicks the control; page navigation must not fire)
 *
 * Previously each page's `openFocus` binding called this and early-returned.
 * Now it is also enforced centrally in `useKeyboardLayer` so the guarantee is
 * branch-wide and page authors cannot forget it.
 */
export function isNativeInteractiveTarget(e: KeyboardEvent): boolean {
  if (!(e.target instanceof Element)) return false;
  // Typing contexts: Enter submits/edits rather than navigating.
  if (e.target instanceof HTMLElement) {
    const tag = e.target.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
    if (e.target.isContentEditable) return true;
    const editableHost = e.target.closest("[contenteditable]");
    if (editableHost !== null && editableHost.getAttribute("contenteditable") !== "false") return true;
  }
  // Activation targets: Enter natively activates the control.
  return e.target.closest("button, a, [role=button], [aria-sort]") != null;
}
