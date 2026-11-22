export function suppressGlobalKeys(e: KeyboardEvent): boolean {
  if (e.defaultPrevented) return true;
  const t = e.target;
  if (t instanceof HTMLElement) {
    if (isEditable(t)) return true;
    if (t.closest('[role="dialog"], [data-keys-trap]') !== null) return true;
  }
  return document.querySelector('[role="dialog"][aria-modal="true"]') !== null;
}
export function isEditable(t: EventTarget | null): boolean {
  if (!(t instanceof HTMLElement)) return false;
  const tag = t.tagName;
  if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return true;
  if (t.isContentEditable) return true;
  const editableHost = t.closest("[contenteditable]");
  return editableHost !== null && editableHost.getAttribute("contenteditable") !== "false";
}
export function inDialog(e: KeyboardEvent): boolean {
  const t = e.target;
  if (t instanceof Element && t.closest('[role="dialog"], [data-keys-trap]') !== null) return true;
  return document.querySelector('[role="dialog"][aria-modal="true"]') !== null;
}
export function helpKeyFromEmptyField(e: KeyboardEvent): boolean {
  if (e.key !== "?" || e.defaultPrevented || e.metaKey || e.ctrlKey || e.altKey || inDialog(e)) return false;
  const t = e.target;
  return t instanceof HTMLInputElement && t.value === "";
}