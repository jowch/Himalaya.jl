export function callChoices(root: ParentNode = document): HTMLElement[] {
  return Array.from(root.querySelectorAll<HTMLElement>("[data-call-choices] button:not(:disabled)"));
}

export function inCallChoices(target: EventTarget | null): boolean {
  return target instanceof Element && target.closest("[data-call-choices]") !== null;
}

export function cycleChoice(choices: readonly HTMLElement[], current: EventTarget | null, back: boolean): HTMLElement | undefined {
  const n = choices.length;
  if (n === 0) return undefined;
  const i = current instanceof Node ? choices.findIndex((c) => c.contains(current)) : -1;
  if (i < 0) return back ? choices[n - 1] : choices[0];
  return choices[(i + (back ? n - 1 : 1)) % n];
}
