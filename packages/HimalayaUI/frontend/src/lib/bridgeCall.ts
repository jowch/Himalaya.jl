import type { AssignmentState, CorpusSample } from "../api";

export type KeepCall = number | "both";

type Listed = Pick<CorpusSample, "id" | "name" | "assignment_state" | "call_phases" | "phase">;

function state(s: Listed): AssignmentState {
  return s.assignment_state ?? "indexed";
}

function phases(s: Listed): string[] {
  return [...(s.call_phases ?? (s.phase ? [s.phase] : []))].sort();
}

export function isCalled(s: Listed): boolean {
  return state(s) === "form_factor" || (state(s) === "indexed" && phases(s).length > 0);
}

export function callText(s: Listed): string {
  if (state(s) === "form_factor") return "form factor";
  return phases(s).join(" + ") || "no call";
}

export function callsConflict(a: Listed, b: Listed): boolean {
  return isCalled(a) && isCalled(b) && (state(a) !== state(b) || phases(a).join("\n") !== phases(b).join("\n"));
}

export interface KeepChoice {
  value: string;
  label: string;
  keep: KeepCall;
}

export function keepChoices(here: Listed, there: Listed, name: (s: Listed) => string = (s) => s.name): KeepChoice[] | undefined {
  if (!callsConflict(here, there)) return undefined;
  const out: KeepChoice[] = [
    { value: String(here.id), label: `${name(here)}: ${callText(here)}`, keep: here.id },
    { value: String(there.id), label: `${name(there)}: ${callText(there)}`, keep: there.id },
  ];
  if (state(here) === "indexed" && state(there) === "indexed") {
    out.push({ value: "both", label: `Both: ${[...new Set([...phases(here), ...phases(there)])].sort().join(" + ")}`, keep: "both" });
  }
  return out;
}

export function callStatesSplit(a: Listed, b: Listed): boolean {
  return isCalled(a) && isCalled(b) && state(a) !== state(b);
}
