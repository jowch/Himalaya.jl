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

function sameCall(a: Listed, b: Listed): boolean {
  return isCalled(a) && isCalled(b) && !callsConflict(a, b);
}

function unionCall(a: Listed, b: Listed): Listed {
  return { ...a, assignment_state: "indexed", call_phases: [...new Set([...phases(a), ...phases(b)])].sort() };
}

export function groupKeepChoices(here: Listed, mates: readonly Listed[], name: (s: Listed) => string = (s) => s.name): KeepChoice[] | undefined {
  const calls: Listed[][] = [];
  for (const s of [here, ...mates]) {
    if (!isCalled(s)) continue;
    const same = calls.find((c) => sameCall(c[0]!, s));
    if (same) same.push(s);
    else calls.push([s]);
  }
  if (calls.length < 2) return undefined;
  const out: KeepChoice[] = calls.map((c) => ({ value: String(c[0]!.id), label: `${c.map(name).join(", ")}: ${callText(c[0]!)}`, keep: c[0]!.id }));
  if (calls.every((c) => state(c[0]!) === "indexed")) {
    const all = [...new Set(calls.flatMap((c) => phases(c[0]!)))].sort().join(" + ");
    out.push({ value: "both", label: `${calls.length === 2 ? "Both" : "All"}: ${all}`, keep: "both" });
  }
  return out;
}

export interface BridgeStep {
  to: number;
  keep?: KeepCall;
}

export function bridgeSteps(here: Listed, mates: readonly Listed[], keep: KeepCall | undefined): BridgeStep[] {
  const kept = typeof keep === "number" ? [here, ...mates].find((s) => s.id === keep) : undefined;
  let group: Listed = here;
  return mates.map((m) => {
    if (!callsConflict(group, m)) {
      if (!isCalled(group) && isCalled(m)) group = m;
      return { to: m.id };
    }
    if (keep === "both") {
      group = unionCall(group, m);
      return { to: m.id, keep: "both" };
    }
    if (kept && sameCall(kept, m)) {
      group = m;
      return { to: m.id, keep: m.id };
    }
    return { to: m.id, keep: here.id };
  });
}
