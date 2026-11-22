import type { HelpSection } from "../../workspace/ui/HelpOverlay";
import { MOD_KEY } from "../../lib/modKey";

const GENERAL: HelpSection = {
  title: "General",
  rows: [
    { keys: [MOD_KEY, "K"], what: "Commands and go to an experiment or sample" },
    { keys: ["/"], what: "The same, when no field has the cursor" },
    { keys: ["?"], what: "This list" },
  ],
};

const HOME: HelpSection = {
  title: "Experiments",
  rows: [
    { keys: ["↑", "↓", "←", "→"], what: "Move between experiments" },
    { keys: ["Enter"], what: "Open the experiment" },
  ],
};

export function shellHelp(pathname: string): HelpSection[] {
  if (pathname === "/experiments") return [HOME, GENERAL];
  return [GENERAL];
}
