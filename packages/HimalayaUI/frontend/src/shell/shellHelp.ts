import type { HelpSection } from "../ui/HelpOverlay";
import { MOD_KEY } from "../lib/modKey";

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

const NEW_EXPERIMENT: HelpSection = {
  title: "Add a beamtime",
  rows: [
    { keys: ["↑", "↓"], what: "Move between folders (Tab and Shift+Tab too)" },
    { keys: ["Enter"], what: "Add the beamtime or setup, or open a plain folder" },
    { keys: ["→"], what: "Open the folder" },
    { keys: ["←"], what: "Go up one folder (Esc too)" },
    { keys: ["/"], what: "Type or paste a path" },
    { keys: [MOD_KEY, "Enter"], what: "Add the experiment" },
  ],
};

export function shellHelp(pathname: string): HelpSection[] {
  if (pathname === "/experiments") return [HOME, GENERAL];
  if (pathname === "/experiments/new") return [NEW_EXPERIMENT, GENERAL];
  return [GENERAL];
}
