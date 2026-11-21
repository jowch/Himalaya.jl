import type { HelpSection } from "../../workspace/ui/HelpOverlay";
import { MOD_KEY } from "../../lib/modKey";

const GENERAL: HelpSection = {
  title: "Anywhere",
  rows: [
    { keys: [MOD_KEY, "K"], what: "Find an experiment or sample" },
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
  title: "New experiment",
  rows: [
    { keys: ["Tab"], what: "Complete the folder path" },
    { keys: ["↑", "↓"], what: "Move through the folder suggestions" },
    { keys: [MOD_KEY, "Enter"], what: "Add the experiment" },
  ],
};

export function shellHelp(pathname: string): HelpSection[] {
  if (pathname === "/experiments/new") return [NEW_EXPERIMENT, GENERAL];
  if (pathname === "/experiments") return [HOME, GENERAL];
  return [GENERAL];
}
