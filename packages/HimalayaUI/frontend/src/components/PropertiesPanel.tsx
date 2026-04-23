import { useState } from "react";
import { ExposuresTab } from "./ExposuresTab";
import { PeaksTab } from "./PeaksTab";
import { TagsTab } from "./TagsTab";
import { NotesTab } from "./NotesTab";

type TabKey = "exposures" | "peaks" | "tags" | "notes";

interface TabDef {
  key: TabKey;
  label: string;
  Component: () => JSX.Element;
}

const TABS: readonly TabDef[] = [
  { key: "exposures", label: "Exposures", Component: ExposuresTab },
  { key: "peaks",     label: "Peaks",     Component: PeaksTab     },
  { key: "tags",      label: "Tags",      Component: TagsTab      },
  { key: "notes",     label: "Notes",     Component: NotesTab     },
] as const;

export function PropertiesPanel(): JSX.Element {
  const [active, setActive] = useState<TabKey>("exposures");
  const ActiveBody = TABS.find((t) => t.key === active)!.Component;

  return (
    <div className="flex flex-col h-full">
      <div role="tablist" className="flex gap-1 border-b border-border pb-1 mb-2">
        {TABS.map((t) => {
          const selected = t.key === active;
          return (
            <button
              key={t.key}
              id={`tab-${t.key}`}
              role="tab"
              aria-selected={selected}
              aria-controls={`tabpanel-${t.key}`}
              className={
                "px-3 py-1 rounded-md text-[13px] focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent " +
                (selected
                  ? "bg-bg-elevated text-fg border-b-2 border-accent"
                  : "text-fg-muted hover:text-fg hover:bg-bg-hover")
              }
              onClick={() => setActive(t.key)}
            >
              {t.label}
            </button>
          );
        })}
      </div>
      <div
        role="tabpanel"
        id={`tabpanel-${active}`}
        aria-labelledby={`tab-${active}`}
        className="flex-1 overflow-auto"
      >
        <ActiveBody />
      </div>
    </div>
  );
}
