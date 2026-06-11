import { Kicker, MetaList, TagList, KbLegend } from "../ui";
import type { MetaEntry, Tag, Shortcut } from "../ui";
import { Verdict } from "./Verdict";
import { RepresentativeBox } from "./RepresentativeBox";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface LoupeSidePanelProps {
  /** "This exposure" metadata rows. */
  meta: MetaEntry[];
  /** Verdict state + toggles (tri-state: dropped / kept / unscreened). */
  dropped: boolean;
  kept?: boolean;
  onToggleDrop?: () => void;
  onToggleKeep?: () => void;
  /** Representative state + setter. */
  isRepresentative: boolean;
  /** The sample's representative exposure (any frame) is dropped. */
  representativeDropped?: boolean;
  onSetRepresentative?: () => void;
  /** Sample tags + edit handlers. */
  tags: Tag[];
  onAddTag?: (t: Tag) => void;
  onRemoveTag?: (t: Tag) => void;
  /** Key legend — defaults to the canonical loupe keys. */
  shortcuts?: Shortcut[];
  /** PLACEMENT ONLY. */
  className?: string;
}

export const LOUPE_KEYS: Shortcut[] = [
  { keyLabel: "← →", description: "flip frames" },
  { keyLabel: "X", description: "drop / restore" },
  { keyLabel: "K", description: "keep / restore" },
  { keyLabel: "R", description: "set representative" },
  { keyLabel: "Esc", description: "back to the sheet" },
];

export function LoupeSidePanel({
  meta,
  dropped,
  kept,
  onToggleDrop,
  onToggleKeep,
  isRepresentative,
  representativeDropped,
  onSetRepresentative,
  tags,
  onAddTag,
  onRemoveTag,
  shortcuts = LOUPE_KEYS,
  className,
}: LoupeSidePanelProps): JSX.Element {
  return (
    <aside
      data-testid="loupe-side-panel"
      className={cx("flex flex-col gap-[18px]", className)}
    >
      {/* Block 1: This exposure */}
      <div>
        <Kicker tone="faint" as="h3" className="mb-2">
          This exposure
        </Kicker>
        <MetaList entries={meta} />
      </div>

      {/* Block 2: Verdict */}
      <Verdict
        dropped={dropped}
        {...(kept !== undefined ? { kept } : {})}
        {...(onToggleDrop ? { onToggle: onToggleDrop } : {})}
        {...(onToggleKeep ? { onToggleKeep } : {})}
      />

      {/* Block 3: RepresentativeBox */}
      <RepresentativeBox
        isRepresentative={isRepresentative}
        {...(representativeDropped !== undefined ? { representativeDropped } : {})}
        {...(onSetRepresentative ? { onSetRepresentative } : {})}
      />

      {/* Block 4: Sample tags */}
      <div>
        <Kicker tone="faint" as="h3" className="mb-2">
          Sample tags
        </Kicker>
        <TagList
          tags={tags}
          editable
          persistentAdd
          {...(onAddTag ? { onAdd: onAddTag } : {})}
          {...(onRemoveTag ? { onRemove: onRemoveTag } : {})}
        />
      </div>

      {/* Block 5: Keys */}
      <div>
        <Kicker tone="faint" as="h3" className="mb-2">
          Keys
        </Kicker>
        <KbLegend shortcuts={shortcuts} className="flex-col gap-1" />
      </div>
    </aside>
  );
}
