import type { RefObject } from "react";
import { Kicker, MetaList, TagList, KbLegend, Button } from "../ui";
import type { MetaEntry, Tag, Shortcut } from "../ui";
import type { LoupeTag } from "../pages/loupeAdapters";
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
  /** Sample tags + edit handlers. Accepts `LoupeTag[]` (carries `id`) for the
   *  id-exact remove path (LO-TAGDUP), or plain `Tag[]` for read-only contexts. */
  tags: LoupeTag[] | Tag[];
  onAddTag?: (t: Tag) => void;
  /** Id-exact remove: called with the tag's `id` (LO-TAGDUP fix). */
  onRemoveTag?: (id: number) => void;
  /** Called when the user clicks the "Manage" affordance next to "Sample tags".
   *  When absent the button is not rendered. */
  onManageTags?: () => void;
  /** Ref for the "Manage" button — forwarded to ManageTagsModal for focus-restore. */
  manageTagsTriggerRef?: RefObject<HTMLButtonElement>;
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
  onManageTags,
  manageTagsTriggerRef,
  shortcuts = LOUPE_KEYS,
  className,
}: LoupeSidePanelProps): JSX.Element {
  return (
    <aside
      data-testid="loupe-side-panel"
      className={cx("flex flex-col gap-[18px]", className)}
    >
      {/* Block 1: This frame (LO-TERM: the loupe says "frame", not "exposure").
          LO-EXPSPARSE: the only always-present row is the position ("N of M"),
          which the header subtitle AND the BigFrame caption already show, so a
          lone-row block reads as an unfinished section. Render it only once a
          SECOND, real row exists (e.g. a dropped frame's rejection reason);
          otherwise suppress the kicker entirely. (Reinstate the always-on block
          when the API serves integration/collected — see toMetaEntries.) */}
      {meta.length > 1 && (
        <div>
          <Kicker tone="soft" as="h2" className="mb-2">
            This frame
          </Kicker>
          <MetaList entries={meta} />
        </div>
      )}

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
        <div className="flex items-center justify-between mb-2">
          <Kicker tone="soft" as="h2">
            Sample tags
          </Kicker>
          {onManageTags && (
            <Button
              ref={manageTagsTriggerRef}
              variant="ghost"
              onClick={onManageTags}
              data-testid="manage-tags-trigger"
              className="px-2 py-0.5"
            >
              Manage
            </Button>
          )}
        </div>
        <TagList
          tags={tags}
          editable
          persistentAdd
          existingKeys={tags.map((t) => t.key)}
          {...(onAddTag ? { onAdd: onAddTag } : {})}
          {...(onRemoveTag ? { onRemoveById: onRemoveTag } : {})}
        />
      </div>

      {/* Block 5: Keys */}
      <div>
        <Kicker tone="soft" as="h2" className="mb-2">
          Keys
        </Kicker>
        <KbLegend shortcuts={shortcuts} className="flex-col gap-1" />
      </div>
    </aside>
  );
}
