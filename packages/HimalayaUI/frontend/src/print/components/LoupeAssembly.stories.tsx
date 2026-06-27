import type { Meta, StoryObj } from "@storybook/react-vite";
import { useMemo, useState } from "react";
import { BigFrame } from "./BigFrame";
import { LoupeSidePanel } from "./LoupeSidePanel";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
import type { MetaEntry, Tag } from "../ui";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";
import thumb93 from "../fixtures/thumbs/93.png?url";

/**
 * Page simulation (NOT a component): assembles `BigFrame` + a filmstrip +
 * `LoupeSidePanel` into the loupe view as `sample-table.html` lays it out — a
 * 2-col `1fr / 286px` body with the big detector frame and the exposure strip
 * on the left, the aside on the right. The Layer-4 loupe page (the plate shell,
 * back button, serif head) is deferred; this story owns the cross-component
 * state (selected frame, binary verdicts, representative) the page will own.
 * Verdicts mirror the LoupePage keyboard semantics: X drops, K keeps.
 */

interface FrameDatum {
  id: number;
  src: string;
  frameNo: string;
  meta: MetaEntry[];
}

const FRAMES: FrameDatum[] = [
  { id: 65, src: thumb65, frameNo: "65", meta: [
    { key: "frame", value: "65" }, { key: "exposure", value: "0.40 s" },
    { key: "detector", value: "Pilatus 1M" }, { key: "energy", value: "12.4 keV" },
  ] },
  { id: 66, src: thumb66, frameNo: "66", meta: [
    { key: "frame", value: "66" }, { key: "exposure", value: "0.40 s" },
    { key: "detector", value: "Pilatus 1M" }, { key: "energy", value: "12.4 keV" },
  ] },
  { id: 67, src: thumb67, frameNo: "67", meta: [
    { key: "frame", value: "67" }, { key: "exposure", value: "0.40 s" },
    { key: "detector", value: "Pilatus 1M" }, { key: "energy", value: "12.4 keV" },
  ] },
  { id: 93, src: thumb93, frameNo: "93", meta: [
    { key: "frame", value: "93" }, { key: "exposure", value: "0.80 s" },
    { key: "detector", value: "Pilatus 1M" }, { key: "energy", value: "12.4 keV" },
  ] },
];

const TAGS: Tag[] = [{ key: "LL37" }, { key: "temp", value: "37C" }];

function LoupeView(): JSX.Element {
  const [selectedId, setSelectedId] = useState(65);
  const [dropped, setDropped] = useState<Record<number, boolean>>({ 93: true });
  const [repId, setRepId] = useState<number | null>(65);

  const current = FRAMES.find((f) => f.id === selectedId) ?? FRAMES[0]!;
  const isDropped = dropped[current.id] ?? false;
  const verdictWord = isDropped ? "dropped" : "kept";

  const exposures = useMemo<GalleryExposure[]>(
    () => FRAMES.map((f) => ({
      id: f.id,
      src: f.src,
      frameNo: f.frameNo,
      rejected: dropped[f.id] ?? false,
      representative: repId === f.id,
    })),
    [dropped, repId],
  );

  // Binary verdict: Drop is a toggle (cull ↔ un-cull).
  const onToggleDrop = (): void =>
    setDropped((prev) => ({ ...prev, [current.id]: !(prev[current.id] ?? false) }));

  return (
    <div className="bg-paper p-7">
      <div className="grid grid-cols-[1fr_286px] gap-7 max-w-[1100px] mx-auto">
        <div>
          <BigFrame
            src={current.src}
            caption={`frame ${current.frameNo} · ${verdictWord}`}
            rejected={isDropped}
          />
          <ThumbnailGallery
            exposures={exposures}
            selectedId={selectedId}
            onSelect={setSelectedId}
            size="lg"
            align="center"
            className="mt-3"
          />
        </div>
        <LoupeSidePanel
          meta={current.meta}
          dropped={isDropped}
          isRepresentative={repId === current.id}
          tags={TAGS}
          onToggleDrop={onToggleDrop}
          onSetRepresentative={() => setRepId(current.id)}
          onAddTag={() => {}}
          onRemoveTag={() => {}}
        />
      </div>
    </div>
  );
}

const meta = {
  title: "components/LoupeAssembly",
  component: LoupeView,
} satisfies Meta<typeof LoupeView>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Loupe: Story = {};
