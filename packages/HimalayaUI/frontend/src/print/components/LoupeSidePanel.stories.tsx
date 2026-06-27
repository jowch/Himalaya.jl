import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { LoupeSidePanel } from "./LoupeSidePanel";
import type { LoupeTag } from "../pages/loupeAdapters";

const mockMeta = [
  { key: "frame", value: "65" },
  { key: "exposure", value: "0.40 s" },
  { key: "detector", value: "Pilatus 1M" },
  { key: "energy", value: "12.4 keV" },
];

const mockTags: LoupeTag[] = [
  { id: 1, key: "LL37", source: "manual" },
  { id: 2, key: "temp", value: "37C", source: "manual" },
  { id: 3, key: "buffer", value: "PBS", source: "manual" },
];

const meta = {
  title: "components/LoupeSidePanel",
  component: LoupeSidePanel,
  args: {
    meta: mockMeta,
    dropped: false,
    isRepresentative: false,
    tags: mockTags,
  },
} satisfies Meta<typeof LoupeSidePanel>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Kept — the default verdict (not dropped). */
export const Kept: Story = {
  args: {
    dropped: false,
    isRepresentative: false,
    onToggleDrop: () => {},
    onSetRepresentative: () => {},
    onAddTag: () => {},
    onRemoveTag: () => {},
  },
  render: (args) => (
    <div className="bg-paper p-5 w-[300px]">
      <LoupeSidePanel {...args} />
    </div>
  ),
};

export const DroppedRepresentative: Story = {
  args: {
    dropped: true,
    isRepresentative: true,
    onToggleDrop: () => {},
    onSetRepresentative: () => {},
    onAddTag: () => {},
    onRemoveTag: () => {},
  },
  render: (args) => (
    <div className="bg-paper p-5 w-[300px]">
      <LoupeSidePanel {...args} />
    </div>
  ),
};

function InteractivePanel(): JSX.Element {
  // Binary verdict mirroring the LoupePage keyboard semantics:
  // X toggles the drop (cull ↔ un-cull).
  const [dropped, setDropped] = useState(false);
  const [isRep, setIsRep] = useState(false);
  const [tags, setTags] = useState<LoupeTag[]>(mockTags);
  return (
    <div className="bg-paper p-5 w-[300px]">
      <LoupeSidePanel
        meta={mockMeta}
        dropped={dropped}
        isRepresentative={isRep}
        tags={tags}
        onToggleDrop={() => setDropped((d) => !d)}
        onSetRepresentative={() => setIsRep(true)}
        onAddTag={(t) => setTags((prev) => [...prev, { id: Date.now(), key: t.key, source: "manual", ...(t.value !== undefined ? { value: t.value } : {}) }])}
        onRemoveTag={(id) => setTags((prev) => prev.filter((x) => x.id !== id))}
      />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <InteractivePanel />,
};
