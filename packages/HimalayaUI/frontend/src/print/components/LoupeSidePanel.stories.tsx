import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { LoupeSidePanel } from "./LoupeSidePanel";
import type { Tag } from "../ui";

const mockMeta = [
  { key: "frame", value: "65" },
  { key: "exposure", value: "0.40 s" },
  { key: "detector", value: "Pilatus 1M" },
  { key: "energy", value: "12.4 keV" },
];

const mockTags: Tag[] = [
  { key: "LL37" },
  { key: "temp", value: "37C" },
  { key: "buffer", value: "PBS" },
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
  const [dropped, setDropped] = useState(false);
  const [isRep, setIsRep] = useState(false);
  const [tags, setTags] = useState<Tag[]>(mockTags);
  return (
    <div className="bg-paper p-5 w-[300px]">
      <LoupeSidePanel
        meta={mockMeta}
        dropped={dropped}
        isRepresentative={isRep}
        tags={tags}
        onToggleDrop={() => setDropped((d) => !d)}
        onSetRepresentative={() => setIsRep(true)}
        onAddTag={(t) => setTags((prev) => [...prev, t])}
        onRemoveTag={(t) =>
          setTags((prev) =>
            prev.filter((x) => !(x.key === t.key && x.value === t.value)),
          )
        }
      />
    </div>
  );
}

export const Interactive: Story = {
  render: () => <InteractivePanel />,
};
