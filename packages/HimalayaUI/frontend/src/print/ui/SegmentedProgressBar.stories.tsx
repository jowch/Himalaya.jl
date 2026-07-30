import type { Meta, StoryObj } from "@storybook/react-vite";
import { SegmentedProgressBar } from "./SegmentedProgressBar";
import { deriveSegments } from "../../lib/ingestStages";

const meta = {
  title: "ui/SegmentedProgressBar",
  component: SegmentedProgressBar,
  decorators: [
    (Story) => (
      <div className="w-72">
        <Story />
      </div>
    ),
  ],
} satisfies Meta<typeof SegmentedProgressBar>;

export default meta;
type Story = StoryObj<typeof meta>;

/** Mid-discovery: first segment filling, later stages untouched. */
export const Discovery: Story = {
  args: {
    segments: deriveSegments("discovery", 780, 1100)!,
    label: "Scan progress",
  },
};

/** Mid-analysis: discovery held full (dimmed accent), analyzing active. */
export const Analyzing: Story = {
  args: {
    segments: deriveSegments("analyzing", 92, 604)!,
    label: "Scan progress",
  },
};

/** Final stage — the shape the bar holds just before `ingest_complete`. */
export const Thumbnails: Story = {
  args: {
    segments: deriveSegments("thumbnails", 310, 604)!,
    label: "Scan progress",
  },
};

/** Clean rescan: nothing new to analyze, so the stage closes as 0-of-0 and reads
 *  COMPLETE rather than stalled at 0%. */
export const NothingToDo: Story = {
  args: {
    segments: deriveSegments("analyzing", 0, 0)!,
    label: "Scan progress",
  },
};
