import type { Meta, StoryObj } from "@storybook/react-vite";
import { SpecCell } from "./SpecCell";
import { KeptCell } from "./KeptCell";
import { StatusCell } from "./StatusCell";

// ---------------------------------------------------------------------------
// SpecCell
// ---------------------------------------------------------------------------

const specMeta = {
  title: "components/SpecCell",
  component: SpecCell,
  args: {
    name: "Sample A",
    sampleId: "s-001",
  },
} satisfies Meta<typeof SpecCell>;

export default specMeta;
type SpecStory = StoryObj<typeof specMeta>;

export const Screened: SpecStory = {
  args: { screened: true },
};

export const Unscreened: SpecStory = {
  args: { screened: false },
};

// ---------------------------------------------------------------------------
// KeptCell
// ---------------------------------------------------------------------------

const keptMeta = {
  title: "components/KeptCell",
  component: KeptCell,
  args: {
    kept: 4,
    total: 5,
  },
} satisfies Meta<typeof KeptCell>;

type KeptStory = StoryObj<typeof keptMeta>;

export const WithDropped: KeptStory = {
  args: { dropped: 1 },
};

export const AllKept: KeptStory = {
  args: { kept: 5, total: 5, dropped: 0 },
};

// ---------------------------------------------------------------------------
// StatusCell
// ---------------------------------------------------------------------------

const statusMeta = {
  title: "components/StatusCell",
  component: StatusCell,
} satisfies Meta<typeof StatusCell>;

type StatusStory = StoryObj<typeof statusMeta>;

export const IndexedPn3m: StatusStory = {
  args: { phase: "Pn3m" },
};

export const IndexedIm3m: StatusStory = {
  args: { phase: "Im3m" },
};

export const Unset: StatusStory = {
  args: { phase: null },
};
