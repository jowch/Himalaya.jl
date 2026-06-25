import type { Meta, StoryObj } from "@storybook/react-vite";
import { SeriesCard } from "./SeriesCard";
import { FULL, TRANSITION } from "../waterfall/waterfall.fixtures";

const noop = () => {};

const meta: Meta<typeof SeriesCard> = {
  title: "components/SeriesCard",
  component: SeriesCard,
  parameters: { layout: "padded" },
};
export default meta;

type Story = StoryObj<typeof SeriesCard>;

export const Transition: Story = {
  render: () => (
    <div style={{ width: 360 }}>
      <SeriesCard
        rows={TRANSITION}
        segments={[
          { phase: "Ia3d" },
          { phase: "Im3m", coexistWith: ["Ia3d"] },
          { phase: "Im3m" },
        ]}
        figLabel="Fig. 1"
        title="LL37 titration of lipid 1-2"
        sampleCount={3}
        variable="LL37 : lipid ratio"
        provenance="SSRL · Apr 2026"
        editedLabel="2 days ago"
        author="JC"
        notice={{ tone: "new", count: 1 }}
        onClick={noop}
      />
    </div>
  ),
};

export const Full: Story = {
  render: () => (
    <div style={{ width: 360 }}>
      <SeriesCard
        rows={FULL}
        segments={[
          { phase: "Ia3d" },
          { phase: "Ia3d" },
          { phase: "Im3m", coexistWith: ["Ia3d"] },
          { phase: "Im3m" },
          { phase: "Lamellar" },
        ]}
        figLabel="Fig. 3"
        title="Lipid 1-3 concentration sweep"
        sampleCount={5}
        variable="lipid concentration"
        provenance="ALS · Jul 2025"
        editedLabel="1 week ago"
        author="MK"
        onClick={noop}
      />
    </div>
  ),
};

export const Draft: Story = {
  render: () => (
    <div style={{ width: 360 }}>
      <SeriesCard
        rows={TRANSITION}
        segments={[
          { phase: "Ia3d" },
          { phase: "Im3m", coexistWith: ["Ia3d"] },
          { phase: "Im3m" },
        ]}
        figLabel="Recipe"
        title="Lipid 1-2 dose response"
        sampleCount={3}
        variable="LL37 : lipid ratio"
        provenance="SSRL · Apr 2026"
        editedLabel="just now"
        author="JC"
        draft
        notice={{ tone: "draft" }}
        onClick={noop}
      />
    </div>
  ),
};

export const CrossExperiment: Story = {
  render: () => (
    <div style={{ width: 360 }}>
      <SeriesCard
        rows={FULL}
        segments={[
          { phase: "Pn3m" },
          { phase: "Pn3m" },
          { phase: "Im3m", coexistWith: ["Pn3m"] },
          { phase: "Im3m" },
          { phase: "Lamellar" },
        ]}
        figLabel="Fig. 5"
        title="Cross-experiment lipid comparison"
        sampleCount={5}
        variable="temperature"
        provenance={<span>&#8644; April + July &middot; q normalized</span>}
        editedLabel="3 days ago"
        author="JC"
        onClick={noop}
      />
    </div>
  ),
};
