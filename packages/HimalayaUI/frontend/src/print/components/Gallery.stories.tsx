import type { Meta, StoryObj } from "@storybook/react-vite";
import { Gallery } from "./Gallery";
import { SeriesCard } from "./SeriesCard";
import { EmptyState } from "../ui";
import { FULL, TRANSITION } from "../waterfall/waterfall.fixtures";

const noop = () => {};

const meta: Meta<typeof Gallery> = {
  title: "components/Gallery",
  component: Gallery,
  parameters: { layout: "fullscreen" },
};
export default meta;

type Story = StoryObj<typeof Gallery>;

// ---------------------------------------------------------------------------
// Wall — masonry of 5 cards, varying row counts (3 TRANSITION + 2 FULL)
// so distinct heights are visible in the multi-column layout.
// ---------------------------------------------------------------------------

export const Wall: Story = {
  render: () => (
    <div className="p-8">
      <Gallery>
        {/* Card 1: 3-row TRANSITION — shortest */}
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
          notice={{ tone: "new", count: 2 }}
          onClick={noop}
        />

        {/* Card 2: 5-row FULL — tallest */}
        <SeriesCard
          rows={FULL}
          segments={[
            { phase: "Ia3d" },
            { phase: "Ia3d" },
            { phase: "Im3m", coexistWith: ["Ia3d"] },
            { phase: "Im3m" },
            { phase: "Lamellar" },
          ]}
          figLabel="Fig. 2"
          title="Lipid 1-3 concentration sweep"
          sampleCount={5}
          variable="lipid concentration"
          provenance="ALS · Jul 2025"
          editedLabel="1 week ago"
          author="MK"
          onClick={noop}
        />

        {/* Card 3: 3-row TRANSITION */}
        <SeriesCard
          rows={TRANSITION}
          segments={[
            { phase: "Pn3m" },
            { phase: "Pn3m", coexistWith: ["Im3m"] },
            { phase: "Im3m" },
          ]}
          figLabel="Fig. 3"
          title="Pn3m → Im3m temperature ramp"
          sampleCount={3}
          variable="temperature (°C)"
          provenance="SSRL · Mar 2026"
          editedLabel="5 days ago"
          author="JC"
          draft
          notice={{ tone: "draft" }}
          onClick={noop}
        />

        {/* Card 4: 5-row FULL */}
        <SeriesCard
          rows={FULL}
          segments={[
            { phase: "Pn3m" },
            { phase: "Pn3m" },
            { phase: "Ia3d", coexistWith: ["Pn3m"] },
            { phase: "Ia3d" },
            { phase: "Lamellar" },
          ]}
          figLabel="Fig. 4"
          title="Cross-experiment lipid comparison"
          sampleCount={5}
          variable="temperature"
          provenance={<span>&#8644; April + July &middot; q normalized</span>}
          editedLabel="3 days ago"
          author="JC"
          onClick={noop}
        />

        {/* Card 5: 3-row TRANSITION */}
        <SeriesCard
          rows={TRANSITION}
          segments={[
            { phase: "Hexagonal" },
            { phase: "Hexagonal", coexistWith: ["Ia3d"] },
            { phase: "Ia3d" },
          ]}
          figLabel="Fig. 5"
          title="Hexagonal → cubic phase boundary"
          sampleCount={3}
          variable="monoolein fraction"
          provenance="ALS · Jun 2026"
          editedLabel="just now"
          author="MK"
          notice={{ tone: "new", count: 1 }}
          onClick={noop}
        />
      </Gallery>
    </div>
  ),
};

// ---------------------------------------------------------------------------
// Empty — no children, shows the empty-state node.
// ---------------------------------------------------------------------------

export const Empty: Story = {
  render: () => (
    <div className="p-8">
      <Gallery
        empty={
          <EmptyState
            title="No series match"
            body="Adjust the search or filter to find a series."
          />
        }
      >
        {/* intentionally empty */}
      </Gallery>
    </div>
  ),
};
