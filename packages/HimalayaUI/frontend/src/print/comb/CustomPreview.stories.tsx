import type { Meta, StoryObj } from "@storybook/react-vite";
import { CustomPreview } from "./CustomPreview";
import { PN3M, IM3M } from "./comb.fixtures";
import type { CombSeries } from "./combModel";

const meta: Meta<typeof CustomPreview> = {
  title: "comb/CustomPreview",
  component: CustomPreview,
};
export default meta;
type Story = StoryObj<typeof CustomPreview>;

const frame = "rounded border border-hair bg-plate p-4";

const observedOf = (s: CombSeries): number[] =>
  s.teeth.filter((t) => t.observed).map((t) => t.q);

const WITH_ABSENT: CombSeries = {
  ...PN3M,
  teeth: PN3M.teeth.map((t, i) => (i === 2 ? { ...t, observed: false } : t)),
};

export const Default: Story = {
  render: () => (
    <div className={frame} style={{ width: 540 }}>
      <CustomPreview series={PN3M} observed={observedOf(PN3M)} />
    </div>
  ),
};

export const Cubic: Story = {
  render: () => (
    <div className={frame} style={{ width: 540 }}>
      <CustomPreview series={IM3M} observed={observedOf(IM3M)} />
    </div>
  ),
};

export const WithAbsent: Story = {
  render: () => (
    <div className={frame} style={{ width: 540 }}>
      <CustomPreview series={WITH_ABSENT} observed={observedOf(WITH_ABSENT)} />
    </div>
  ),
};
