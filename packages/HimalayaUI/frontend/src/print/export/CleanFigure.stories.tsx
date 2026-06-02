import type { Meta, StoryObj } from "@storybook/react-vite";
import { CleanFigure } from "./CleanFigure";
import { TRANSITION, FULL } from "../waterfall/waterfall.fixtures";

const meta: Meta<typeof CleanFigure> = {
  title: "export/CleanFigure",
  component: CleanFigure,
};
export default meta;
type Story = StoryObj<typeof CleanFigure>;

export const Default: Story = {
  render: () => (
    <div style={{ width: 600 }}>
      <CleanFigure
        rows={TRANSITION}
        title="LL37 Titration: Pn3m → Lamellar"
        footer="Lipid 1-2 · SSRL Apr 2026 · q normalized · intensity offset for clarity"
      />
    </div>
  ),
};

export const FullSeries: Story = {
  render: () => (
    <div style={{ width: 600 }}>
      <CleanFigure
        rows={FULL}
        title="Sample 9 Series: Ia3d → Im3m"
        footer="Lipid 1-2 · SSRL Apr 2026 · q normalized · intensity offset for clarity"
      />
    </div>
  ),
};
