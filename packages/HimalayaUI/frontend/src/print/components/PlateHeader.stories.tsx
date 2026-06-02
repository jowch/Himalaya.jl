import type { Meta, StoryObj } from "@storybook/react-vite";
import { Button } from "../ui";
import { PlateHeader } from "./PlateHeader";

const meta = {
  title: "components/PlateHeader",
  component: PlateHeader,
  args: {
    kicker: "Integration",
    title: "Lipid 1‑2 + LL37 1:0.5",
    subtitle: "smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03",
  },
} satisfies Meta<typeof PlateHeader>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {
  args: {
    children: (
      <>
        <Button variant="ghost">Export</Button>
        <Button variant="ghost">Share</Button>
      </>
    ),
  },
};

export const NoKicker: Story = {
  args: {
    kicker: undefined,
    children: (
      <Button variant="ghost">Export</Button>
    ),
  },
};

export const NoSubtitle: Story = {
  args: {
    subtitle: undefined,
    children: (
      <Button variant="ghost">Export</Button>
    ),
  },
};

export const TitleOnly: Story = {
  args: {
    kicker: undefined,
    subtitle: undefined,
  },
};

export const AsH1: Story = {
  args: {
    as: "h1",
    children: (
      <Button variant="ghost">Export</Button>
    ),
  },
};
