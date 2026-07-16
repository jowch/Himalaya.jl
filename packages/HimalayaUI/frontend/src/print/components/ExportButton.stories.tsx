import type { Meta, StoryObj } from "@storybook/react-vite";
import { ExportButton } from "./ExportButton";

const noop = (): void => {};

const meta: Meta<typeof ExportButton> = {
  title: "components/ExportButton",
  component: ExportButton,
  args: {
    onCopy: noop,
    onDownloadPng: noop,
    onDownloadSvg: noop,
    ariaContext: "the LL37 series",
  },
};
export default meta;

type Story = StoryObj<typeof ExportButton>;

/** Default — Copy enabled; click the chevron to reveal Download PNG / SVG. */
export const Default: Story = {};

/** Clipboard unavailable (e.g. insecure origin) — Copy disabled, downloads remain. */
export const CopyDisabled: Story = { args: { copyDisabled: true } };

/** A render is in flight — every action is blocked. */
export const Pending: Story = { args: { pending: true, copyDisabled: true, pngDisabled: true } };

/** Page-level gate (data not ready) — fully disabled. */
export const FullyDisabled: Story = { args: { disabled: true, copyDisabled: true, pngDisabled: true } };
