import { definePreview } from "@storybook/react-vite";
import "../src/styles.css";
import { DetectorLutFilters } from "../src/print/detector";

export default definePreview({
  parameters: {
    layout: "centered",
  },
  // DetectorImage recolors via `filter: url(#detector-lut-…)`; mount the shared
  // SVG filter defs globally so detector stories render in colour (and a missing
  // #id can't blank the image). Mirrors the app-root mount in App.tsx.
  decorators: [
    (Story) => (
      <>
        <DetectorLutFilters />
        <Story />
      </>
    ),
  ],
});
