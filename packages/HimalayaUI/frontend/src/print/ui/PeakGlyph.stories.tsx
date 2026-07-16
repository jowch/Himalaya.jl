import type { Meta, StoryObj } from "@storybook/react-vite";
// PeakGlyph is intentionally NOT in the barrel index.ts — import directly.
import { PeakGlyph } from "./PeakGlyph";
import { peakGlyph } from "./peakMark";

// PeakGlyph paints into a parent SVG coordinate space, so each story wraps it
// in a small SVG and centres the glyph. The descriptor is built from the §5.1
// atoms via peakGlyph(); colour is a resolved string (callers never pass a
// CSS-var READ to the builder).
const PHASE_COLOR = "#c2603f"; // resolved terracotta accent stand-in

function Frame({
  children,
}: {
  children: React.ReactNode;
}): JSX.Element {
  return (
    <svg width={48} height={48} viewBox="0 0 48 48">
      {children}
    </svg>
  );
}

const meta = {
  title: "ui/PeakGlyph",
  component: PeakGlyph,
  render: (args) => (
    <Frame>
      <PeakGlyph {...args} />
    </Frame>
  ),
  args: {
    x: 24,
    y: 32,
    descriptor: peakGlyph({ source: "auto", color: PHASE_COLOR }),
  },
} satisfies Meta<typeof PeakGlyph>;

export default meta;
type Story = StoryObj<typeof meta>;

// auto → filled downward triangle (points at the peak).
export const TriangleAuto: Story = {};

// manual → diamond silhouette.
export const Diamond: Story = {
  args: { descriptor: peakGlyph({ source: "manual", color: PHASE_COLOR }) },
};

// predicted-absent → hollow caret.
export const Caret: Story = {
  args: {
    descriptor: peakGlyph({
      source: "auto",
      predictedAbsent: true,
      color: PHASE_COLOR,
    }),
  },
};

// hot (q-link) → grown glyph emphasized with a darker, heavier accent OUTLINE
// stroke on the glyph itself (no surrounding ring).
export const Hot: Story = {
  args: {
    descriptor: peakGlyph({ source: "auto", hot: true, color: PHASE_COLOR }),
  },
};

// excluded → ghosted (hollow) auto triangle.
export const Excluded: Story = {
  args: {
    descriptor: peakGlyph({ source: "auto", excluded: true, color: PHASE_COLOR }),
  },
};
