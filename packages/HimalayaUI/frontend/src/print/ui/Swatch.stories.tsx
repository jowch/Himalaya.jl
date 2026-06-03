import type { Meta, StoryObj } from "@storybook/react-vite";
import { Swatch } from "./Swatch";

const meta = {
  title: "ui/Swatch",
  component: Swatch,
  args: { phase: "Pn3m" },
} satisfies Meta<typeof Swatch>;

export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};

const PHASES = ["Pn3m", "Im3m", "Ia3d", "Lamellar"];

export const AllShapes: Story = {
  render: () => (
    <div className="inline-flex flex-col gap-2">
      <div className="inline-flex items-center gap-2">
        {PHASES.map((p) => (
          <Swatch key={p} phase={p} shape="square" />
        ))}
      </div>
      <div className="inline-flex items-center gap-2">
        {PHASES.map((p) => (
          <Swatch key={p} phase={p} shape="circle" />
        ))}
      </div>
    </div>
  ),
};

// Two sizes: sm (9px, the default — reading-row dot / legend) and md (11px,
// the series-plot member swatch).
export const Sizes: Story = {
  render: () => (
    <div className="inline-flex items-center gap-4">
      <span className="inline-flex items-center gap-2">
        <Swatch phase="Pn3m" shape="circle" size="sm" />
        <span className="text-data text-ink-soft">sm · 9px</span>
      </span>
      <span className="inline-flex items-center gap-2">
        <Swatch phase="Pn3m" shape="circle" size="md" />
        <span className="text-data text-ink-soft">md · 11px</span>
      </span>
    </div>
  ),
};

// Coexistence: a 135° gradient blending two phases — the self-decoding swatch
// for a series-plot member with two phases present (e.g. Pn3m + Lamellar).
const COEXIST_PAIRS: Array<[string, string]> = [
  ["Pn3m", "Lamellar"],
  ["Im3m", "Hexagonal"],
  ["Ia3d", "Pn3m"],
];

export const Coexistence: Story = {
  render: () => (
    <div className="inline-flex flex-col gap-2">
      {COEXIST_PAIRS.map(([a, b]) => (
        <span key={`${a}-${b}`} className="inline-flex items-center gap-2">
          <Swatch phase={a} coexistWith={b} shape="circle" size="md" />
          <span className="text-data text-ink-soft">
            {a} + {b}
          </span>
        </span>
      ))}
    </div>
  ),
};

// Form factor: an empty (transparent + dashed) swatch — the member has no
// Bragg peaks, so there's no phase color to show.
export const FormFactor: Story = {
  render: () => (
    <div className="inline-flex items-center gap-4">
      <span className="inline-flex items-center gap-2">
        <Swatch phase="" empty shape="circle" size="md" />
        <span className="text-data text-ink-soft">circle</span>
      </span>
      <span className="inline-flex items-center gap-2">
        <Swatch phase="" empty shape="square" size="md" />
        <span className="text-data text-ink-soft">square</span>
      </span>
    </div>
  ),
};
