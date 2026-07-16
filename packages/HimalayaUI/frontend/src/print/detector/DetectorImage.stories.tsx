import type { Meta, StoryObj } from "@storybook/react-vite";
import { DetectorImage } from "./DetectorImage";
import { buildPrintDetectorLut } from "./detectorLut";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof DetectorImage> = { title: "detector/DetectorImage", component: DetectorImage };
export default meta;
type Story = StoryObj<typeof DetectorImage>;

// Real detector geometry for fixture 37: 981×1043 (portrait). The frame matches
// that aspect so object-fit:contain fills it edge-to-edge (no letterbox).
const frame = "aspect-[981/1043] w-[420px] overflow-hidden rounded border border-frame-edge bg-frame-edge";

export const Full: Story = {
  render: () => <div className={frame}><DetectorImage src={thumb37} size="full" className="h-full w-full" /></div>,
};
export const Thumb: Story = {
  render: () => (
    <div className="h-8 w-8 overflow-hidden rounded-sm border border-frame-edge bg-frame-edge">
      <DetectorImage src={thumb37} size="thumb" className="h-full w-full" />
    </div>
  ),
};
export const MissingImage: Story = {
  render: () => <div className={frame}><DetectorImage src={null} size="full" className="h-full w-full" /></div>,
};

// The saturated "beauty" ramp — kept available but NOT the default, since warm
// data clashes with the (often warm-hued) phase rings.
export const WarmExposure: Story = {
  render: () => <div className={frame}><DetectorImage src={thumb37} size="full" lutVariant="warm" className="h-full w-full" /></div>,
};

// Both colormaps as horizontal swatches — neutral (default) above, warm below.
// The detector image is raw measurement, so the default ramp stays near-neutral
// and the phase rings carry the colour.
export const Colormap: Story = {
  render: () => {
    const neutral = buildPrintDetectorLut();
    const warm = buildPrintDetectorLut("warm");
    const swatch = (lut: Uint8ClampedArray, label: string) => (
      <div>
        <div className="text-meta text-ink-faint mb-1">{label}</div>
        <div className="flex h-10 w-[512px] overflow-hidden rounded border border-hair">
          {Array.from({ length: 256 }, (_, i) => (
            <div key={i} className="h-full flex-1"
                 style={{ background: `rgb(${lut[i*3]}, ${lut[i*3+1]}, ${lut[i*3+2]})` }} />
          ))}
        </div>
      </div>
    );
    return (
      <div className="flex flex-col gap-3">
        {swatch(neutral, "neutral (default — image stays a warm gray)")}
        {swatch(warm, "warm (optional beauty exposure)")}
      </div>
    );
  },
};
