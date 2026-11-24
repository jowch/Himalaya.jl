import type { Meta, StoryObj } from "@storybook/react-vite";
import { CollectionRibbon } from "./CollectionRibbon";
import type { CollectionTimeline, ReviewProgress } from "../api";

function beamtime(samples: number, reviewedShare: number, nightAt: number | null): { timeline: CollectionTimeline; review: ReviewProgress } {
  const rows: CollectionTimeline["samples"] = [];
  let t = 0;
  for (let i = 0; i < samples; i++) {
    if (i > 0 && i % 13 === 0) t += 1500;
    if (nightAt !== null && i === nightAt) t += 8 * 3600;
    rows.push([t, t + 520, i < samples * reviewedShare ? 1 : 0]);
    t += 540;
  }
  const timeline = { start: "2026-04-25T17:40:00", samples: rows, beam_off: nightAt === null ? [] : [rows[40]![1]] };
  return { timeline, review: { reviewed: rows.filter(([, , r]) => r === 1).length, total: samples } };
}

const meta = {
  title: "ui/CollectionRibbon",
  component: CollectionRibbon,
  decorators: [(Story) => <div className="w-96 p-4 bg-plate"><Story /></div>],
} satisfies Meta<typeof CollectionRibbon>;

export default meta;
type Story = StoryObj<typeof meta>;

export const PartlyReviewed: Story = { args: beamtime(275, 0.6, 120) };
export const AllReviewed: Story = { args: beamtime(300, 1, null) };
export const NothingReviewed: Story = { args: beamtime(13, 0, null) };
export const NoImageTimes: Story = { args: { timeline: null, review: { reviewed: 0, total: 4 } } };
