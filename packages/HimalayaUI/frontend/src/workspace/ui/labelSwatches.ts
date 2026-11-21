import type { LabelColor } from "../../api";

export const LABEL_SWATCH: Record<LabelColor, { ink: string; wash: string; solid: string; name: string }> = {
  violet: { ink: "text-label-violet-ink", wash: "bg-label-violet-wash", solid: "bg-label-violet-ink", name: "Violet" },
  blue: { ink: "text-label-blue-ink", wash: "bg-label-blue-wash", solid: "bg-label-blue-ink", name: "Blue" },
  teal: { ink: "text-label-teal-ink", wash: "bg-label-teal-wash", solid: "bg-label-teal-ink", name: "Teal" },
  rose: { ink: "text-label-rose-ink", wash: "bg-label-rose-wash", solid: "bg-label-rose-ink", name: "Rose" },
  stone: { ink: "text-label-stone-ink", wash: "bg-label-stone-wash", solid: "bg-label-stone-ink", name: "Stone" },
};
