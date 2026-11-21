import type { TagColor } from "../../api";

export const TAG_SWATCH: Record<TagColor, { ink: string; wash: string; solid: string; name: string }> = {
  violet: { ink: "text-tag-violet-ink", wash: "bg-tag-violet-wash", solid: "bg-tag-violet-ink", name: "Violet" },
  blue: { ink: "text-tag-blue-ink", wash: "bg-tag-blue-wash", solid: "bg-tag-blue-ink", name: "Blue" },
  teal: { ink: "text-tag-teal-ink", wash: "bg-tag-teal-wash", solid: "bg-tag-teal-ink", name: "Teal" },
  rose: { ink: "text-tag-rose-ink", wash: "bg-tag-rose-wash", solid: "bg-tag-rose-ink", name: "Rose" },
  stone: { ink: "text-tag-stone-ink", wash: "bg-tag-stone-wash", solid: "bg-tag-stone-ink", name: "Stone" },
};
