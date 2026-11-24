import { cx } from "../lib/cx";

export type ReviewMark = "reviewed" | "open" | "none";

export const MAX_SEPARATE_SEGMENTS = 60;

const FILL: Record<ReviewMark, string> = {
  reviewed: "bg-success",
  open: "bg-hair-strong",
  none: "bg-hair",
};

export interface RibbonRun {
  mark: ReviewMark;
  count: number;
}

export function ribbonRuns(marks: readonly ReviewMark[]): RibbonRun[] {
  if (marks.length <= MAX_SEPARATE_SEGMENTS) return marks.map((mark) => ({ mark, count: 1 }));
  const runs: RibbonRun[] = [];
  for (const mark of marks) {
    const last = runs[runs.length - 1];
    if (last && last.mark === mark) last.count++;
    else runs.push({ mark, count: 1 });
  }
  return runs;
}

function samplesWord(n: number): string {
  return n === 1 ? "sample" : "samples";
}

export function SampleProgress({ marks }: { marks: readonly ReviewMark[] }): JSX.Element {
  const total = marks.filter((m) => m !== "none").length;
  const reviewed = marks.filter((m) => m === "reviewed").length;
  const shown = total > 0;
  const runs = ribbonRuns(marks);
  const separate = marks.length <= MAX_SEPARATE_SEGMENTS;
  return (
    <div className="flex flex-col gap-2">
      <div className="flex items-baseline justify-between gap-2">
        <h2 className="text-sm font-semibold text-ink-soft">Samples</h2>
        <span
          className={"text-xs text-ink-faint" + (shown ? "" : " invisible")}
          data-testid={shown ? "workspace-progress" : undefined}
          data-reserved={shown ? undefined : ""}
          aria-hidden={shown ? undefined : true}
        >
          {reviewed} of {total} reviewed
        </span>
      </div>
      <div
        data-testid="sample-ribbon"
        data-separate={separate ? "" : undefined}
        className={cx("flex h-2.5 w-full overflow-hidden bg-paper-sunk", separate && "gap-px")}
        role={shown ? "img" : undefined}
        aria-label={shown ? `${reviewed} of ${total} ${samplesWord(total)} reviewed` : undefined}
        aria-hidden={shown ? undefined : true}
      >
        {runs.map((r, i) => (
          <span key={i} data-mark={r.mark} data-count={r.count} className={cx("min-w-0", FILL[r.mark])} style={{ flexGrow: r.count, flexBasis: 0 }} />
        ))}
      </div>
    </div>
  );
}
