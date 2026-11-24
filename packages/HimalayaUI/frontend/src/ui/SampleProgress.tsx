export type ReviewMark = "reviewed" | "open" | "none";

const VIEW_W = 1000;
const VIEW_H = 10;

const FILL: Record<ReviewMark, string> = {
  reviewed: "fill-success",
  open: "fill-hair-strong",
  none: "fill-hair",
};

function samplesWord(n: number): string {
  return n === 1 ? "sample" : "samples";
}

export function SampleProgress({ marks }: { marks: readonly ReviewMark[] }): JSX.Element {
  const total = marks.filter((m) => m !== "none").length;
  const reviewed = marks.filter((m) => m === "reviewed").length;
  const width = VIEW_W / Math.max(marks.length, 1);
  const shown = total > 0;
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
      <svg
        className="block h-2.5 w-full"
        viewBox={`0 0 ${VIEW_W} ${VIEW_H}`}
        preserveAspectRatio="none"
        shapeRendering="crispEdges"
        role={shown ? "img" : undefined}
        aria-label={shown ? `${reviewed} of ${total} ${samplesWord(total)} reviewed` : undefined}
        aria-hidden={shown ? undefined : true}
      >
        <rect className="fill-paper-sunk" x={0} y={0} width={VIEW_W} height={VIEW_H} />
        {marks.map((m, i) => (
          <rect key={i} data-mark={m} className={FILL[m]} x={i * width} y={0} width={Math.max(width - 2, 1)} height={VIEW_H} />
        ))}
      </svg>
    </div>
  );
}
