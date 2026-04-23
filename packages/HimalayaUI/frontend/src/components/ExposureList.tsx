import { useAppState } from "../state";
import { useExposures } from "../queries";

export function ExposureList(): JSX.Element {
  const activeSampleId   = useAppState((s) => s.activeSampleId);
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);

  const q = useExposures(activeSampleId);

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }

  if (q.isPending) return <p className="text-fg-muted">Loading exposures…</p>;
  if (q.error)     return <p className="text-error">Error: {(q.error as Error).message}</p>;

  const exposures = q.data ?? [];
  if (exposures.length === 0) {
    return <p className="text-fg-muted italic">No exposures for this sample.</p>;
  }

  return (
    <ul className="list-none flex flex-col gap-1">
      {exposures.map((e) => {
        const active = e.id === activeExposureId;
        return (
          <li
            key={e.id}
            data-exposure-id={e.id}
            data-active={active}
            className={
              "flex items-center gap-2 px-2 py-1 rounded-md cursor-pointer border-l-2 " +
              (active ? "bg-bg-elevated border-accent" : "border-transparent hover:bg-bg-hover")
            }
            onClick={() => setActiveExposure(e.id)}
          >
            <span className="font-mono text-[13px]">{e.filename ?? "(derived)"}</span>
            <span className="text-fg-muted text-[12px]">{e.kind}</span>
          </li>
        );
      })}
    </ul>
  );
}
