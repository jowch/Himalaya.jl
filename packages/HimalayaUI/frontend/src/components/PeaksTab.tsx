import { useAppState } from "../state";
import { usePeaks } from "../queries";

function num(value: number | null, digits: number): string {
  return value != null ? value.toFixed(digits) : "—";
}

export function PeaksTab(): JSX.Element {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const q = usePeaks(activeExposureId);

  if (activeExposureId === undefined) {
    return <p className="text-fg-muted italic">No exposure selected.</p>;
  }
  if (q.isPending) return <p className="text-fg-muted">Loading peaks…</p>;
  if (q.error) {
    return <p className="text-error">Error: {(q.error as Error).message}</p>;
  }

  const peaks = q.data ?? [];
  if (peaks.length === 0) {
    return <p className="text-fg-muted italic">No peaks for this exposure.</p>;
  }

  return (
    <table className="w-full text-[13px] font-mono">
      <thead className="text-left text-fg-muted uppercase text-[11px] tracking-wide">
        <tr>
          <th className="py-1 pr-3">q</th>
          <th className="py-1 pr-3">prominence</th>
          <th className="py-1 pr-3">sharpness</th>
          <th className="py-1">source</th>
        </tr>
      </thead>
      <tbody>
        {peaks.map((p) => (
          <tr key={p.id} data-peak-id={p.id} className="border-t border-border">
            <td className="py-1 pr-3">{num(p.q, 4)}</td>
            <td className="py-1 pr-3">{num(p.prominence, 2)}</td>
            <td className="py-1 pr-3">{num(p.sharpness, 2)}</td>
            <td className="py-1">{p.source}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}
