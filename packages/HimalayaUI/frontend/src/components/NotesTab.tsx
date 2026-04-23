import { useEffect, useRef, useState } from "react";
import { useAppState } from "../state";
import { useSamples, useUpdateSample } from "../queries";

const EXPERIMENT_ID = 1;

export function NotesTab(): JSX.Element {
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const samplesQ = useSamples(EXPERIMENT_ID);
  const updateSample = useUpdateSample(EXPERIMENT_ID, activeSampleId ?? 0);

  const sample = samplesQ.data?.find((s) => s.id === activeSampleId);
  const serverNotes = sample?.notes ?? "";

  const [value, setValue] = useState(serverNotes);
  const lastServerRef = useRef(serverNotes);

  useEffect(() => {
    if (serverNotes !== lastServerRef.current) {
      setValue(serverNotes);
      lastServerRef.current = serverNotes;
    }
  }, [serverNotes]);

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }

  function onBlur(): void {
    if (value === lastServerRef.current) return;
    updateSample.mutate({ notes: value });
    lastServerRef.current = value;
  }

  return (
    <textarea
      aria-label="Sample notes"
      className="w-full h-full min-h-24 bg-bg-elevated border border-border rounded-md p-2 text-[13px] focus:outline focus:outline-1 focus:outline-accent focus:border-accent resize-none"
      value={value}
      onChange={(e) => setValue(e.target.value)}
      onBlur={onBlur}
      placeholder="Notes about this sample…"
    />
  );
}
