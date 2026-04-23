import { useState } from "react";
import { useAppState } from "../state";
import { useSamples, useAddSampleTag, useRemoveSampleTag } from "../queries";

const EXPERIMENT_ID = 1;

export function TagsTab(): JSX.Element {
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const samplesQ = useSamples(EXPERIMENT_ID);
  const addTag    = useAddSampleTag(EXPERIMENT_ID, activeSampleId ?? 0);
  const removeTag = useRemoveSampleTag(EXPERIMENT_ID, activeSampleId ?? 0);

  const [key, setKey]     = useState("");
  const [value, setValue] = useState("");

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }

  const sample = samplesQ.data?.find((s) => s.id === activeSampleId);
  const tags = sample?.tags ?? [];

  const inputClass =
    "bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent";

  function onSubmit(ev: React.FormEvent): void {
    ev.preventDefault();
    const k = key.trim();
    const v = value.trim();
    if (!k || !v) return;
    addTag.mutate({ key: k, value: v }, {
      onSuccess: () => { setKey(""); setValue(""); },
    });
  }

  return (
    <div className="flex flex-col gap-3">
      <ul className="flex flex-wrap gap-1">
        {samplesQ.isPending ? (
          <li className="text-fg-muted italic text-[13px]">Loading…</li>
        ) : tags.length === 0 ? (
          <li className="text-fg-muted italic text-[13px]">No tags.</li>
        ) : tags.map((t) => (
          <li
            key={t.id}
            data-tag-id={t.id}
            className="flex items-center gap-1 px-2 py-1 bg-bg-elevated border border-border rounded-md text-[13px]"
          >
            <span className="font-mono">{t.key}</span>
            <span className="text-fg-muted">=</span>
            <span>{t.value}</span>
            {t.source === "manifest" && (
              <span className="text-fg-muted text-[11px] ml-1">(manifest)</span>
            )}
            <button
              className="ml-1 text-fg-muted hover:text-error rounded-md px-1 focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent"
              aria-label={`Remove ${t.key}`}
              onClick={() => removeTag.mutate(t.id)}
            >
              ×
            </button>
          </li>
        ))}
      </ul>

      <form onSubmit={onSubmit} className="flex gap-2 items-center">
        <input
          className={inputClass + " w-28"}
          type="text"
          placeholder="key"
          value={key}
          onChange={(e) => setKey(e.target.value)}
        />
        <input
          className={inputClass + " w-32"}
          type="text"
          placeholder="value"
          value={value}
          onChange={(e) => setValue(e.target.value)}
        />
        <button
          type="submit"
          className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
          aria-label="Add tag"
          disabled={addTag.isPending}
        >
          {addTag.isPending ? "Adding…" : "Add"}
        </button>
      </form>
    </div>
  );
}
