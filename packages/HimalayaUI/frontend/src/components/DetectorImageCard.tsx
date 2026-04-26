import { useState } from "react";
import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposure: Exposure;
  onSetStatus: (status: "accepted" | "rejected" | null) => void;
  onSetIndexing: () => void;
  onAddTag: (key: string, value: string) => void;
}

export function DetectorImageCard({
  exposure,
  onSetStatus,
  onSetIndexing,
  onAddTag,
}: Props): JSX.Element {
  const [rejectMode, setRejectMode] = useState(false);
  const [rejectNote, setRejectNote] = useState("");

  const isRejected = exposure.status === "rejected";
  const isAccepted = exposure.status === "accepted";
  const isIndexing = exposure.selected;

  function handleRejectConfirm() {
    onSetStatus("rejected");
    if (rejectNote.trim()) {
      onAddTag("rejection_reason", rejectNote.trim());
    }
    setRejectMode(false);
    setRejectNote("");
  }

  const existingNote = exposure.tags.find(
    (t) => t.key === "rejection_reason",
  )?.value;

  return (
    <div className="card flex flex-col h-full min-h-0 p-3 gap-3">
      <div className="text-xs text-fg-muted truncate shrink-0">
        {exposure.filename ?? `Exposure #${exposure.id}`}
      </div>

      {/* Portrait image — takes available height */}
      <div className="flex-1 min-h-0 flex items-center justify-center">
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          size="full"
          className="max-h-full max-w-full object-contain"
        />
      </div>

      {/* Rejection note (read mode) */}
      {isRejected && existingNote && !rejectMode && (
        <p className="text-[10px] text-fg-muted italic shrink-0">
          {existingNote}{" "}
          <button
            className="underline"
            onClick={() => {
              setRejectNote(existingNote);
              setRejectMode(true);
            }}
          >
            edit
          </button>
        </p>
      )}

      {/* Rejection note input */}
      {rejectMode && (
        <div className="flex gap-1 shrink-0">
          <input
            className="flex-1 text-xs bg-bg border border-border rounded px-2 py-1"
            placeholder="Rejection reason (flare, missed sample…)"
            value={rejectNote}
            onChange={(e) => setRejectNote(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && handleRejectConfirm()}
            autoFocus
          />
          <button
            className="text-xs px-2 py-1 border border-red-400 text-red-400 rounded"
            onClick={handleRejectConfirm}
          >
            Confirm
          </button>
        </div>
      )}

      {/* Controls */}
      <div className="flex flex-col gap-1.5 shrink-0">
        <button
          onClick={() => onSetStatus(isAccepted ? null : "accepted")}
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            isAccepted
              ? "border-green-400 text-green-400 bg-green-400/10"
              : "border-border text-fg-muted hover:border-green-400 hover:text-green-400",
          ].join(" ")}
        >
          ✓ Accept
        </button>

        <button
          onClick={() => {
            if (isRejected) {
              onSetStatus(null);
            } else {
              setRejectMode(true);
            }
          }}
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            isRejected
              ? "border-red-400 text-red-400 bg-red-400/10"
              : "border-border text-fg-muted hover:border-red-400 hover:text-red-400",
          ].join(" ")}
        >
          ✗ Reject
        </button>

        <button
          disabled={isRejected}
          onClick={onSetIndexing}
          aria-label="Use for indexing"
          className={[
            "w-full text-xs py-1.5 rounded border transition-colors",
            "disabled:opacity-40 disabled:cursor-not-allowed",
            isIndexing
              ? "border-accent text-accent bg-accent/10"
              : "border-border text-fg-muted hover:border-accent hover:text-accent",
          ].join(" ")}
        >
          ⊙ Use for indexing
        </button>
      </div>
    </div>
  );
}
