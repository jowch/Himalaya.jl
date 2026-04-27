import { useEffect, useState } from "react";
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

  // Reset rejection note state when switching to a different exposure
  useEffect(() => {
    setRejectMode(false);
    setRejectNote("");
  }, [exposure.id]);

  const isRejected = exposure.status === "rejected";
  const isAccepted = exposure.status === "accepted";
  const isIndexing = exposure.selected;

  const existingNote = exposure.tags.find(
    (t) => t.key === "rejection_reason",
  )?.value;

  function handleRejectPillClick() {
    if (isRejected) {
      onSetStatus(null);
    } else {
      setRejectMode(true);
    }
  }

  function handleRejectConfirm() {
    onSetStatus("rejected");
    if (rejectNote.trim()) {
      onAddTag("rejection_reason", rejectNote.trim());
    }
    setRejectMode(false);
    setRejectNote("");
  }

  const filename = exposure.filename ?? `Exposure #${exposure.id}`;

  // Status header strip
  const headerCls = isAccepted
    ? "bg-success/8 border-b border-success/20"
    : isRejected
      ? "bg-error/8 border-b border-error/20"
      : "border-b border-border";

  const barCls = isAccepted
    ? "bg-success"
    : isRejected
      ? "bg-error"
      : "bg-transparent";

  // Segmented pill classes
  const pillBase =
    "px-3 py-1.5 text-xs transition-colors border-r border-border last:border-r-0";

  const pendingCls = [
    pillBase,
    !isAccepted && !isRejected
      ? "bg-bg-subtle text-fg"
      : "text-fg-muted hover:text-fg hover:bg-bg-hover",
  ].join(" ");

  const acceptCls = [
    pillBase,
    isAccepted
      ? "bg-success/15 text-success"
      : "text-fg-muted hover:text-fg hover:bg-bg-hover",
  ].join(" ");

  const rejectCls = [
    pillBase,
    isRejected
      ? "bg-error/15 text-error"
      : "text-fg-muted hover:text-fg hover:bg-bg-hover",
  ].join(" ");

  return (
    <div className="flex flex-col h-full min-h-0">
      {/* Status header strip */}
      <div className={`flex items-center gap-2 shrink-0 ${headerCls}`}>
        <div className={`w-1 self-stretch shrink-0 ${barCls}`} />
        <div className="flex items-center gap-2 py-2 pr-3 min-w-0 flex-1">
          {isAccepted && (
            <span className="text-[10px] font-medium text-success shrink-0">
              ✓ Accepted
            </span>
          )}
          {isRejected && (
            <span className="text-[10px] font-medium text-error shrink-0">
              ✗ Rejected
            </span>
          )}
          {isRejected && existingNote && !rejectMode && (
            <span className="text-[10px] text-fg-muted italic truncate">
              {existingNote}{" "}
              <button
                className="underline shrink-0"
                onClick={() => {
                  setRejectNote(existingNote);
                  setRejectMode(true);
                }}
              >
                edit
              </button>
            </span>
          )}
          <span className="text-[10px] text-fg-muted truncate ml-auto">
            {filename}
          </span>
        </div>
      </div>

      {/* Detector image — takes all available height */}
      <div className="flex-1 min-h-0 flex items-center justify-center p-2">
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          size="full"
          className="max-h-full max-w-full object-contain"
        />
      </div>

      {/* Rejection note input (only in rejectMode) */}
      {rejectMode && (
        <div className="flex gap-1 shrink-0 px-3 pb-2">
          <input
            className="flex-1 text-xs bg-bg border border-border rounded px-2 py-1"
            placeholder="Rejection reason (flare, missed sample…)"
            value={rejectNote}
            onChange={(e) => setRejectNote(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && handleRejectConfirm()}
            autoFocus
          />
          <button
            className="text-xs px-2 py-1 border border-error text-error rounded"
            onClick={handleRejectConfirm}
          >
            Confirm
          </button>
        </div>
      )}

      {/* Control row */}
      <div className="flex items-center gap-2 shrink-0 px-3 pb-3">
        {/* Segmented Accept/Reject/Pending control */}
        <div className="flex rounded-md overflow-hidden border border-border">
          <button className={pendingCls} onClick={() => onSetStatus(null)}>
            ○ Pending
          </button>
          <button
            className={acceptCls}
            onClick={() => onSetStatus(isAccepted ? null : "accepted")}
          >
            ✓ Accept
          </button>
          <button className={rejectCls} onClick={handleRejectPillClick}>
            ✗ Reject
          </button>
        </div>

        {/* Indexing button */}
        <button
          disabled={isRejected}
          onClick={onSetIndexing}
          aria-label="Use for indexing"
          className={[
            "ml-auto text-xs px-3 py-1.5 rounded border transition-colors",
            "disabled:opacity-40 disabled:cursor-not-allowed",
            isIndexing
              ? "border-accent text-accent bg-accent/10"
              : "border-border text-fg-muted hover:border-accent hover:text-accent",
          ].join(" ")}
        >
          ⊙ Index
        </button>
      </div>
    </div>
  );
}
