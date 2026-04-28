import { useEffect, useState } from "react";
import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposure: Exposure;
  onSetStatus: (status: "accepted" | "rejected" | null) => void;
  onSetIndexing: () => void;
  onAddTag: (key: string, value: string) => void;
}

type RejectStep = "idle" | "picking" | "custom";

const BUILT_IN_REJECT_REASONS = ["Flare"] as const;

export function DetectorImageCard({
  exposure,
  onSetStatus,
  onSetIndexing,
  onAddTag,
}: Props): JSX.Element {
  const [rejectStep, setRejectStep] = useState<RejectStep>("idle");
  const [customNote, setCustomNote] = useState("");

  // Reset rejection state when switching to a different exposure
  useEffect(() => {
    setRejectStep("idle");
    setCustomNote("");
  }, [exposure.id]);

  const isRejected = exposure.status === "rejected";
  const isAccepted = exposure.status === "accepted";
  const isIndexing = exposure.selected;

  const existingNote = exposure.tags.find(
    (t) => t.key === "rejection_reason",
  )?.value;

  function submitReject(note: string) {
    onSetStatus("rejected");
    const trimmed = note.trim();
    if (trimmed) onAddTag("rejection_reason", trimmed);
    setRejectStep("idle");
    setCustomNote("");
  }

  function handleRejectPillClick() {
    if (isRejected) {
      // Un-reject (preserves existing toggle behavior)
      onSetStatus(null);
    } else {
      setRejectStep("picking");
    }
  }

  function handleEditExisting() {
    setCustomNote(existingNote ?? "");
    setRejectStep("custom");
  }

  function cancelReject() {
    setRejectStep("idle");
    setCustomNote("");
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

  // Segmented pill classes (idle state)
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

  // Reject-flow chip styling — visually distinct from idle pills
  const chipBase =
    "px-3 py-1.5 text-xs transition-colors border-r border-border last:border-r-0";
  const cancelChipCls = `${chipBase} text-fg-muted hover:text-fg hover:bg-bg-hover`;
  const reasonChipCls = `${chipBase} text-error hover:bg-error/10`;

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
          {isRejected && existingNote && rejectStep === "idle" && (
            <span className="text-[10px] text-fg-muted italic truncate">
              {existingNote}{" "}
              <button
                className="underline shrink-0"
                onClick={handleEditExisting}
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
          imageVersion={exposure.image_version}
          size="full"
          className="max-h-full max-w-full object-contain"
        />
      </div>

      {/* Control row — morphs based on rejectStep */}
      <div className="flex items-center gap-2 shrink-0 px-3 pb-3">
        {rejectStep === "idle" && (
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
        )}

        {rejectStep === "picking" && (
          <div className="flex rounded-md overflow-hidden border border-error/40">
            <button className={cancelChipCls} onClick={cancelReject}>
              ← Cancel
            </button>
            {BUILT_IN_REJECT_REASONS.map((reason) => (
              <button
                key={reason}
                className={reasonChipCls}
                onClick={() => submitReject(reason)}
              >
                {reason}
              </button>
            ))}
            <button
              className={reasonChipCls}
              onClick={() => setRejectStep("custom")}
            >
              Other
            </button>
          </div>
        )}

        {rejectStep === "custom" && (
          <div className="flex items-stretch rounded-md overflow-hidden border border-error/40">
            <button className={cancelChipCls} onClick={cancelReject}>
              ← Cancel
            </button>
            <input
              className="text-xs bg-bg px-2 py-1 border-r border-border min-w-[180px] focus:outline-none"
              placeholder="Rejection reason…"
              value={customNote}
              onChange={(e) => setCustomNote(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") submitReject(customNote);
                if (e.key === "Escape") cancelReject();
              }}
              autoFocus
            />
            <button
              className={reasonChipCls}
              onClick={() => submitReject(customNote)}
            >
              Confirm
            </button>
          </div>
        )}

        {/* Indexing button — always visible on the right */}
        <button
          disabled={isRejected || rejectStep !== "idle"}
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
