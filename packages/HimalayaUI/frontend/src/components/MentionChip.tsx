import { useState, useCallback } from "react";
import { useAppState } from "../state";
import { useExperiment } from "../queries";
import { latticeUnitFromQUnits, inverseSquareUnits, formatKappa } from "../lib/units";
import type { ResolvedMention } from "../hooks/useMentionResolution";
import { CUBIC_PHASES } from "../phases";

interface ChipProps {
  resolved: ResolvedMention | "loading" | "dead";
  originalText?: string;
  /**
   * Phase 10: only meaningful for `comparison` mentions. The 8-char hash
   * embedded in the source token (`[[comparison:N@hhhhhhhh]]`). When present
   * and divergent from the resolved comparison's `content_hash[:8]`, the
   * chip shows a "(changed)" annotation and `data-hash-drift="true"` for
   * E2E selectors. When omitted, the drift indicator is suppressed.
   */
  tokenHash?: string;
}

const CHIP_STYLES: Record<string, string> = {
  peak:       "text-[#7cb8e8] bg-[#15222e] border-[#4a7aaa]",
  index:      "text-[#b5a0d8] bg-[#1e1828] border-[#7a60a8]",
  exposure:   "text-[#88c0a8] bg-[#162018] border-[#508070]",
  sample:     "text-[#c0b878] bg-[#201c10] border-[#887840]",
  experiment: "text-[#c0b878] bg-[#201c10] border-[#887840]",
  comparison: "text-[#e0a878] bg-[#241c14] border-[#a07848]",
  dead:       "text-[#484848] bg-[#181818] border-[#333333]",
  loading:    "text-[#484848] bg-[#181818] border-[#333333]",
};

const CHIP_HOVER_STYLES: Record<string, string> = {
  peak:       "hover:text-[#a8d4f8] hover:bg-[#1c3045] hover:border-[#7cb8e8]",
  index:      "hover:text-[#cdb8f0] hover:bg-[#2a2040] hover:border-[#b5a0d8]",
  exposure:   "hover:text-[#a8e0c0] hover:bg-[#1c3028] hover:border-[#88c0a8]",
  sample:     "hover:text-[#dcd090] hover:bg-[#282215] hover:border-[#c0b878]",
  experiment: "hover:text-[#dcd090] hover:bg-[#282215] hover:border-[#c0b878]",
  comparison: "hover:text-[#f4c498] hover:bg-[#332518] hover:border-[#e0a878]",
  dead:       "hover:text-[#606060] hover:bg-[#1e1e1e] hover:border-[#444444]",
  loading:    "hover:text-[#606060] hover:bg-[#1e1e1e] hover:border-[#444444]",
};

function chipLabel(resolved: ResolvedMention): string {
  switch (resolved.type) {
    case "peak":       return `q = ${resolved.data.q.toFixed(3)}`;
    case "index":      return `${resolved.data.phase} · ${(resolved.data.score ?? 0).toFixed(2)}`;
    case "exposure":   return resolved.data.filename ?? `exposure ${resolved.data.id}`;
    case "sample":     return resolved.data.name ?? resolved.data.label ?? `sample ${resolved.data.id}`;
    case "experiment": return resolved.data.name ?? `experiment ${resolved.data.id}`;
    case "comparison": return resolved.data.title;
  }
}

interface TooltipProps {
  resolved: ResolvedMention;
  latticeUnit: string;
  curvatureUnit: string;
  hashDrift: boolean;
}

function TooltipContent({ resolved, latticeUnit, curvatureUnit, hashDrift }: TooltipProps): JSX.Element | null {
  switch (resolved.type) {
    case "peak":
      return (
        <span>
          source <code>{resolved.data.source}</code>
          {resolved.data.prominence != null && (
            <> · prominence <code>{resolved.data.prominence.toFixed(1)}</code></>
          )}
        </span>
      );
    case "index": {
      const q1      = resolved.data.predicted_q[0];
      const isCubic = CUBIC_PHASES.has(resolved.data.phase);
      const ngc     = resolved.data.ngc;
      return (
        <span>
          {q1 != null && <>q₁ <code>{q1.toFixed(3)}</code> · </>}
          {resolved.data.lattice_d != null && (
            <>d <code>{resolved.data.lattice_d.toFixed(2)} {latticeUnit}</code> · </>
          )}
          R² <code>{(resolved.data.r_squared ?? 0).toFixed(3)}</code>
          {isCubic && ngc != null && (
            <> · κ <code>{formatKappa(ngc)} {curvatureUnit}</code></>
          )}
        </span>
      );
    }
    case "exposure":
      return resolved.data.status != null
        ? <span>status <code>{resolved.data.status}</code></span>
        : null;
    case "sample":
    case "experiment":
      return null;
    case "comparison": {
      const memberCount = resolved.data.members.length;
      return (
        <span>
          {memberCount} member{memberCount === 1 ? "" : "s"}
          {hashDrift && (
            <> · <span className="text-[#cca888]">
              this comparison has changed since the citation was made
            </span></>
          )}
        </span>
      );
    }
  }
}

/**
 * Compute the hash drift state for comparison mentions.
 *
 * Returns `false` (no drift) when:
 *   - the mention is not a comparison
 *   - the source token did not carry a hash (legacy `[[comparison:N]]`)
 *   - the resolved comparison is loading or dead
 *   - the token hash matches the first 8 chars of the live `content_hash`
 *
 * Returns `true` only when both halves are present and disagree.
 */
function computeHashDrift(
  resolved: ResolvedMention | "loading" | "dead",
  tokenHash: string | undefined,
): boolean {
  if (tokenHash === undefined) return false;
  if (resolved === "loading" || resolved === "dead") return false;
  if (resolved.type !== "comparison") return false;
  return resolved.data.content_hash.slice(0, 8).toLowerCase() !== tokenHash.toLowerCase();
}

export function MentionChip({ resolved, originalText, tokenHash }: ChipProps): JSX.Element {
  const [isHovered, setIsHovered] = useState(false);
  const setHoveredPeak  = useAppState((s) => s.setHoveredPeak);
  const setHoveredIndex = useAppState((s) => s.setHoveredIndex);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);
  const experimentQ = useExperiment(activeExperimentId ?? 0);
  const qUnits = experimentQ.data?.q_units ?? null;
  const latticeUnit   = latticeUnitFromQUnits(qUnits);
  const curvatureUnit = inverseSquareUnits(qUnits);

  const handleMouseEnter = useCallback(() => {
    setIsHovered(true);
    if (resolved === "loading" || resolved === "dead") return;
    if (resolved.type === "peak")  setHoveredPeak(resolved.data.id);
    if (resolved.type === "index") setHoveredIndex(resolved.data.id);
  }, [resolved, setHoveredPeak, setHoveredIndex]);

  const handleMouseLeave = useCallback(() => {
    setIsHovered(false);
    if (resolved === "loading" || resolved === "dead") return;
    if (resolved.type === "peak")  setHoveredPeak(undefined);
    if (resolved.type === "index") setHoveredIndex(undefined);
  }, [resolved, setHoveredPeak, setHoveredIndex]);

  const stateKey  = resolved === "loading" ? "loading" : resolved === "dead" ? "dead" : resolved.type;
  const baseStyle  = CHIP_STYLES[stateKey]  ?? CHIP_STYLES["dead"]!;
  const hoverStyle = CHIP_HOVER_STYLES[stateKey] ?? CHIP_HOVER_STYLES["dead"]!;

  const label = resolved === "loading" || resolved === "dead"
    ? (originalText ?? "…")
    : chipLabel(resolved);

  const hashDrift = computeHashDrift(resolved, tokenHash);

  // E2E selector hooks. Resolved entities expose data-mention-type/id so
  // tests can grab a chip without scraping rendered text. tokenHash + drift
  // are scoped to comparison chips per Phase 10 spec; for other types the
  // drift attribute is always "false".
  const dataAttrs: Record<string, string> = {
    "data-testid": "mention-chip",
    "data-mention-state": stateKey,
    "data-hash-drift": String(hashDrift),
  };
  if (resolved !== "loading" && resolved !== "dead") {
    dataAttrs["data-mention-type"] = resolved.type;
    dataAttrs["data-mention-id"]   = String(resolved.data.id);
  }

  const tooltip = isHovered && resolved !== "loading"
    ? (resolved === "dead"
        ? <span className="text-[#555555]">no longer exists</span>
        : <TooltipContent
            resolved={resolved}
            latticeUnit={latticeUnit}
            curvatureUnit={curvatureUnit}
            hashDrift={hashDrift}
          />)
    : null;

  return (
    <span
      {...dataAttrs}
      className={`relative inline border-b pb-px px-1 rounded-sm cursor-pointer whitespace-nowrap
                  text-sm transition-colors ${baseStyle} ${hoverStyle}`}
      onMouseEnter={handleMouseEnter}
      onMouseLeave={handleMouseLeave}
    >
      {label}
      {hashDrift && (
        <span
          data-testid="mention-chip-drift"
          className="ml-1 text-xs text-[#cca888] italic"
        >
          (changed)
        </span>
      )}
      {tooltip && (
        <span className="absolute bottom-[calc(100%+6px)] left-1/2 -translate-x-1/2 z-10
                         bg-[#252525] border border-[#3a3a3a] rounded-md px-2 py-1
                         text-xs text-[#999999] whitespace-nowrap shadow-lg pointer-events-none">
          {tooltip}
          <span className="absolute top-full left-1/2 -translate-x-1/2
                           border-4 border-transparent border-t-[#3a3a3a]" />
        </span>
      )}
    </span>
  );
}
