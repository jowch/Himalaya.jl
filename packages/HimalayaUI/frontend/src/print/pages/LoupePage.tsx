import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useLocation, useNavigate, useParams, useSearchParams } from "react-router-dom";
import { PageFrame } from "../components/PageFrame";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSamples,
  useExposures,
  useSetExposureStatus,
  useSelectExposure,
  useAddCorpusSampleTag,
  useRemoveCorpusSampleTag,
  useEditCorpusSampleTag,
  useCorpusSampleTags,
} from "../../queries";
import type { Tag } from "../ui";
import { EmptyState, Button, IconButton } from "../ui";
import { Dock } from "../ui/Dock";
import { resolveSampleOrder, sampleNeighbors } from "../../lib/sample/sampleOrder";
import { announce } from "../../lib/announce";
import { showToast } from "../../lib/toast";
import { useShortcuts } from "../shell/useShortcuts";
import { KbdLegend } from "../shell/KbdLegend";
import { isValidationError } from "../../lib/queue/errors";
import { BigFrame } from "../components/BigFrame";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { LoupeSidePanel } from "../components/LoupeSidePanel";
import { ManageTagsModal } from "../components/ManageTagsModal";
import { PlateHeader } from "../components/PlateHeader";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
} from "./loupeAdapters";
import type { SampleTagPair } from "../../api";

// Boneyard fixture — a real render with mock props so the headless capture CLI
// measures the greenfield loupe body. image_path:null → DetectorImage takes its
// placeholder branch (a plain rectangle at the frame's fixed aspect).
const FIXTURE_EXPOSURE = {
  id: 0, sample_id: 0, filename: "JC000-001.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null,
  image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
// EX-SPACING: the loupe body grid (detector hero + 286px side panel) is a
// skeleton↔real ALIGNMENT INVARIANT — the loading placeholder must use the same
// template as the loaded render, or the load→loaded transition shifts. One
// source (a literal Tailwind class string so the JIT scanner still sees it).
// LO-MOBILE: below the tablet breakpoint the fixed 286px panel + gap crushed the
// detector hero to a few px (390px viewport → ~76px hero). Stack to a single
// column under md; the two-column split holds at tablet width and up.
const LOUPE_BODY_GRID = "grid grid-cols-1 md:grid-cols-[minmax(0,1fr)_286px] gap-7";
const LOUPE_FIXTURE = (
  <div className={LOUPE_BODY_GRID}>
    <div className="min-w-0">
      <BigFrame src={null} caption="frame 1 of 1 · kept" accepted />
      <ThumbnailGallery
        exposures={[{ id: 0, src: null, frameNo: 1, kept: true }]}
        selectedId={0}
        size="lg"
        align="center"
        className="mt-3"
      />
    </div>
    <LoupeSidePanel
      meta={toMetaEntries(FIXTURE_EXPOSURE, [FIXTURE_EXPOSURE])}
      dropped={false}
      kept
      isRepresentative={false}
      tags={[]}
    />
  </div>
);

/**
 * Terminal save failure for the loupe's status verbs. The queue layer already
 * toasts 4xx validation failures (buildValidationMessage) — re-toasting here
 * would double up, so those are skipped. Infrastructure failures (5xx/network)
 * retry behind the Saving banner, but once the retry policy exhausts the
 * mutation settles as error, the banner disappears and the optimistic change
 * rolls back — without this toast that rollback would be silent.
 */
function notifySaveFailed(err: unknown): void {
  if (isValidationError(err)) return;
  showToast("Couldn't save the change. It was rolled back.", "error");
}

/**
 * LoupePage (greenfield) — the sample loupe at /sample/:sampleId/loupe.
 * URL-owned: the sample id is the route param, never Zustand `activeSampleId`.
 * Mounts body-only inside the app shell <Outlet> (T3.2).
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  const exposuresQ = useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);
  const isLoading = corpusQ.isLoading || exposuresQ.isLoading;

  // Active exposure — local state, defaulted by defaultExposureId, reset on
  // sample change so the next sample picks its own default.
  const [activeId, setActiveId] = useState<number | undefined>(undefined);
  useEffect(() => { setActiveId(undefined); }, [sampleId]);
  const computedDefault = defaultExposureId(exposures);

  // SA-F2: ?exposure=<id> preselects the frame the contact sheet was
  // double-clicked on. Read once for the INITIAL selection only — flipping
  // frames (arrows, thumb clicks) stays URL-silent; this page syncs no other
  // state to the URL, so the param is a permalink seed, not live state. An id
  // outside THIS sample's exposure list (unknown or foreign) is dropped
  // silently and the default wins — no error surface.
  const exposureParam = searchParams.get("exposure");
  const requestedId =
    exposureParam !== null && /^\d+$/.test(exposureParam)
      ? Number(exposureParam)
      : undefined;
  // LO-STALEACTIVEID: heal a selection that is unset OR no longer present. The
  // second case fires when another user drops the shown frame via SSE — without
  // it `activeExposure` resolves undefined and the page falsely claims "no
  // exposures" while frames still exist. The URL ?exposure= param seeds the
  // INITIAL selection only (activeId === undefined); a vanished-frame heal goes
  // straight to the computed default, never back to the stale permalink seed.
  const activeStillPresent = activeId !== undefined && exposures.some((e) => e.id === activeId);
  useEffect(() => {
    if (activeStillPresent) return;
    const seed =
      activeId === undefined && requestedId !== undefined && exposures.some((e) => e.id === requestedId)
        ? requestedId
        : computedDefault;
    if (seed !== undefined && seed !== activeId) setActiveId(seed);
  }, [activeStillPresent, activeId, computedDefault, exposures, requestedId]);
  // Resolve the rendered frame resiliently so the one paint before the heal
  // effect runs shows the surviving default, not the empty state. When the list
  // is non-empty `computedDefault` is always a real id, so `activeExposure` is
  // defined whenever a frame exists — the empty state is reserved for the
  // genuinely-frameless case (exposures.length === 0).
  const activeExposure =
    exposures.find((e) => e.id === activeId) ??
    exposures.find((e) => e.id === computedDefault);

  const frameIndex = activeExposure
    ? exposures.findIndex((e) => e.id === activeExposure.id)
    : -1;
  // LO-TERM: one word for a single detector image across the loupe — "frame"
  // (the BigFrame caption, toasts, and meta keys all say frame); never mix in
  // "exposure" on the same screen.
  const exposurePosition =
    frameIndex >= 0 ? `frame ${frameIndex + 1} of ${exposures.length}` : "—";

  const setStatus = useSetExposureStatus(hasValidId ? sampleId : 0);
  const setRepresentative = useSelectExposure(hasValidId ? sampleId : 0);
  const addTag = useAddCorpusSampleTag(hasValidId ? sampleId : 0);
  const removeTag = useRemoveCorpusSampleTag(hasValidId ? sampleId : 0);
  const editTag = useEditCorpusSampleTag(hasValidId ? sampleId : 0);

  // ManageTagsModal state
  const [manageOpen, setManageOpen] = useState(false);
  const manageTagsTriggerRef = useRef<HTMLButtonElement>(null);

  // Corpus-wide tag pairs for key/value suggestions
  const corpusTagsQ = useCorpusSampleTags();
  const corpusTags: SampleTagPair[] = corpusTagsQ.data ?? [];

  // Derive keyOptions: distinct keys ranked by descending total count
  const keyOptions = useMemo(() => {
    const counts = new Map<string, number>();
    for (const t of corpusTags) {
      counts.set(t.key, (counts.get(t.key) ?? 0) + (t.count ?? 1));
    }
    return Array.from(counts.entries())
      .sort((a, b) => b[1] - a[1])
      .map(([label, count]) => ({ label, count }));
  }, [corpusTags]);

  // Derive valueOptionsFor: values for a given key, ranked by count
  const valueOptionsFor = useCallback(
    (key: string) => {
      return corpusTags
        .filter((t) => t.key === key)
        .sort((a, b) => (b.count ?? 1) - (a.count ?? 1))
        .map((t) => ({
          label: t.value,
          ...(t.count !== undefined ? { count: t.count } : {}),
        }));
    },
    [corpusTags],
  );

  const handleDropToggle = useCallback(() => {
    if (!activeExposure) return;
    const dropping = activeExposure.status !== "rejected";
    // LO-LASTFRAME-NOGUARD: dropping the ONLY non-rejected frame leaves the
    // sample with no kept frame at all — a consequential edge the bland "Frame
    // dropped" copy hid. The action stays reversible (X toggles it back), so
    // name both the consequence and the reversal gesture rather than guarding
    // it outright (a sample legitimately may end up all-rejected).
    const lastKept =
      dropping && exposures.filter((e) => e.status !== "rejected").length === 1;
    setStatus.mutate(
      { exposureId: activeExposure.id, status: dropping ? "rejected" : null },
      {
        // Consequential single-frame status change → visible toast, fired on
        // CONFIRMATION (HTTP or own-op SSE), never optimistically: a premature
        // "Frame dropped" would lie if the save later failed and rolled back.
        onSuccess: () =>
          showToast(
            dropping
              ? lastKept
                ? "Last kept frame dropped. This sample now has no kept frame; press X to restore it."
                : "Frame dropped"
              : "Frame restored",
            lastKept ? "warning" : "success",
          ),
        onError: notifySaveFailed,
      },
    );
  }, [activeExposure, exposures, setStatus]);

  // The Keep verb (SA-SCREENED): K toggles accepted ↔ null. On a rejected
  // frame, K accepts directly — last verb wins, no trip through unscreened.
  const handleKeepToggle = useCallback(() => {
    if (!activeExposure) return;
    const keeping = activeExposure.status !== "accepted";
    setStatus.mutate(
      { exposureId: activeExposure.id, status: keeping ? "accepted" : null },
      {
        onSuccess: () =>
          showToast(keeping ? "Frame kept" : "Frame restored", "success"),
        onError: notifySaveFailed,
      },
    );
  }, [activeExposure, setStatus]);

  // Restore: set the active frame's status back to null (unscreened). Mirrors
  // the drop/keep toggle path — same setStatus mutation, status: null payload.
  const handleRestore = useCallback(() => {
    if (!activeExposure) return;
    setStatus.mutate(
      { exposureId: activeExposure.id, status: null },
      {
        onSuccess: () => showToast("Frame restored", "success"),
        onError: notifySaveFailed,
      },
    );
  }, [activeExposure, setStatus]);

  const handleSetRepresentative = useCallback(() => {
    if (!activeExposure) return;
    if (activeExposure.selected) {
      // Already the representative — no mutation, no false success toast.
      // Quiet SR-only acknowledgement keeps the R key honest.
      announce("Already the representative frame");
      return;
    }
    setRepresentative.mutate(activeExposure.id, {
      onSuccess: () => showToast("Set as the representative frame", "success"),
      onError: notifySaveFailed,
    });
  }, [activeExposure, setRepresentative]);

  const handleAddTag = useCallback((t: Tag) => {
    addTag.mutate({ key: t.key, value: t.value ?? "" });
  }, [addTag]);

  const handleRemoveTag = useCallback((tagId: number) => {
    removeTag.mutate(tagId);
  }, [removeTag]);

  const flip = useCallback((delta: number) => {
    if (activeId === undefined || exposures.length === 0) return;
    const idx = exposures.findIndex((e) => e.id === activeId);
    if (idx < 0) return;
    const next = Math.min(Math.max(idx + delta, 0), exposures.length - 1);
    setActiveId(exposures[next]!.id);
    // High-frequency, immediately-visible navigation → SR-only (not a toast).
    announce(`Frame ${next + 1} of ${exposures.length}`);
  }, [activeId, exposures]);

  // Thumbnail-click selection routes through the SAME announcement as the
  // keyboard flip so mouse and keyboard navigation are consistent for SR users.
  const selectFrame = useCallback((id: number) => {
    setActiveId(id);
    const idx = exposures.findIndex((e) => e.id === id);
    if (idx >= 0) announce(`Frame ${idx + 1} of ${exposures.length}`);
  }, [exposures]);

  const goBack = useCallback(() => {
    const experimentParam = searchParams.get("experiment");
    // LO-FOCUSRET (WCAG 2.4.3): carry the originating sample id back so the
    // sheet restores focus to that row instead of dropping it to <body>.
    const opts = hasValidId ? { state: { focusSampleId: sampleId } } : undefined;
    navigate(experimentParam ? `/experiments?experiment=${experimentParam}` : "/experiments", opts);
  }, [navigate, searchParams, hasValidId, sampleId]);

  // ── LO-NEXT: prev/next-SAMPLE navigation ──────────────────────────────────
  // Culling N samples should not cost N sheet round-trips. The sheet hands its
  // visible, sorted+filtered sample order through router state when it opens
  // the loupe; we walk THAT list so prev/next match exactly what the user saw.
  // A direct URL (permalink, no state) falls back to the experiment-scoped corpus
  // order, which also matches the sheet's default (ingest order).
  const location = useLocation();
  const experimentParam2 = searchParams.get("experiment");
  const experimentId =
    experimentParam2 !== null && /^\d+$/.test(experimentParam2)
      ? Number(experimentParam2)
      : undefined;
  // Shared with the app shell sample stepper (resolveSampleOrder), so the
  // topbar's ‹ › and the loupe's own `[`/`]` keyboard nav always agree.
  const orderedSampleIds = useMemo(
    () =>
      resolveSampleOrder(
        corpusQ.data ?? [],
        experimentId,
        sampleId,
        (location.state as { sampleOrder?: number[] } | null)?.sampleOrder,
      ),
    [location.state, corpusQ.data, experimentId, sampleId],
  );
  const { prevId: prevSampleId, nextId: nextSampleId } = sampleNeighbors(
    orderedSampleIds,
    sampleId,
  );
  const gotoSample = useCallback(
    (id: number): void => {
      const params = new URLSearchParams();
      if (experimentId !== undefined) params.set("experiment", String(experimentId));
      const qs = params.toString();
      // Carry the SAME order forward (so the next step still has the list) and
      // drop ?exposure= so each sample opens at its own default frame.
      navigate(`/sample/${id}/loupe${qs ? `?${qs}` : ""}`, {
        state: { sampleOrder: orderedSampleIds },
      });
    },
    [navigate, experimentId, orderedSampleIds],
  );

  // Loupe keyboard, via the shared shortcut library. The vocabulary nests:
  // ← → step the FRAME, [ ] step the SAMPLE; x/k/r screen the active frame;
  // Esc backs out to the sheet. Modifier chords (⌘R reload, ⌘X cut) never match
  // the registry's bare combos, and typing/modal suppression is handled inside
  // useShortcuts — so the per-key guards the old ad-hoc handler spelled out are
  // now owned uniformly by the library.
  useShortcuts({
    // ←/→ = frame (detail) axis — renamed from prevExposure/nextExposure in T2.5.
    prevDetail: () => flip(-1),
    nextDetail: () => flip(1),
    drop: () => handleDropToggle(),
    keep: () => handleKeepToggle(),
    representative: () => handleSetRepresentative(),
    restore: () => handleRestore(),
    prevSample: () => {
      if (prevSampleId !== undefined) gotoSample(prevSampleId);
    },
    nextSample: () => {
      if (nextSampleId !== undefined) gotoSample(nextSampleId);
    },
    dismiss: () => goBack(),
  });

  if (!corpusQ.isLoading && !sample) {
    return (
      <PageFrame width="loupe" className="px-8 py-7">
        <div data-testid="loupe-page">
          <div data-testid="loupe-not-found">
            <EmptyState
              as="h1"
              title="Sample not found"
              body="Nothing in the corpus matches this address."
              action={
                <Button variant="outline" onClick={goBack}>
                  Back to the sheet
                </Button>
              }
            />
          </div>
        </div>
      </PageFrame>
    );
  }

  const isDropped = activeExposure?.status === "rejected";
  const isKept = activeExposure?.status === "accepted";
  // Honest tri-state caption: a null status is unscreened, never "kept".
  const verdictWord = isDropped ? "dropped" : isKept ? "kept" : "unscreened";
  // Sample-level truth, not frame-level: the backend's Index-stage resolution
  // never consults status, so a dropped representative still carries forward.
  const representativeDropped = exposures.some(
    (e) => e.selected && e.status === "rejected",
  );

  return (
    <PageFrame width="loupe" className="px-8 py-7">
      <div data-testid="loupe-page">
        {/* The inter-sample stepper lives in the Dock (§3.3); the loupe's own
            `[`/`]` keyboard nav still steps via gotoSample, sharing
            resolveSampleOrder. */}
        <div className="mb-3.5 flex items-center gap-3">
          <button data-testid="loupe-back" onClick={goBack} className="text-sm font-semibold text-print-accent hover:underline">
            ← Back to the sheet
          </button>
        </div>
        <PlateHeader
          as="h1"
          title={sample?.name ?? "—"}
          subtitle={`${sample?.name ?? "—"} · ${exposurePosition}`}
          className="mb-5"
        />
        <Skeleton name="loupe" className="block" loading={isLoading} stagger={50} transition={200}
          fixture={LOUPE_FIXTURE}
          fallback={<div data-testid="loupe-skeleton" className="p-8 text-sm italic text-ink-soft">Loading sample…</div>}>
          <div className={LOUPE_BODY_GRID}>
            {sample && activeExposure ? (
              <>
                {/* LO-STRIP-OVERHANG: cap the whole detector stack (frame +
                    filmstrip + legend) at the frame's max width and centre it, so
                    the filmstrip aligns to the frame edges instead of overhanging
                    it by ~a hundred px each side. */}
                <div className="min-w-0 w-full max-w-[640px] mx-auto">
                  <BigFrame
                    src={buildExposureImageUrl(activeExposure)}
                    caption={`frame ${frameIndex + 1} of ${exposures.length} · ${verdictWord}`}
                    rejected={isDropped}
                    accepted={!!isKept}
                  />
                  <ThumbnailGallery
                    exposures={toGalleryExposures(exposures)}
                    {...(activeExposure ? { selectedId: activeExposure.id } : {})}
                    onSelect={selectFrame}
                    size="lg"
                    align="center"
                    className="mt-3"
                  />
                  {/* LO-KBDLEGEND: a registry-driven key legend so the loupe's
                      shortcuts are discoverable on the surface itself. Curated to
                      the screen verbs ALONE — the frame arrows are conventional
                      (and their registry label says "exposure", which LO-TERM
                      keeps off the loupe), and the sample steps ([ ]) are
                      advertised by the Dock's stepper buttons (LO-STEPDEDUP),
                      so repeating them here is redundant. */}
                  <KbdLegend
                    ids={["drop", "keep", "representative"]}
                    testId="loupe-kbd-legend"
                    className="mt-4"
                  />
                </div>
                <LoupeSidePanel
                  meta={toMetaEntries(activeExposure, exposures)}
                  dropped={!!isDropped}
                  kept={!!isKept}
                  isRepresentative={activeExposure.selected}
                  representativeDropped={representativeDropped}
                  tags={toLoupeTags(sample.tags)}
                  onToggleDrop={handleDropToggle}
                  onToggleKeep={handleKeepToggle}
                  onSetRepresentative={handleSetRepresentative}
                  onAddTag={handleAddTag}
                  onRemoveTag={handleRemoveTag}
                  onManageTags={() => setManageOpen(true)}
                  manageTagsTriggerRef={manageTagsTriggerRef}
                />
                <ManageTagsModal
                  open={manageOpen}
                  sampleName={sample.name ?? ""}
                  tags={toLoupeTags(sample.tags)}
                  keyOptions={keyOptions}
                  valueOptionsFor={valueOptionsFor}
                  onEdit={(tagId, key, value) => editTag.mutate({ tagId, key, value })}
                  onAdd={(key, value) => addTag.mutate({ key, value: value })}
                  onRemove={(tagId) => removeTag.mutate(tagId)}
                  onClose={() => setManageOpen(false)}
                  triggerRef={manageTagsTriggerRef}
                />
              </>
            ) : (
              <div className="col-span-2 p-8 text-sm text-ink-soft">This sample has no exposures.</div>
            )}
          </div>
        </Skeleton>

        {/* ── Contextual bottom dock (Loupe grammar §3.3) ──────────────────────────
            ‹ Corpus · Sample↑↓ · Frame‹› · Drop · Keep · Set representative · Restore · Focus
            Each verb calls the SAME callback the keyboard shortcut uses. */}
        <Dock>
          {/* Up-link back to corpus */}
          <button
            onClick={goBack}
            className="text-meta font-semibold text-print-accent hover:underline mr-1"
            data-testid="dock-up-link"
          >
            ‹ Corpus
          </button>

          <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

          {/* Sample↑↓ stepper */}
          <IconButton
            label="Previous sample"
            tone="ghost"
            disabled={prevSampleId === undefined}
            onClick={() => prevSampleId !== undefined && gotoSample(prevSampleId)}
            data-testid="dock-prev-sample"
          >
            ↑
          </IconButton>
          <IconButton
            label="Next sample"
            tone="ghost"
            disabled={nextSampleId === undefined}
            onClick={() => nextSampleId !== undefined && gotoSample(nextSampleId)}
            data-testid="dock-next-sample"
          >
            ↓
          </IconButton>

          <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

          {/* Frame‹› stepper */}
          <IconButton
            label="Previous frame"
            tone="ghost"
            disabled={activeExposure === undefined || exposures.indexOf(activeExposure) <= 0}
            onClick={() => flip(-1)}
            data-testid="dock-prev-frame"
          >
            ‹
          </IconButton>
          <IconButton
            label="Next frame"
            tone="ghost"
            disabled={activeExposure === undefined || exposures.indexOf(activeExposure) >= exposures.length - 1}
            onClick={() => flip(1)}
            data-testid="dock-next-frame"
          >
            ›
          </IconButton>

          <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

          {/* Cull verbs — coloured outlines (spec §3.3) */}
          <Button
            variant="outlineAccent"
            onClick={handleDropToggle}
            data-testid="dock-drop"
          >
            Drop
          </Button>
          <Button
            variant="outlineSuccess"
            onClick={handleKeepToggle}
            data-testid="dock-keep"
          >
            Keep
          </Button>

          <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

          {/* Set representative (Loupe-only — NOT on Corpus) */}
          <Button
            variant="ghost"
            onClick={handleSetRepresentative}
            data-testid="dock-set-representative"
          >
            Set representative
          </Button>

          {/* Restore */}
          <Button
            variant="ghost"
            onClick={handleRestore}
            data-testid="dock-restore"
          >
            Restore
          </Button>

          <span className="w-px self-stretch bg-hair mx-1" aria-hidden />

          {/* Focus destination */}
          <Button
            variant="accent"
            onClick={() => {
              if (hasValidId) navigate(`/sample/${sampleId}`);
            }}
            data-testid="dock-focus"
          >
            Focus
          </Button>
        </Dock>
      </div>
    </PageFrame>
  );
}
