/**
 * SeriesRecipeEditor — the I3.5b recipe-edit controls injected into the
 * SeriesBuilderRail's edit slot. Drives the Zustand `seriesDraft` (optimistic)
 * and flushes `useSaveSeries` (recipe save) + `useCommitSeriesPlate` (spinner
 * commit). Renders only when a draft for this series is active.
 *
 *   - recipe list: each row removable + reorderable (▲/▼); newly-added rows
 *     carry a negative placeholder id (visually unmarked — a render detail).
 *   - add-sample: a corpus-sample dropdown → `addSeriesSample` (placeholder id).
 *   - ordering-variable + order-rule editors.
 *   - Save recipe (optimistic via draft → PATCH), Commit plate (spinner → POST
 *     commit, carries expected_content_hash), Cancel (discard draft).
 *
 * A double-click ref guards each flush (each mutate mints a fresh client_op_id,
 * so the idempotency layer can't dedupe a racing pair — guard synchronously).
 */
import { useCallback, useEffect, useMemo, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { useCorpusSamples, useSaveSeries, useCommitSeriesPlate } from "../queries";
import { buildSeriesSaveBody } from "../lib/series/buildSeriesSaveBody";
import { buildSeriesCommitBody } from "../lib/series/buildSeriesCommitBody";
import type { SeriesMember, OrderRule } from "../api";

const ORDER_RULES: readonly OrderRule[] = ["ascending", "descending", "manual"];

export function SeriesRecipeEditor({
  seriesId, members,
}: {
  seriesId: number;
  /** The loaded series' current plate — the source for the commit body. */
  members: SeriesMember[];
}): JSX.Element | null {
  const draft = useAppState((s) => s.seriesDraft);
  const addSample = useAppState((s) => s.addSeriesSample);
  const removeSample = useAppState((s) => s.removeSeriesSample);
  const reorderSample = useAppState((s) => s.reorderSeriesSample);
  const setTitle = useAppState((s) => s.setSeriesDraftTitle);
  const setDescription = useAppState((s) => s.setSeriesDraftDescription);
  const setOrderingVariable = useAppState((s) => s.setSeriesOrderingVariable);
  const setOrderRule = useAppState((s) => s.setSeriesOrderRule);
  const discardSeriesDraft = useAppState((s) => s.discardSeriesDraft);

  const navigate = useNavigate();
  const save = useSaveSeries();
  const commit = useCommitSeriesPlate();

  const corpus = useCorpusSamples();
  const sampleNameById = useMemo(() => {
    const m = new Map<number, string>();
    for (const s of corpus.data ?? []) {
      m.set(s.id, s.display_name || s.name || `Sample #${s.id}`);
    }
    return m;
  }, [corpus.data]);

  const saveInFlight = useRef(false);
  const commitInFlight = useRef(false);

  const handleSaveRecipe = useCallback(() => {
    if (draft === null) return;
    if (save.isPending || saveInFlight.current) return;
    saveInFlight.current = true;
    // Recipe save is a PATCH (id present); NO expected_content_hash.
    save.mutate({ id: seriesId, ...buildSeriesSaveBody(draft) });
  }, [draft, save, seriesId]);

  const handleCommit = useCallback(() => {
    if (draft === null) return;
    if (commit.isPending || commitInFlight.current) return;
    commitInFlight.current = true;
    commit.mutate({ id: seriesId, ...buildSeriesCommitBody(draft, members) });
  }, [draft, commit, seriesId, members]);

  const handleCancel = useCallback(() => {
    discardSeriesDraft();
  }, [discardSeriesDraft]);

  // Release the save guard once the recipe save settles (success or error);
  // the draft stays so the user can keep editing / then commit.
  useEffect(() => {
    if ((save.isSuccess || save.error) && saveInFlight.current) {
      saveInFlight.current = false;
    }
  }, [save.isSuccess, save.error]);

  // Commit success → discard the draft + return to the read surface (the
  // committed plate). A 409 surfaces via the conflict bridge →
  // SeriesCommitConflictModal; release the guard so the modal can retry.
  useEffect(() => {
    if (commit.isSuccess && commitInFlight.current) {
      commitInFlight.current = false;
      discardSeriesDraft();
      navigate(`/series/${seriesId}`);
    }
  }, [commit.isSuccess, discardSeriesDraft, navigate, seriesId]);
  useEffect(() => {
    if (commit.error && commitInFlight.current) commitInFlight.current = false;
  }, [commit.error]);

  if (draft === null) return null;

  const addableSamples = (corpus.data ?? []).filter(
    (s) => !draft.recipe.some((r) => r.sample_id === s.id),
  );

  return (
    <div data-testid="series-recipe-editor" className="flex flex-col gap-3">
      <label className="flex flex-col gap-1">
        <span className="text-xs font-semibold text-ink-faint">Title</span>
        <input
          data-testid="recipe-title"
          type="text"
          value={draft.title}
          onChange={(e) => setTitle(e.target.value)}
          className="rounded border border-hair bg-paper px-2 py-1 text-sm text-ink"
        />
      </label>

      <label className="flex flex-col gap-1">
        <span className="text-xs font-semibold text-ink-faint">Description</span>
        <textarea
          data-testid="recipe-description"
          value={draft.description}
          onChange={(e) => setDescription(e.target.value)}
          rows={2}
          className="rounded border border-hair bg-paper px-2 py-1 text-sm text-ink"
        />
      </label>

      <div className="text-xs font-semibold uppercase tracking-wide text-ink">Recipe</div>

      <ol data-testid="recipe-list" className="flex flex-col gap-1">
        {draft.recipe.map((row, i) => (
          <li
            key={row.id}
            data-testid="recipe-row"
            data-sample-id={String(row.sample_id)}
            className="flex items-center gap-1.5 text-sm text-ink"
          >
            <span className="flex-1 truncate">
              {sampleNameById.get(row.sample_id) ?? `Sample #${row.sample_id}`}
            </span>
            <button
              type="button"
              data-testid="recipe-row-up"
              disabled={i === 0}
              onClick={() => reorderSample(i, i - 1)}
              title="Move up"
              className="rounded px-1 text-ink-faint hover:text-ink disabled:opacity-30"
            >
              ▲
            </button>
            <button
              type="button"
              data-testid="recipe-row-down"
              disabled={i === draft.recipe.length - 1}
              onClick={() => reorderSample(i, i + 1)}
              title="Move down"
              className="rounded px-1 text-ink-faint hover:text-ink disabled:opacity-30"
            >
              ▼
            </button>
            <button
              type="button"
              data-testid="recipe-row-remove"
              onClick={() => removeSample(row.id)}
              title="Remove from recipe"
              className="rounded px-1 text-ink-faint hover:text-ink"
            >
              ✕
            </button>
          </li>
        ))}
        {draft.recipe.length === 0 && (
          <li className="text-sm text-ink-faint">No samples in this recipe yet.</li>
        )}
      </ol>

      <label className="flex flex-col gap-1">
        <span className="text-xs font-semibold text-ink-faint">Add sample</span>
        <select
          data-testid="recipe-add-sample"
          value=""
          onChange={(e) => {
            const id = Number(e.target.value);
            if (Number.isFinite(id) && id > 0) addSample(id);
          }}
          className="rounded border border-hair bg-paper px-2 py-1 text-sm text-ink"
        >
          <option value="">Select a sample…</option>
          {addableSamples.map((s) => (
            <option key={s.id} value={s.id}>
              {s.display_name || s.name || `Sample #${s.id}`}
            </option>
          ))}
        </select>
      </label>

      <label className="flex flex-col gap-1">
        <span className="text-xs font-semibold text-ink-faint">Ordering variable</span>
        <input
          data-testid="recipe-ordering-variable"
          type="text"
          value={draft.orderingVariable ?? ""}
          onChange={(e) => setOrderingVariable(e.target.value === "" ? null : e.target.value)}
          className="rounded border border-hair bg-paper px-2 py-1 text-sm text-ink"
        />
      </label>

      <label className="flex flex-col gap-1">
        <span className="text-xs font-semibold text-ink-faint">Order rule</span>
        <select
          data-testid="recipe-order-rule"
          value={draft.orderRule}
          onChange={(e) => setOrderRule(e.target.value as OrderRule)}
          className="rounded border border-hair bg-paper px-2 py-1 text-sm text-ink"
        >
          {ORDER_RULES.map((r) => (
            <option key={r} value={r}>{r}</option>
          ))}
        </select>
      </label>

      <div className="flex items-center gap-2 pt-1">
        <button
          type="button"
          data-testid="recipe-cancel"
          onClick={handleCancel}
          className="rounded border border-hair px-2 py-1 text-sm text-ink"
        >
          Cancel
        </button>
        <button
          type="button"
          data-testid="recipe-save"
          onClick={handleSaveRecipe}
          disabled={save.isPending}
          className="rounded border border-hair px-2 py-1 text-sm text-ink disabled:opacity-50"
        >
          {save.isPending ? "Saving…" : "Save recipe"}
        </button>
        <button
          type="button"
          data-testid="recipe-commit"
          onClick={handleCommit}
          // Disable while a recipe save is in flight: the commit body is built
          // from the loaded `members` (the server-resolved plate), NOT from the
          // edited `draft.recipe`, and the backend commits verbatim with no
          // re-resolution (routes_series.jl:247 → series.jl). Committing during
          // the save round-trip would freeze the STALE plate. Blocking the
          // window closes that foot-gun (round-1 nit).
          disabled={commit.isPending || save.isPending}
          title={save.isPending ? "Finish saving the recipe before committing" : undefined}
          className="rounded border border-print-accent bg-print-accent px-2 py-1 text-sm text-paper disabled:opacity-50"
        >
          {commit.isPending ? "Committing…" : "Commit plate"}
        </button>
      </div>
    </div>
  );
}
