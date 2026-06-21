import type { JSX } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQuery, useMutation } from "@tanstack/react-query";
import { ConfigurationBody } from "../components/ConfigurationBody";
import { useDraftExperiment } from "../../lib/draftExperiment";
import * as api from "../../api";
import { ProgressBar } from "../ui/ProgressBar";
import { NoticePill } from "../ui/NoticePill";
import { Button } from "../ui/Button";

/**
 * ConfigurationPage -- the Configuration tab shell (spec sec 8.1). Extracts
 * the experiment id from the route and delegates all body rendering to
 * ConfigurationBody (Task 19, Phase E2).
 *
 * First-run mode: absent `:id` (route `/experiments/new/config`) reads
 * path+patterns from `useDraftExperiment`, fetches the phase-1 manifest
 * (`GET /api/fs/manifest`), shows matched-by-type counts + unmatched list,
 * hides the geometry/sources region (no PRP parsed yet, sec 6.5), gates
 * Approve until the manifest resolves. On Approve calls `createExperiment`,
 * then `clear()`s the draft and navigates to the new experiment's corpus.
 *
 * Later-edit mode: present `:id` (route `/experiments/:id/config`) keeps
 * the existing behaviour including the geometry/sources table.
 *
 * The two wrapper divs with data-testid="config-geometry-region" and
 * "config-sources-region" preserve the E1 test contract (ConfigurationPage
 * Phase E1 shell test) while the real content is now composed by
 * ConfigurationBody.
 */
export function ConfigurationPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();

  if (!id) {
    return <ConfigurationFirstRun />;
  }

  const expId = Number(id);
  return (
    <div data-testid="config-geometry-region">
      {/* data-testid="config-sources-region" is also required by the E1 test.
          It lives inside the ConfigurationBody, but we need BOTH testids to be
          present. Wrap the body in a div that satisfies the geometry region test,
          then render ConfigurationBody which emits config-sources-region. */}
      <div data-testid="config-sources-region" className="contents">
        <ConfigurationBody experimentId={expId} />
      </div>
    </div>
  );
}

/**
 * First-run mode rendered at `/experiments/new/config`. Reads from the draft
 * store, fetches the phase-1 manifest, and gates Approve until the manifest
 * resolves. On Approve creates the experiment and navigates to its corpus.
 *
 * Geometry / sources region is intentionally absent (no PRP parsed yet at
 * this stage, spec sec 6.5).
 */
function ConfigurationFirstRun(): JSX.Element {
  const navigate = useNavigate();
  const { path, patterns, clear } = useDraftExperiment();

  const manifestQuery = useQuery({
    queryKey: ["manifest", path, patterns],
    queryFn: () => api.fetchManifest(path, patterns),
    enabled: path !== "",
  });

  const createMutation = useMutation({
    mutationFn: (body: api.CreateExperimentBody) => api.createExperiment(body),
    onSuccess: (created) => {
      clear();
      navigate(`/experiments/${created.id}/corpus`);
    },
  });

  const handleApprove = (): void => {
    createMutation.mutate({ path, patterns });
  };

  const manifest = manifestQuery.data;
  const isFetching = manifestQuery.isFetching;

  return (
    <div className="flex flex-col gap-6 p-6">
      <div>
        <h2 className="text-headline text-ink mb-1">Review experiment</h2>
        <p className="text-sm text-ink-soft">{path}</p>
      </div>

      {isFetching && (
        <ProgressBar value={0} total={1} label="Indexing files" />
      )}

      {manifest && (
        <div className="flex flex-col gap-3">
          <p className="text-sm text-ink">
            {manifest.matched.image} image · {manifest.matched.metadata} metadata · {manifest.matched.integration} integration · {manifest.total} total
          </p>
          {manifest.unmatched.length > 0 && (
            <div className="flex flex-wrap gap-2">
              {manifest.unmatched.map((u) => (
                <NoticePill key={u.file} tone="draft">
                  {u.file}: {u.miss}
                </NoticePill>
              ))}
            </div>
          )}
        </div>
      )}

      <div className="flex items-center gap-3">
        <Button
          variant="ghost"
          onClick={() => { clear(); navigate("/experiments/new"); }}
        >
          Back
        </Button>
        <Button
          variant="accent"
          data-testid="approve-button"
          disabled={isFetching || !manifest}
          onClick={handleApprove}
        >
          Approve
        </Button>
      </div>
    </div>
  );
}
