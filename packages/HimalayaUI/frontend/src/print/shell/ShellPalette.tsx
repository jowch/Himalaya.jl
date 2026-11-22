import { useEffect, useMemo, useState } from "react";
import { useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../../state";
import { useCorpusSamples, useExperiments } from "../../queries";
import { useAppTheme } from "../../lib/theme";
import { sampleDisplayName } from "../../lib/sample/displayName";
import { Palette } from "../../workspace/ui/palette/Palette";
import { workspaceHref } from "../../lib/url/workspaceHref";
import { shellPalette } from "./shellPalette";

export function ShellPalette({ openHelp }: { openHelp: () => void }): JSX.Element | null {
  const open = useAppState((s) => s.paletteOpen);
  const close = useAppState((s) => s.closePalette);
  useEffect(() => close, [close]);
  return open ? <OpenShellPalette openHelp={openHelp} /> : null;
}

function OpenShellPalette({ openHelp }: { openHelp: () => void }): JSX.Element {
  const close = useAppState((s) => s.closePalette);
  const username = useAppState((s) => s.username);
  const clearUsername = useAppState((s) => s.clearUsername);
  const { pathname } = useLocation();
  const navigate = useNavigate();
  const theme = useAppTheme();
  const experimentsQ = useExperiments();
  const [searched, setSearched] = useState(false);
  const corpusQ = useCorpusSamples({ enabled: searched });
  const experiments = useMemo(
    () => (experimentsQ.data ?? []).map((e) => ({ id: e.id, name: e.name ?? `Experiment ${e.id}` })),
    [experimentsQ.data],
  );
  const samples = useMemo(() => {
    const names = new Map(experiments.map((e) => [e.id, e.name]));
    return (corpusQ.data ?? []).map((s) => ({
      id: s.id,
      name: sampleDisplayName(s),
      experiment: names.get(s.experiment_id) ?? `Experiment ${s.experiment_id}`,
      href: workspaceHref(s.experiment_id, s.id),
    }));
  }, [corpusQ.data, experiments]);
  const palette = shellPalette({
    pathname,
    experiments,
    samples,
    go: (path) => navigate(path),
    theme: { current: theme.theme, choose: theme.choose },
    user: { name: username, switch: clearUsername },
    help: openHelp,
  });
  const onType = (q: string): void => {
    if (q.trim()) setSearched(true);
  };
  return <Palette groups={palette.groups} pages={palette.pages} onClose={close} onType={onType} />;
}
