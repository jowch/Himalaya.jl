import { useEffect } from "react";
import { useLocation, useNavigate } from "react-router";
import { useAppState } from "../state";
import { useAppTheme } from "../lib/theme";
import { Palette } from "../ui/palette/Palette";
import { shellPalette } from "./shellPalette";
import { useCorpusPalette } from "../hooks/useCorpusPalette";

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
  const { experiments, samples, onType } = useCorpusPalette();
  const palette = shellPalette({
    pathname,
    experiments,
    samples,
    go: (path) => navigate(path),
    theme: { current: theme.theme, choose: theme.choose },
    user: { name: username, switch: clearUsername },
    help: openHelp,
  });
  return <Palette groups={palette.groups} pages={palette.pages} onClose={close} onType={onType} />;
}
