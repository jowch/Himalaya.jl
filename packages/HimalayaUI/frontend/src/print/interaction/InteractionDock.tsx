import { Dock, DockUpLink, DockStepper, Button, KbKey } from "../ui";
import type { Action } from "./types";
import { useInteraction } from "./registry";

function enabledOf(a: Action): boolean {
  return a.enabled ? a.enabled() : true;
}

export function InteractionDock(): JSX.Element | null {
  const cursor = useInteraction((s) => s.cursor);
  const actions = useInteraction((s) => s.actions);
  const extraSteppers = useInteraction((s) => s.extraSteppers);
  const dockExtra = useInteraction((s) => s.dockExtra);

  if (cursor === null && actions.length === 0 && extraSteppers.length === 0 && !dockExtra) return null;

  const back = actions.find((a) => a.id === "back");
  const primary = actions.find((a) => a.dock === "primary");
  // Mode-gated actions (mode !== undefined) are HIDDEN entirely when enabled()
  // is false — a verb that only lives in "selection" mode shouldn't sit greyed in
  // browse mode. Non-mode actions still render (and grey) when disabled. NOTE: if
  // a future action is disabled for a NON-mode reason (e.g. transiently loading)
  // and declares a mode, this would hide rather than grey it — revisit then.
  const buttons = actions.filter((a) => a.dock === true && a.id !== "back" && (a.mode === undefined || enabledOf(a)));
  const stepper = cursor ? cursor.stepperProps() : null;

  return (
    <Dock>
      {back && <DockUpLink label={back.label} onClick={() => back.run()} />}

      {/* Extra steppers (e.g. sample axis) render BEFORE the cursor stepper. */}
      {extraSteppers.map((s) => (
        <DockStepper
          key={s.testIdBase}
          label={s.label}
          axis={s.axis}
          testIdBase={s.testIdBase}
          count={s.count}
          onPrev={s.onPrev}
          onNext={s.onNext}
          prevDisabled={s.prevDisabled}
          nextDisabled={s.nextDisabled}
        />
      ))}

      {stepper && (
        <DockStepper
          label={stepper.label}
          axis={stepper.axis}
          testIdBase={stepper.testIdBase}
          count={stepper.count}
          onPrev={stepper.onPrev}
          onNext={stepper.onNext}
          prevDisabled={stepper.prevDisabled}
          nextDisabled={stepper.nextDisabled}
        />
      )}

      {/* dockExtra: non-interactive slot between stepper and action buttons.
          Renders page-local content (e.g. cursor identity readout) without
          needing an Action wrapper. InteractionDock owns spacing. */}
      {dockExtra}

      <div className="flex-1" />

      {buttons.map((a) => (
        <Button
          key={a.id}
          variant="outline"
          disabled={!enabledOf(a)}
          onClick={(e) => a.run(e)}
          data-testid={`dock-action-${a.id}`}
        >
          {a.label}
          {a.keys?.[0] && <KbKey className="ml-1.5">{a.keys[0]}</KbKey>}
        </Button>
      ))}

      {primary && (
        <Button
          variant="accent"
          disabled={!enabledOf(primary)}
          onClick={(e) => primary.run(e)}
          data-testid="dock-primary"
        >
          {primary.label}
          <KbKey variant="frost" className="ml-1.5">↵</KbKey>
        </Button>
      )}
    </Dock>
  );
}
