import { Dock, DockUpLink, DockStepper, Button, KbKey } from "../ui";
import type { Action } from "./types";
import { useInteraction } from "./registry";

function enabledOf(a: Action): boolean {
  return a.enabled ? a.enabled() : true;
}

export function InteractionDock(): JSX.Element | null {
  const cursor = useInteraction((s) => s.cursor);
  const actions = useInteraction((s) => s.actions);

  if (cursor === null && actions.length === 0) return null;

  const back = actions.find((a) => a.id === "back");
  const primary = actions.find((a) => a.dock === "primary");
  const buttons = actions.filter((a) => a.dock === true && a.id !== "back");
  const stepper = cursor ? cursor.stepperProps() : null;

  return (
    <Dock>
      {back && <DockUpLink label={back.label} onClick={() => back.run()} />}

      {stepper && (
        <DockStepper
          label={stepper.label}
          axis={stepper.axis}
          testIdBase={stepper.testIdBase}
          count={stepper.count}
          onPrev={stepper.onPrev}
          onNext={() => {
            stepper.onNext?.();
            cursor?.moveBy(1);
          }}
          prevDisabled={stepper.prevDisabled}
          nextDisabled={stepper.nextDisabled}
        />
      )}

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
