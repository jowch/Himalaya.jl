/**
 * Proves that `useListCursor` satisfies the shared cursor parity contract.
 *
 * If any invariant fails here it is a real finding — do NOT weaken the
 * assertion; investigate the cursor implementation instead.
 */
import { vi } from "vitest";
import { useListCursor } from "../../src/print/interaction/useListCursor";
import { DockStepper } from "../../src/print/ui";
import { runCursorContract } from "./cursorContract";

runCursorContract("useListCursor (base list)", () => {
  const onActivate = vi.fn();
  const IDS = [10, 20, 30];

  return {
    ui: (capture) => {
      function Wrapper(): JSX.Element {
        const cursor = useListCursor({
          ids: IDS,
          onActivate,
          stepperTestIdBase: "item",
        });

        // Call capture in the render body — not in a useEffect — so `m` is
        // always the cursor from the most recent render.  This is the fix for
        // the stale-cursor hazard described in cursorContract.tsx.
        capture({ cursor, ids: IDS, onActivate });

        return (
          <div role="grid" aria-multiselectable={true} data-interaction-scope="">
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>
                row {id}
              </div>
            ))}
            <DockStepper {...cursor.stepperProps()} />
          </div>
        );
      }

      return <Wrapper />;
    },
  };
});
