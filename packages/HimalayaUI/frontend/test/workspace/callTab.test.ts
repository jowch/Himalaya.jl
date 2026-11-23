import { describe, it, expect } from "vitest";
import { callChoices, cycleChoice, inCallChoices } from "../../src/workspace/callTab";

function list(): { root: HTMLElement; buttons: HTMLButtonElement[] } {
  const root = document.createElement("div");
  root.innerHTML = `<button>outside</button><div data-call-choices=""><button>1</button><button disabled>2</button><button><span>3</span></button><button>0</button></div>`;
  document.body.append(root);
  return { root, buttons: Array.from(root.querySelectorAll("button")) };
}

describe("callTab", () => {
  it("collects only the enabled buttons inside the phase choices", () => {
    const { root } = list();
    expect(callChoices(root).map((b) => b.textContent)).toEqual(["1", "3", "0"]);
    root.remove();
  });

  it("steps forward and back, wrapping at both ends", () => {
    const { root } = list();
    const [a, b, c] = callChoices(root) as [HTMLElement, HTMLElement, HTMLElement];
    expect(cycleChoice([a, b, c], a, false)).toBe(b);
    expect(cycleChoice([a, b, c], c, false)).toBe(a);
    expect(cycleChoice([a, b, c], a, true)).toBe(c);
    expect(cycleChoice([a, b, c], b.querySelector("span"), false)).toBe(c);
    root.remove();
  });

  it("enters at the first choice going forward and the last going back", () => {
    const { root, buttons } = list();
    const choices = callChoices(root);
    expect(cycleChoice(choices, buttons[0]!, false)).toBe(choices[0]);
    expect(cycleChoice(choices, document.body, true)).toBe(choices[2]);
    expect(cycleChoice([], document.body, false)).toBeUndefined();
    expect(inCallChoices(buttons[0]!)).toBe(false);
    expect(inCallChoices(choices[1]!.querySelector("span"))).toBe(true);
    root.remove();
  });
});
