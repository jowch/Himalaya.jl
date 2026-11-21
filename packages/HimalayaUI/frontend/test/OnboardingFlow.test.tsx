import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { OnboardingFlow } from "../src/print/shell/OnboardingFlow";
import { useAppState } from "../src/state";
import * as api from "../src/api";

async function walkTutorial(): Promise<{ text: string; keys: string[] }> {
  const user = userEvent.setup();
  renderWithProviders(<OnboardingFlow />);
  await user.type(await screen.findByPlaceholderText("First name"), "Nora");
  await user.type(screen.getByPlaceholderText("handle"), "nora");
  await user.click(screen.getByRole("button", { name: "Continue" }));
  const text: string[] = [];
  const keys: string[] = [];
  for (;;) {
    const tutorial = await screen.findByTestId("onboarding-tutorial");
    text.push(tutorial.textContent ?? "");
    keys.push(...within(tutorial).queryAllByTestId("kbkey").map((k) => k.textContent ?? ""));
    const next = within(tutorial).queryByRole("button", { name: "Next" });
    if (!next) break;
    await user.click(next);
  }
  return { text: text.join("\n"), keys };
}

describe("OnboardingFlow tutorial", () => {
  beforeEach(() => {
    useAppState.getState().clearUsername();
    vi.spyOn(api, "listUsers").mockResolvedValue([]);
    vi.spyOn(api, "createUser").mockResolvedValue({ username: "nora", first_name: "Nora", last_name: null } as never);
  });

  it("describes the workspace, not the deleted corpus, grouping review or Focus pages", async () => {
    const { text } = await walkTutorial();
    expect(text).not.toMatch(/corpus|grouping|focus/i);
    expect(text).toMatch(/unreviewed sample/i);
  });

  it("teaches the triage keys as keycaps", async () => {
    const { keys } = await walkTutorial();
    expect(keys).toEqual(expect.arrayContaining(["X", "1–9", "0", "C", "Z", "?"]));
  });
});
