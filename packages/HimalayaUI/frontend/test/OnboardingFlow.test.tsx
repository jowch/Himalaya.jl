import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { OnboardingFlow } from "../src/shell/OnboardingFlow";
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
    keys.push(...[...tutorial.querySelectorAll("kbd")].map((k) => k.textContent ?? ""));
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
    vi.spyOn(api, "createUser").mockImplementation(async (handle) => ({ username: handle, first_name: null, last_name: null }) as never);
    useAppState.getState().setTutorialSeen(false);
    useAppState.getState().setTutorialPendingFor(undefined);
  });

  afterEach(() => vi.restoreAllMocks());

  it("describes the workspace, not the deleted corpus, grouping review or Focus pages", async () => {
    const { text } = await walkTutorial();
    expect(text).not.toMatch(/corpus|grouping|focus/i);
    expect(text).toMatch(/unreviewed sample/i);
  });

  it("teaches the triage keys as keycaps", async () => {
    const { keys } = await walkTutorial();
    expect(keys).toEqual(expect.arrayContaining(["X", "1–9", "0", "C", "Z", "?"]));
  });

  it("after a reload mid-tutorial, preselects the user it was started for and shows it again", async () => {
    const user = userEvent.setup();
    const first = renderWithProviders(<OnboardingFlow />);
    await user.type(await screen.findByPlaceholderText("First name"), "Nora");
    await user.type(screen.getByPlaceholderText("handle"), "nora");
    await user.click(screen.getByRole("button", { name: "Continue" }));
    await screen.findByTestId("onboarding-tutorial");
    first.unmount();

    vi.mocked(api.listUsers).mockResolvedValue([
      { username: "ada", first_name: "Ada", last_name: null } as never,
      { username: "nora", first_name: "Nora", last_name: null } as never,
    ]);
    renderWithProviders(<OnboardingFlow />);
    await screen.findByRole("option", { name: /nora/ });
    expect(screen.getByRole("combobox")).toHaveValue("nora");
    await user.click(screen.getByRole("button", { name: "Continue" }));
    expect(await screen.findByTestId("onboarding-tutorial")).toBeInTheDocument();
  });

  it("does not show the tutorial to an existing user signing in on a fresh browser", async () => {
    vi.mocked(api.listUsers).mockResolvedValue([{ username: "ada", first_name: "Ada", last_name: null } as never]);
    const user = userEvent.setup();
    renderWithProviders(<OnboardingFlow />);
    await screen.findByRole("option", { name: /ada/ });
    await user.click(screen.getByRole("button", { name: "Continue" }));
    expect(useAppState.getState().username).toBe("ada");
    expect(screen.queryByTestId("onboarding-tutorial")).not.toBeInTheDocument();
  });

  it("keeps the arrow keys working after clicking Back", async () => {
    const user = userEvent.setup();
    renderWithProviders(<OnboardingFlow />);
    await user.type(await screen.findByPlaceholderText("First name"), "Nora");
    await user.type(screen.getByPlaceholderText("handle"), "nora");
    await user.click(screen.getByRole("button", { name: "Continue" }));
    await screen.findByTestId("onboarding-tutorial");
    await user.click(screen.getByRole("button", { name: "Next" }));
    await user.click(screen.getByRole("button", { name: /Back/ }));
    expect(screen.getByText(/1 of 3/)).toBeInTheDocument();
    await user.keyboard("{ArrowRight}");
    expect(screen.getByText(/2 of 3/)).toBeInTheDocument();
  });

  it("steps slides with the arrow keys without clicking first", async () => {
    const user = userEvent.setup();
    renderWithProviders(<OnboardingFlow />);
    await user.type(await screen.findByPlaceholderText("First name"), "Nora");
    await user.type(screen.getByPlaceholderText("handle"), "nora");
    await user.click(screen.getByRole("button", { name: "Continue" }));
    await screen.findByTestId("onboarding-tutorial");
    await user.keyboard("{ArrowRight}");
    expect(screen.getByText(/2 of 3/)).toBeInTheDocument();
    await user.keyboard("{ArrowLeft}");
    expect(screen.getByText(/1 of 3/)).toBeInTheDocument();
  });
});
