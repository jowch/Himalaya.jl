import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { ChatCard } from "../src/components/ChatCard";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const MESSAGES: api.SampleMessage[] = [
  { id: 1, sample_id: 3, author_id: 1, author: "alice", body: "looks cubic",    created_at: "2026-04-24 10:00:00" },
  { id: 2, sample_id: 3, author_id: 2, author: "bob",   body: "Im3m a=19.3 nm", created_at: "2026-04-24 10:02:00" },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: 3,
    username: "carol",
    activeExperimentId: 1,
  });
});

describe("<ChatCard>", () => {
  it("renders a hint when no sample is selected", () => {
    useAppState.setState({ activeSampleId: undefined });
    renderWithProviders(<ChatCard />);
    expect(screen.getByText(/pick a sample/i)).toBeInTheDocument();
  });

  it("renders messages for the active sample", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue(MESSAGES);
    renderWithProviders(<ChatCard />);
    expect(await screen.findByText("looks cubic")).toBeInTheDocument();
    expect(screen.getByText("Im3m a=19.3 nm")).toBeInTheDocument();
    expect(screen.getByText("alice")).toBeInTheDocument();
    expect(screen.getByText("bob")).toBeInTheDocument();
  });

  it("renders empty state when no messages yet", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    renderWithProviders(<ChatCard />);
    expect(await screen.findByText(/no notes yet/i)).toBeInTheDocument();
  });

  it("pressing Enter submits via usePostSampleMessage", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    const postSpy = vi.spyOn(api, "postSampleMessage").mockResolvedValue({
      id: 99, sample_id: 3, author_id: 3, author: "carol", body: "hello there",
      created_at: "2026-04-24 10:05:00",
    });
    const user = userEvent.setup();
    renderWithProviders(<ChatCard />);
    const compose = await screen.findByTestId("chat-compose");
    await user.click(compose);
    await user.keyboard("hello there");
    await user.keyboard("{Enter}");
    await waitFor(() => {
      expect(postSpy).toHaveBeenCalledWith(3, "hello there", { username: "carol" });
    });
  });

  it("Shift+Enter inserts newline and does NOT submit", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    const postSpy = vi.spyOn(api, "postSampleMessage");
    const user = userEvent.setup();
    renderWithProviders(<ChatCard />);
    const compose = await screen.findByTestId("chat-compose");
    await user.click(compose);
    await user.keyboard("line1");
    await user.keyboard("{Shift>}{Enter}{/Shift}");
    await user.keyboard("line2");
    expect(postSpy).not.toHaveBeenCalled();
    expect((compose as HTMLTextAreaElement).value).toContain("line1");
    expect((compose as HTMLTextAreaElement).value).toContain("line2");
  });

  it("does not submit an empty/whitespace-only message", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    const postSpy = vi.spyOn(api, "postSampleMessage");
    const user = userEvent.setup();
    renderWithProviders(<ChatCard />);
    const compose = await screen.findByTestId("chat-compose");
    await user.click(compose);
    await user.keyboard("   ");
    await user.keyboard("{Enter}");
    expect(postSpy).not.toHaveBeenCalled();
  });

  it("renders 'deleted user' when author is null", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([
      { id: 7, sample_id: 3, author_id: null, author: null,
        body: "orphaned note", created_at: "2026-04-24 11:00:00" },
    ]);
    renderWithProviders(<ChatCard />);
    expect(await screen.findByText("orphaned note")).toBeInTheDocument();
    expect(screen.getByText(/deleted user/i)).toBeInTheDocument();
  });

  it("does not submit when username is not set", async () => {
    useAppState.setState({ username: undefined });
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    const postSpy = vi.spyOn(api, "postSampleMessage");
    const user = userEvent.setup();
    renderWithProviders(<ChatCard />);
    const compose = await screen.findByTestId("chat-compose");
    await user.click(compose);
    await user.keyboard("hello");
    await user.keyboard("{Enter}");
    expect(postSpy).not.toHaveBeenCalled();
  });
});
