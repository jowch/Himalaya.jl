import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { ChatCard, MessageList } from "../src/components/ChatCard";
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
      expect(postSpy).toHaveBeenCalledWith(3, "hello there", expect.objectContaining({ username: "carol" }));
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

  it("renders a peak mention chip inline", async () => {
    const PEAK: api.Peak = {
      id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
      sharpness: 0.3, source: "auto", excluded: false,
    };
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([
      { id: 1, sample_id: 3, author_id: 1, author: "alice",
        body: "see [[peak:42]]", created_at: "2026-04-24 10:00:00" },
    ]);
    vi.spyOn(api, "getPeak").mockResolvedValue(PEAK);
    renderWithProviders(<ChatCard />);
    expect(await screen.findByText(/q = 1\.223/)).toBeInTheDocument();
  });

  // Regression #124: gate `<Skeleton>` on `query.isLoading`, NOT `isPending`.
  // `isPending` is true for disabled queries (`enabled: false`) — gating on it
  // would flash a skeleton when the consumer caller hasn't yet picked an
  // entity. The discriminating state is `{ isPending: true, isLoading: false }`,
  // which is unreachable through the public `<ChatCard>` flow (the early-return
  // for undefined entityId masks it), so we exercise `<MessageList>` directly.
  it("MessageList: skeleton hidden when isLoading=false (#124)", () => {
    renderWithProviders(<MessageList messages={[]} isLoading={false} />);
    expect(screen.getByText(/no notes yet/i)).toBeInTheDocument();
    expect(screen.queryByText(/loading…/i)).toBeNull();
  });

  it("MessageList: skeleton shown when isLoading=true (#124)", () => {
    renderWithProviders(<MessageList messages={[]} isLoading={true} />);
    expect(screen.queryByText(/no notes yet/i)).toBeNull();
  });

  it("renders dead chip when mention entity returns 404", async () => {
    vi.spyOn(api, "listSampleMessages").mockResolvedValue([
      { id: 2, sample_id: 3, author_id: 1, author: "alice",
        body: "old ref [[index:999]]", created_at: "2026-04-24 10:00:00" },
    ]);
    vi.spyOn(api, "getIndex").mockRejectedValue(new api.ApiError(404, "not found", null));
    renderWithProviders(<ChatCard />);
    await waitFor(() => {
      const chip = document.querySelector("[data-mention-state='dead']");
      expect(chip).not.toBeNull();
    });
  });
});

// ─── Comparison context (Phase 10, Task 10.1) ─────────────────────────────
//
// `ChatCard` accepts an explicit entity context so it can host BOTH the
// per-sample chat (Index/Inspect pages) and the per-comparison chat
// (review-mode Compare page). The hook layer routes to the right
// API by inspecting `entityType`; the legacy zero-prop call site still
// reads `activeSampleId` from Zustand for backward compatibility.

describe("<ChatCard entityType='comparison'>", () => {
  const COMPARISON_MSGS: api.ComparisonMessage[] = [
    {
      id: 11, comparison_id: 7, author_id: 1, author: "alice",
      body: "first comparison thought", created_at: "2026-05-02 09:00:00",
    },
    {
      id: 12, comparison_id: 7, author_id: 2, author: "bob",
      body: "second comparison thought", created_at: "2026-05-02 09:05:00",
    },
  ];

  beforeEach(() => {
    useAppState.setState({
      activeSampleId: undefined,
      username: "carol",
      activeExperimentId: 1,
    });
  });

  it("calls listComparisonMessages and renders messages", async () => {
    const listSpy = vi.spyOn(api, "listComparisonMessages").mockResolvedValue(COMPARISON_MSGS);
    const sampleSpy = vi.spyOn(api, "listSampleMessages").mockResolvedValue([]);
    renderWithProviders(<ChatCard entityType="comparison" entityId={7} />);
    expect(await screen.findByText("first comparison thought")).toBeInTheDocument();
    expect(screen.getByText("second comparison thought")).toBeInTheDocument();
    expect(listSpy).toHaveBeenCalledWith(7);
    expect(sampleSpy).not.toHaveBeenCalled();
  });

  it("posts via postComparisonMessage on Enter", async () => {
    vi.spyOn(api, "listComparisonMessages").mockResolvedValue([]);
    const postSpy = vi.spyOn(api, "postComparisonMessage").mockResolvedValue({
      id: 99, comparison_id: 7, author_id: 3, author: "carol",
      body: "compare note", created_at: "2026-05-02 10:05:00",
    });
    const sampleSpy = vi.spyOn(api, "postSampleMessage");
    const user = userEvent.setup();
    renderWithProviders(<ChatCard entityType="comparison" entityId={7} />);
    const compose = await screen.findByTestId("chat-compose");
    await user.click(compose);
    await user.keyboard("compare note");
    await user.keyboard("{Enter}");
    await waitFor(() => {
      expect(postSpy).toHaveBeenCalledWith(7, "compare note", expect.objectContaining({ username: "carol" }));
    });
    expect(sampleSpy).not.toHaveBeenCalled();
  });

  it("renders empty state with comparison-specific copy", async () => {
    vi.spyOn(api, "listComparisonMessages").mockResolvedValue([]);
    renderWithProviders(<ChatCard entityType="comparison" entityId={7} />);
    expect(await screen.findByText(/no notes yet/i)).toBeInTheDocument();
  });
});

describe("<ChatCard entityType='sample'>", () => {
  it("explicit sample context routes to listSampleMessages", async () => {
    const sampleSpy = vi.spyOn(api, "listSampleMessages").mockResolvedValue([
      { id: 1, sample_id: 5, author_id: 1, author: "alice",
        body: "hi sample", created_at: "2026-04-24 10:00:00" },
    ]);
    const compSpy = vi.spyOn(api, "listComparisonMessages").mockResolvedValue([]);
    useAppState.setState({ activeSampleId: undefined, username: "carol", activeExperimentId: 1 });
    renderWithProviders(<ChatCard entityType="sample" entityId={5} />);
    expect(await screen.findByText("hi sample")).toBeInTheDocument();
    expect(sampleSpy).toHaveBeenCalledWith(5);
    expect(compSpy).not.toHaveBeenCalled();
  });
});
