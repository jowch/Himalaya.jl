import { useEffect, useMemo, useRef } from "react";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import {
  useSampleMessages, usePostSampleMessage,
  useComparisonMessages, usePostComparisonMessage,
} from "../queries";
import type { SampleMessage, ComparisonMessage } from "../api";
import { HintText } from "./ui";
import { parseMentions, type MentionToken } from "../lib/renderMentions";
import { useMentionResolution } from "../hooks/useMentionResolution";
import { MentionChip } from "./MentionChip";
import { MentionCompose } from "./MentionCompose";

const CHAT_CARD_FIXTURE = (
  <div className="flex-1 overflow-y-auto px-3 py-3 flex flex-col gap-3 min-h-0">
    <div className="flex flex-col gap-0.5 min-w-0">
      <div className="flex items-baseline gap-2">
        <span className="text-meta">jchen</span>
        <span className="text-fg-dim text-xs">10:15</span>
      </div>
      <p className="text-base font-sans text-fg-muted leading-snug">
        Looks like a clean Pn3m. d-spacing consistent with 70% DOPE.
      </p>
    </div>
    <div className="flex flex-col gap-0.5 min-w-0">
      <div className="flex items-baseline gap-2">
        <span className="text-meta">jchen</span>
        <span className="text-fg-dim text-xs">11:02</span>
      </div>
      <p className="text-base font-sans text-fg-muted leading-snug">
        Re-ran with tighter q-range — score improved to 0.94.
      </p>
    </div>
  </div>
);

/** Either flavour of chat row — both are message-like with author + body. */
type AnyMessage = SampleMessage | ComparisonMessage;

interface ChatCardProps {
  /**
   * Phase 10: explicit entity context. When provided, the hook layer routes
   * to the matching API + cache key. When omitted, the legacy path reads
   * `activeSampleId` from Zustand for backward compatibility with the
   * existing Index/Inspect call sites.
   */
  entityType?: "sample" | "comparison";
  entityId?: number;
}

/**
 * ChatCard — chat thread for a sample or a comparison.
 *
 * Renders a scrollable list of messages with the compose textarea pinned at
 * the bottom of the card. Enter submits; Shift-Enter inserts a newline.
 *
 * The component supports two contexts:
 *   - `<ChatCard />` — sample-scoped (legacy): reads `activeSampleId` from
 *     Zustand. Used on the Index and Inspect pages.
 *   - `<ChatCard entityType="comparison" entityId={N} />` — comparison-scoped:
 *     used on review-mode `ComparePage`.
 *
 * Both flavours route to the right API + cache key via the queries layer.
 * Edit-mode Compare does NOT mount the chat (chat is review-only per spec).
 */
export function ChatCard(props: ChatCardProps): JSX.Element {
  const activeSampleId = useAppState((s) => s.activeSampleId);

  // Resolve effective context: explicit props win, fall back to active sample.
  const entityType: "sample" | "comparison" = props.entityType ?? "sample";
  const entityId: number | undefined =
    props.entityId !== undefined
      ? props.entityId
      : entityType === "sample"
        ? activeSampleId
        : undefined;

  if (entityId === undefined) {
    return (
      <Frame>
        <div className="flex-1 flex items-center justify-center p-4">
          <HintText>
            {entityType === "comparison"
              ? "Pick a comparison to start a conversation."
              : "Pick a sample to start a conversation."}
          </HintText>
        </div>
      </Frame>
    );
  }

  return entityType === "comparison" ? (
    <ComparisonThread comparisonId={entityId} />
  ) : (
    <SampleThread sampleId={entityId} />
  );
}

function SampleThread({ sampleId }: { sampleId: number }): JSX.Element {
  const username = useAppState((s) => s.username);
  const messagesQ = useSampleMessages(sampleId);
  const postMsg   = usePostSampleMessage(sampleId);
  return (
    <Frame>
      <MessageList messages={messagesQ.data ?? []} isLoading={messagesQ.isLoading} />
      <MentionCompose
        disabled={username === undefined || postMsg.isPending}
        onSubmit={(body) => postMsg.mutate(body)}
      />
    </Frame>
  );
}

function ComparisonThread({ comparisonId }: { comparisonId: number }): JSX.Element {
  const username = useAppState((s) => s.username);
  const messagesQ = useComparisonMessages(comparisonId);
  const postMsg   = usePostComparisonMessage(comparisonId);
  return (
    <Frame>
      <MessageList messages={messagesQ.data ?? []} isLoading={messagesQ.isLoading} />
      <MentionCompose
        disabled={username === undefined || postMsg.isPending}
        onSubmit={(body) => postMsg.mutate(body)}
      />
    </Frame>
  );
}

function Frame({ children }: { children: React.ReactNode }): JSX.Element {
  return (
    <div
      data-testid="chat-card"
      className="flex flex-col h-full min-h-0 overflow-hidden"
    >
      {children}
    </div>
  );
}

interface MessageListProps {
  messages: AnyMessage[];
  isLoading: boolean;
}

function MessageList({ messages, isLoading }: MessageListProps): JSX.Element {
  const scrollRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const el = scrollRef.current;
    if (el) el.scrollTop = el.scrollHeight;
  }, [messages.length]);

  return (
    <Skeleton
      name="chat-card"
      className="flex-1 min-h-0 flex flex-col"
      loading={isLoading}
      stagger={50}
      transition={200}
      fixture={CHAT_CARD_FIXTURE}
      fallback={<div className="flex-1 overflow-y-auto px-3 py-3"><HintText>Loading…</HintText></div>}
    >
      {messages.length === 0 ? (
        <div className="flex-1 overflow-y-auto px-3 py-3">
          <HintText>No notes yet. Start a conversation below.</HintText>
        </div>
      ) : (
        <div
          ref={scrollRef}
          className="flex-1 overflow-y-auto px-3 py-3 flex flex-col gap-3 min-h-0"
          data-testid="chat-message-list"
        >
          {messages.map((m) => (
            <MessageRow key={m.id} msg={m} />
          ))}
        </div>
      )}
    </Skeleton>
  );
}

function MessageRow({ msg }: { msg: AnyMessage }): JSX.Element {
  const authorLabel   = msg.author ?? "deleted user";
  const authorDeleted = msg.author == null;
  const segments      = useMemo(() => parseMentions(msg.body), [msg.body]);
  const mentions      = useMemo(
    () => segments.filter((s): s is MentionToken => s.kind === "mention"),
    [segments],
  );
  const resolutionMap = useMentionResolution(mentions);

  return (
    <div className="flex flex-col gap-0.5 min-w-0" data-testid={`chat-message-${msg.id}`}>
      <div className="flex items-baseline gap-2">
        <span className={authorDeleted ? "text-meta text-fg-dim italic" : "text-meta"}>
          {authorLabel}
        </span>
        <span className="text-fg-dim text-xs">{formatTime(msg.created_at)}</span>
      </div>
      <p className="text-base font-sans text-fg-muted leading-snug break-words whitespace-pre-wrap">
        {segments.map((seg, i) => {
          if (seg.kind === "text") return <span key={i}>{seg.text}</span>;
          const key   = `${seg.type}:${seg.id}`;
          const token = seg.hash !== undefined
            ? `[[${seg.type}:${seg.id}@${seg.hash}]]`
            : `[[${seg.type}:${seg.id}]]`;
          const entry = resolutionMap.get(key) ?? "loading";
          // Phase 10 — pass the parsed token hash through to the chip so
          // comparison drift can render. Other mention types ignore it.
          const extraProps = seg.hash !== undefined ? { tokenHash: seg.hash } : {};
          return (
            <MentionChip
              key={i}
              resolved={entry}
              originalText={token}
              {...extraProps}
            />
          );
        })}
      </p>
    </div>
  );
}

function formatTime(iso: string): string {
  // Accept "YYYY-MM-DD HH:MM:SS" or ISO-8601. Display HH:MM for today, date+time otherwise.
  const clean = iso.includes("T") ? iso : iso.replace(" ", "T") + "Z";
  const d = new Date(clean);
  if (Number.isNaN(d.getTime())) return iso;
  const hm = d.toLocaleTimeString(undefined, { hour: "2-digit", minute: "2-digit", hour12: false });
  const today = new Date();
  if (d.toDateString() === today.toDateString()) return hm;
  return `${d.toLocaleDateString(undefined, { month: "short", day: "numeric" })} ${hm}`;
}
