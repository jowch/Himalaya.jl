import { useEffect, useMemo, useRef } from "react";
import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import { useSampleMessages, usePostSampleMessage } from "../queries";
import type { SampleMessage } from "../api";
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

/**
 * ChatCard — per-sample notebook/chat log.
 *
 * Renders a scrollable list of messages with the compose textarea pinned at
 * the bottom of the card. Enter submits; Shift-Enter inserts a newline.
 */
export function ChatCard(): JSX.Element {
  const sampleId = useAppState((s) => s.activeSampleId);
  const username = useAppState((s) => s.username);

  const messagesQ = useSampleMessages(sampleId);
  const postMsg   = usePostSampleMessage(sampleId ?? 0);

  if (sampleId === undefined) {
    return (
      <Frame>
        <div className="flex-1 flex items-center justify-center p-4">
          <HintText>Pick a sample to start a conversation.</HintText>
        </div>
      </Frame>
    );
  }

  return (
    <Frame>
      <MessageList messages={messagesQ.data ?? []} isPending={messagesQ.isPending} />
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
  messages: SampleMessage[];
  isPending: boolean;
}

function MessageList({ messages, isPending }: MessageListProps): JSX.Element {
  const scrollRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const el = scrollRef.current;
    if (el) el.scrollTop = el.scrollHeight;
  }, [messages.length]);

  return (
    <Skeleton
      name="chat-card"
      loading={isPending}
      stagger={50}
      transition={200}
      fixture={CHAT_CARD_FIXTURE}
      fallback={<div className="flex-1 overflow-y-auto px-3 py-3"><HintText>Loading…</HintText></div>}
    >
      {messages.length === 0 ? (
        <div className="flex-1 overflow-y-auto px-3 py-3">
          <HintText>No notes yet. Start a conversation about this sample below.</HintText>
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

function MessageRow({ msg }: { msg: SampleMessage }): JSX.Element {
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
          const token = `[[${seg.type}:${seg.id}]]`;
          const entry = resolutionMap.get(key) ?? "loading";
          return (
            <MentionChip key={i} resolved={entry} originalText={token} />
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

