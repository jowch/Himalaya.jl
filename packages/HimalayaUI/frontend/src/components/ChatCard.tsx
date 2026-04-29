import { useEffect, useRef, useState } from "react";
import { useAppState } from "../state";
import { useSampleMessages, usePostSampleMessage } from "../queries";
import type { SampleMessage } from "../api";
import { HintText } from "./ui";

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
      <Compose
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
    // Scroll to bottom when messages change
    const el = scrollRef.current;
    if (el) el.scrollTop = el.scrollHeight;
  }, [messages.length]);

  if (isPending) {
    return (
      <div className="flex-1 overflow-y-auto px-3 py-3">
        <HintText>Loading…</HintText>
      </div>
    );
  }
  if (messages.length === 0) {
    return (
      <div className="flex-1 overflow-y-auto px-3 py-3">
        <HintText>No notes yet. Start a conversation about this sample below.</HintText>
      </div>
    );
  }

  return (
    <div
      ref={scrollRef}
      className="flex-1 overflow-y-auto px-3 py-3 flex flex-col gap-3 min-h-0"
      data-testid="chat-message-list"
    >
      {messages.map((m) => (
        <MessageRow key={m.id} msg={m} />
      ))}
    </div>
  );
}

function MessageRow({ msg }: { msg: SampleMessage }): JSX.Element {
  const authorLabel = msg.author ?? "deleted user";
  const authorDeleted = msg.author == null;
  return (
    <div className="flex flex-col gap-0.5 min-w-0" data-testid={`chat-message-${msg.id}`}>
      <div className="flex items-baseline gap-2">
        <span className={authorDeleted
          ? "text-meta text-fg-dim italic"
          : "text-meta"}>
          {authorLabel}
        </span>
        <span className="text-fg-dim text-xs">{formatTime(msg.created_at)}</span>
      </div>
      <p className="text-base font-sans text-fg-muted leading-snug break-words whitespace-pre-wrap">
        {msg.body}
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

interface ComposeProps {
  disabled: boolean;
  onSubmit: (body: string) => void;
}

function Compose({ disabled, onSubmit }: ComposeProps): JSX.Element {
  const [text, setText] = useState("");
  const ref = useRef<HTMLTextAreaElement>(null);

  const trySubmit = (): void => {
    const trimmed = text.trim();
    if (!trimmed || disabled) return;
    onSubmit(trimmed);
    setText("");
  };

  const onKeyDown = (e: React.KeyboardEvent<HTMLTextAreaElement>): void => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      trySubmit();
    }
  };

  return (
    <div className="flex-shrink-0 border-t border-border bg-bg px-2.5 py-2">
      <textarea
        ref={ref}
        value={text}
        onChange={(e) => setText(e.target.value)}
        onKeyDown={onKeyDown}
        rows={2}
        placeholder={disabled ? "Sign in to post…" : "Write a note…"}
        data-testid="chat-compose"
        className="w-full resize-none bg-transparent text-fg text-base font-sans
                   placeholder:text-fg-dim outline-0 border-0"
      />
      <div className="flex items-center justify-between text-xs text-fg-dim">
        <span>
          <kbd className="border border-border rounded px-1">⏎</kbd> send
          {" · "}
          <kbd className="border border-border rounded px-1">⇧⏎</kbd> newline
        </span>
      </div>
    </div>
  );
}
