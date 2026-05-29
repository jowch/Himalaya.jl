import { useEffect, useRef, useState } from "react";
import * as api from "../api";
import type { User } from "../api";
import { useAppState } from "../state";
import { Button } from "./ui";
import { useFocusTrap } from "../hooks/useFocusTrap";

/**
 * OnboardingFlow — shown when no username is in persisted state.
 *
 * Step 1 (name):
 *   - Drop-down of existing users; "+ New user…" opens a name input.
 *   - Picking an existing user short-circuits the tutorial.
 *
 * Step 2 (tutorial):
 *   - Four slides introducing title-button, three cards, keyboard shortcuts,
 *     and active-set interaction. Shown only to brand-new users whose
 *     `tutorialSeen` flag is still false.
 */
const TUTORIAL_SLIDES: readonly { title: string; body: string }[] = [
  {
    title: "Scope lives up top",
    body:
      "The title-button above the plot shows which experiment and sample you're looking at. " +
      "Click it — or press / — to change either.",
  },
  {
    title: "Three panels, one screen",
    body:
      "Chat is on the left — jot observations as you go, and type @ to mention a peak, " +
      "index, exposure, or sample inline. The plot is in the middle. " +
      "Index choices are on the right: hover a candidate to preview its peaks.",
  },
  {
    title: "Move quickly between samples",
    body:
      "Press , and . to jump to the previous or next sample. " +
      "Use ← → to switch between the Inspect, Index, and Compare tabs. " +
      "Press / any time to open the picker.",
  },
  {
    title: "Curate the active set",
    body:
      "Click + on a candidate to add it to the active set; click − to take it back out. " +
      "Changes save automatically.",
  },
];

type Phase = "name" | "tutorial";

/** Allowed handle characters: letters, digits, underscore. No @, no spaces. */
const HANDLE_RE = /^[A-Za-z0-9_]+$/;

interface PendingUser {
  username: string;
  firstName?: string | undefined;
  lastName?: string | undefined;
}

export function OnboardingFlow(): JSX.Element | null {
  const username      = useAppState((s) => s.username);
  const tutorialSeen  = useAppState((s) => s.tutorialSeen);
  const setUser       = useAppState((s) => s.setUser);
  const setTutorialSeen = useAppState((s) => s.setTutorialSeen);

  const [phase, setPhase] = useState<Phase>("name");
  const [users, setUsers] = useState<User[]>([]);
  const [selection, setSelection] = useState<string>("__new__");
  const [newHandle, setNewHandle] = useState("");
  const [newFirst, setNewFirst]   = useState("");
  const [newLast, setNewLast]     = useState("");
  const [error, setError] = useState<string | null>(null);
  const [slide, setSlide] = useState(0);
  // Hold the chosen user locally until the tutorial is dismissed; committing
  // to the store would unmount this component and skip the tutorial.
  const [pending, setPending] = useState<PendingUser | null>(null);

  useEffect(() => {
    if (username !== undefined) return;
    setError(null);
    setPhase("name");
    setSlide(0);
    setPending(null);
    void (async () => {
      try {
        const list = await api.listUsers();
        setUsers(list);
        setSelection(list.length === 0 ? "__new__" : list[0]!.username);
      } catch (e) {
        setError(`Failed to load users: ${(e as Error).message}`);
      }
    })();
  }, [username]);

  if (username !== undefined && phase !== "tutorial") return null;

  const onSubmitName = async (): Promise<void> => {
    setError(null);
    const isNew = selection === "__new__";
    // Tolerate a leading @ even though the visual prefix means it shouldn't be there.
    const handle = isNew ? newHandle.trim().replace(/^@/, "") : selection;
    if (!handle) { setError("Handle required"); return; }
    if (isNew && !HANDLE_RE.test(handle)) {
      setError("Handle may only contain letters, digits, and underscores");
      return;
    }
    const fields: { first_name?: string; last_name?: string } = {};
    if (isNew && newFirst.trim()) fields.first_name = newFirst.trim();
    if (isNew && newLast.trim())  fields.last_name  = newLast.trim();
    try {
      const created = await api.createUser(handle, fields, { username: handle });
      const next: PendingUser = {
        username: created.username,
        ...(created.first_name ? { firstName: created.first_name } : {}),
        ...(created.last_name  ? { lastName:  created.last_name  } : {}),
      };
      // Brand-new user AND hasn't seen tutorial? → stash, show tutorial.
      if (isNew && !tutorialSeen) {
        setPending(next);
        setPhase("tutorial");
      } else {
        setUser(next);
      }
    } catch (e) {
      setError(`Failed: ${(e as Error).message}`);
    }
  };

  const closeTutorial = (): void => {
    setTutorialSeen(true);
    if (pending !== null) setUser(pending);
    setPending(null);
    setPhase("name");
  };

  const onTutorialKey = (e: React.KeyboardEvent<HTMLDivElement>): void => {
    if (e.key === "Escape") { closeTutorial(); return; }
    if (e.key === "ArrowRight") { setSlide((s) => Math.min(TUTORIAL_SLIDES.length - 1, s + 1)); return; }
    if (e.key === "ArrowLeft")  { setSlide((s) => Math.max(0, s - 1)); }
  };

  return (
    <div
      data-testid="onboarding-overlay"
      className="fixed inset-0 z-50 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm
                 flex items-center justify-center"
      role="presentation"
    >
      {phase === "name" && (
        <NameStep
          users={users}
          selection={selection}
          onSelection={setSelection}
          newHandle={newHandle}
          onNewHandle={setNewHandle}
          newFirst={newFirst}
          onNewFirst={setNewFirst}
          newLast={newLast}
          onNewLast={setNewLast}
          error={error}
          onSubmit={onSubmitName}
        />
      )}
      {phase === "tutorial" && (
        <TutorialStep
          slideIdx={slide}
          onPrev={() => setSlide((s) => Math.max(0, s - 1))}
          onNext={() => setSlide((s) => Math.min(TUTORIAL_SLIDES.length - 1, s + 1))}
          onDone={closeTutorial}
          onKeyDown={onTutorialKey}
        />
      )}
    </div>
  );
}

// ── Name step ─────────────────────────────────────────────────────────────

interface NameStepProps {
  users: User[];
  selection: string;
  onSelection: (s: string) => void;
  newHandle: string;
  onNewHandle: (s: string) => void;
  newFirst: string;
  onNewFirst: (s: string) => void;
  newLast: string;
  onNewLast: (s: string) => void;
  error: string | null;
  onSubmit: () => void;
}

function NameStep({
  users, selection, onSelection,
  newHandle, onNewHandle, newFirst, onNewFirst, newLast, onNewLast,
  error, onSubmit,
}: NameStepProps): JSX.Element {
  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, true);
  const inputClass =
    "w-full bg-paper border border-hair-strong rounded-md px-2 py-1 " +
    "focus:outline focus:outline-1 focus:outline-accent focus:border-accent";

  const displayName = (u: User): string => {
    const full = [u.first_name, u.last_name].filter(Boolean).join(" ");
    return full ? `${full} (@${u.username})` : u.username;
  };

  const previewLabel = (first: string, last: string, handle: string): string => {
    const name = [first.trim(), last.trim()].filter(Boolean).join(" ");
    const h = handle.trim().replace(/^@/, "");
    if (!name && !h) return " "; // reserve a line of space so layout doesn't jump
    if (!h)          return name;
    if (!name)       return `@${h}`;
    return `${name} @${h}`;
  };

  return (
    <div
      ref={dialogRef}
      data-testid="onboarding-name"
      role="dialog"
      aria-modal="true"
      className="bg-plate border border-hair-strong rounded-lg p-6
                 min-w-[360px] max-w-[480px] flex flex-col gap-4"
      onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); onSubmit(); } }}
    >
      <h2 className="text-base font-semibold text-ink">Who are you?</h2>
      <p className="text-ink-soft text-base">
        Your name is stored with every change so others can see what you've done.
      </p>
      <select
        className={inputClass}
        value={selection}
        onChange={(e) => onSelection(e.target.value)}
        data-testid="onboarding-user-select"
      >
        {users.map((u) => (
          <option key={u.id} value={u.username}>{displayName(u)}</option>
        ))}
        <option value="__new__">+ New user…</option>
      </select>
      {selection === "__new__" && (
        <div className="flex flex-col gap-2">
          <input
            className={inputClass}
            type="text"
            placeholder="First name"
            value={newFirst}
            onChange={(e) => onNewFirst(e.target.value)}
            autoFocus
            data-testid="onboarding-new-first"
          />
          <input
            className={inputClass}
            type="text"
            placeholder="Last name"
            value={newLast}
            onChange={(e) => onNewLast(e.target.value)}
            data-testid="onboarding-new-last"
          />
          <div className="flex items-stretch border border-hair-strong rounded-md
                          bg-paper overflow-hidden
                          focus-within:outline focus-within:outline-1
                          focus-within:outline-accent focus-within:border-accent">
            <span
              aria-hidden="true"
              className="px-2 py-1 text-ink-soft bg-plate border-r border-hair-strong
                         select-none"
            >@</span>
            <input
              className="flex-1 bg-transparent px-2 py-1 focus:outline-none"
              type="text"
              placeholder="handle"
              value={newHandle}
              onChange={(e) => onNewHandle(e.target.value.replace(/^@/, ""))}
              data-testid="onboarding-new-handle"
            />
          </div>
          <p
            data-testid="onboarding-preview"
            className="text-ink-soft text-sm pl-1 min-h-[1.25rem]"
          >
            {previewLabel(newFirst, newLast, newHandle)}
          </p>
        </div>
      )}
      {error && <p className="text-error text-base">{error}</p>}
      <div className="flex justify-end">
        <Button variant="primary" onClick={onSubmit} data-testid="onboarding-continue">
          Continue
        </Button>
      </div>
    </div>
  );
}

// ── Tutorial step ─────────────────────────────────────────────────────────

interface TutorialStepProps {
  slideIdx: number;
  onPrev: () => void;
  onNext: () => void;
  onDone: () => void;
  onKeyDown: (e: React.KeyboardEvent<HTMLDivElement>) => void;
}

function TutorialStep({
  slideIdx, onPrev, onNext, onDone, onKeyDown,
}: TutorialStepProps): JSX.Element {
  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, true);
  const slide = TUTORIAL_SLIDES[slideIdx]!;
  const isLast = slideIdx === TUTORIAL_SLIDES.length - 1;
  return (
    <div
      ref={dialogRef}
      data-testid="onboarding-tutorial"
      role="dialog"
      aria-modal="true"
      tabIndex={-1}
      onKeyDown={onKeyDown}
      className="bg-plate border border-hair-strong rounded-lg p-7
                 min-w-[420px] max-w-[520px] flex flex-col gap-4 outline-0"
    >
      <div className="text-xs uppercase tracking-widest text-ink-faint">
        Welcome · {slideIdx + 1} of {TUTORIAL_SLIDES.length}
      </div>
      <h2 className="text-lg font-semibold text-ink">{slide.title}</h2>
      <p className="text-ink-soft text-base leading-relaxed">{slide.body}</p>
      <div className="flex items-center justify-between pt-2">
        <div className="flex gap-1">
          {TUTORIAL_SLIDES.map((_, i) => (
            <span
              key={i}
              className={
                "w-1.5 h-1.5 rounded-full " +
                (i === slideIdx ? "bg-accent" : "bg-hair-strong")
              }
            />
          ))}
        </div>
        <div className="flex gap-2">
          {slideIdx > 0 && (
            <Button variant="ghost" onClick={onPrev} data-testid="tutorial-prev">
              Back
            </Button>
          )}
          {!isLast ? (
            <Button variant="primary" onClick={onNext} data-testid="tutorial-next">
              Next
            </Button>
          ) : (
            <Button variant="primary" onClick={onDone} data-testid="tutorial-done">
              Got it
            </Button>
          )}
        </div>
      </div>
    </div>
  );
}
