import { useEffect, useState } from "react";
import * as api from "../api";
import type { User } from "../api";

export interface UserModalProps {
  open: boolean;
  onSelect: (username: string) => void;
  onClose: () => void;
}

export function UserModal({ open, onSelect, onClose }: UserModalProps): JSX.Element | null {
  const [users, setUsers]     = useState<User[]>([]);
  const [selection, setSel]   = useState<string>("__new__");
  const [newName, setNewName] = useState("");
  const [error, setError]     = useState<string | null>(null);

  useEffect(() => {
    if (!open) return;
    setError(null);
    void (async () => {
      try {
        const list = await api.listUsers();
        setUsers(list);
        setSel(list.length === 0 ? "__new__" : list[0]!.username);
      } catch (e) {
        setError(`Failed to load users: ${(e as Error).message}`);
      }
    })();
  }, [open]);

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if (e.key === "Escape" && open) onClose();
    }
    window.addEventListener("keydown", onKey);
    return () => { window.removeEventListener("keydown", onKey); };
  }, [open, onClose]);

  if (!open) return null;

  async function submit(): Promise<void> {
    setError(null);
    const username = selection === "__new__" ? newName.trim() : selection;
    if (!username) {
      setError("Username required");
      return;
    }
    try {
      await api.createUser(username, { username });
      onSelect(username);
    } catch (e) {
      setError(`Failed: ${(e as Error).message}`);
    }
  }

  const inputClass =
    "w-full bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent";

  return (
    <div
      className="fixed inset-0 bg-black/60 flex items-center justify-center z-50"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div className="bg-bg-elevated border border-border rounded-md p-6 min-w-[360px] max-w-[480px] flex flex-col gap-4" role="dialog" aria-modal="true">
        <h2 className="text-base font-semibold">Who are you?</h2>
        <p className="text-fg-muted">
          Your name is stored with every change you make so others can see what you've done.
        </p>
        <select
          className={inputClass}
          value={selection}
          onChange={(e) => setSel(e.target.value)}
        >
          {users.map((u) => (
            <option key={u.id} value={u.username}>{u.username}</option>
          ))}
          <option value="__new__">+ New user…</option>
        </select>
        {selection === "__new__" && (
          <input
            className={inputClass}
            type="text"
            placeholder="Enter username"
            value={newName}
            onChange={(e) => setNewName(e.target.value)}
            autoFocus
          />
        )}
        {error && <p className="text-error text-[13px]">{error}</p>}
        <div className="flex justify-end gap-2">
          <button
            className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
            onClick={submit}
          >
            Continue
          </button>
        </div>
      </div>
    </div>
  );
}
