# Ingestion redesign — work ledger

Live worklist for the `ingestion-redesign` branch. This is the moving target; the spec
(`2026-06-22-app-shell-unification-canonical.md`) is the fixed design. When they disagree, the
spec wins on *what to build*, this ledger wins on *what is left and in what order*.

## How agents use this file

- **Pull the topmost unblocked item**, build it, gate it green, commit it, check it off here (`[x]`
  + the commit sha).
- **Append discovered issues** to "Inbox" below as you find them (a bug, a missing dependency, a
  gap in the spec). One line: what + spec anchor (or `no anchor`) + what it blocks or is blocked by.
- **Re-sort after appending.** Keep the numbered list in dependency order: an item never sits above
  something it depends on. Move anything from Inbox into its correct slot in the same pass; leave
  Inbox empty when you can place the item, keep it in Inbox only if its placement is unclear.
- Anchors are `file:line` or `§N` into the spec. Verify the anchor still resolves before relying on it.

## Workstreams (dependency-ordered — do top-down)

### 1. Funnel honesty — completes the just-committed coverage work
- [x] **1a. Backend `image` miss type** (`b78b51fb`) — near-miss extension pass in `/api/fs/manifest`, emits `image` miss + `near` filename; `ManifestUnmatched.near?` + `MissType` extended.
- [x] **1b. D4 zero-coverage hard block** (`2802e29b`) — `zeroCoverageType` first-zero-leg → `No <type> matched` headline + disabled Approve + hidden `/0`; stale fixture fixed + block test added.
- [x] **1c. Autocomplete analysis-dir field** (`2802e29b`) — edit field now `DirectoryPickerField` against `/api/fs/suggest`.

### 2. Scan-failed surface — `p1-failed` (§5.5; `ScanFailedPage.tsx` is a stub)
- [x] **2a. B1 manifest fix** (`bb3f83ac`) — failed-state `fetchManifest` now threads the stored leaf `analysis_dir`; setup_file omitted (no geometry on this surface).
- [x] **2b. p1-failed fidelity** (`47ac79e3`) — warn kicker, "658 of 682 parsed" sentence, two-card layout (Didn't-parse w/ nearest-file pairing + miss groups incl. Image mismatch from 1a; adaptive Test-patterns card w/ live ✓N/N), `Apply all in Configuration →`, two-stage in-place partial-ingest confirm, frosted `↵` chip. *Depends on 1a, 2a.*
  - *Scoping (frontier decision 2026-06-22):* "nearest existing file" pairing is rendered for `image` misses via the backend `near` field (1a). For metadata/integration misses the manifest returns no directory listing, so true nearest-file edit-distance pairing is not computable client-side; those rows show the stem + the expected-missing filename. See Inbox item for the optional backend enrichment.

### 3. Keyboard gaps (§8 status)
- [x] **3a. Series** ↑/↓/Enter/A/⌘Enter (`17fd7874`) — `addSample`/`confirm` ids + Builder wiring.
- [x] **3b. Focus** `P` (`addPeak`) + Enter-apply (`17fd7874`) — `addPeak` id; Enter page-interpreted as apply-candidate.
- [x] **3c-bugfix. Scoping `⌘⇧Z` no longer fires undo** (`8492c196`) — decline Shift-held ⌘Z (Scoping has no redo stack).
- [x] **3d. Enter-on-native-interactive decline gate** (`17fd7874`) — shared `isNativeInteractiveTarget` early-return in every `openFocus`.
- *3c-redo moved into §5 (M6) — see 5d.*

### 4. Dock fidelity (§7, §13)
- [x] **4a. Series dock** (`e018f162`) — §7 grammar: Sample stepper + overlay-matched swatch identity + frosted Focus ↵.
- [x] **4b. Frosted `↵` key-chip** on Focus destinations (Corpus / Loupe / Series; needs `KbKey` frost variant).
- [x] **4c. `CullBar` parameterization** (§13) — optional `actions:{label,onClick,variant}[]` renders in place of the fixed Drop/Keep/Restore/Clear; `variant` ∈ {`accent`,`success`,`ghostInverse`}.

### 5. M6 Configuration edit mode (later-⚙ path)
> **Ledger correction (2026-06-22, on inspection):** the later-⚙ surface (`ConfigurationBody` +
> `GeometryLedger` + `SourcesCard`) ALREADY implements the core: per-field inline geometry override with
> `useUndoStack` undo (`GeometryLedger` aligned `w-32 label | value | SourceChip` grid = the spec's
> `104px | value | source-tag` table), and `handleSourceEdit` does the gated PATCH-then-forced-rescan for
> the three pattern keys (`PATTERN_KEYS`). So 5a + 5b-core are effectively shipped; what remains is a
> **fidelity audit vs `p1-config`** (live, folds into item 7) + the **redo half** (5d).
- [x] **5a. Per-field geometry override + aligned key/value table** — shipped in `GeometryLedger`. Possible fidelity nuance (§5.3 "geometry and sources as one aligned table" — today two aligned cards); defer the merge decision to the live audit (item 7).
- [x] **5b-core. Gated rescan-on-save + undo** — shipped in `ConfigurationBody.handleSourceEdit` + `useUndoStack`. `p1-config` fidelity audit folds into item 7; redo is 5d.
- [ ] **5c. Dock polish sweep** — Configuration has no dock; re-scope at the live audit (likely experiment-shell header/dock polish). *Needs the live walk to characterize.*
- [ ] **5d. Real redo stack (was 3c-redo)** — add a redo stack to `useUndoStack` (forward-appliers per history-entry type), wire Series-Scoping `⌘⇧Z` redo + Configuration `⌘Z`/`⌘⇧Z`. Undo halves exist (Scoping `HistoryEntry`, ConfigurationBody `UndoEntry`); redo needs each entry's forward applier. *Real feature, not wiring.*

### 6. Deferred backend / edges
- [ ] **6a. D6 `series_samples.position` compaction** after merge (§11.4) — renumber `0..n-1` per `series_id` in the merge block + test (shipped leaves gaps).
- [ ] **6b. `last_scanned_at` backend field** for the corpus `N min ago` status (§5.6).
- [ ] **6c. G4 crash-mid-scan** `ingest_status='scanning'` corpus state.
- [ ] **6d. Cross-version SSE deploy sequencing** (§11.5).
- [ ] **6e. toml / `config.jl` removal** (toml deprecated going forward).
- [ ] **6f. e2e** tying an SSE frame to the corpus ProgressBar.
- [ ] **6g. Backend nearest-file enrichment for metadata/integration misses** (§5.5) — widen the 1a near-miss pass to metadata/integration so the manifest emits `near` for case/extension misses on those types (e.g. `RY_BL3_S1780.PRP` vs `*.prp`), letting 2b render true nearest-file pairing for non-image rows. Ratified 1a scope was image-only; this is a deliberate non-blocking follow-up.

### 7. Heavy end-to-end live walk (do LAST)
- [ ] **7. Real ingest** — Approve → full phase-② scan over real `tot_files` data (138 exposures) → grouping unfold → corpus render. The funnel/preview is verified; a real ingest has not been walked. *Depends on 1–2 landing.*

## Inbox (unplaced discoveries)

_Append `- [ ] <what> · <anchor> · blocks/blocked-by <item>` here, then re-sort into the list above._

(empty — placed nearest-file enrichment as 6g)

## Done

_Move checked items here with their commit sha as they land._

(none yet)
