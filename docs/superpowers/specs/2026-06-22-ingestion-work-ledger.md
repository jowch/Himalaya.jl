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

### 7. Heavy end-to-end live walk (DONE — surfaced P0/P1/P2)
- [x] **7. Real ingest live walk** (mini-1p7m → `mini-clean.db`, 138 exposures, 35 samples, 3 loads) — WALKED end to end. ✅ Funnel (picker pre-flight ✓✓ → config coverage 138 Image / 138/138 Metadata / 138/138 Integration, no D4 block → leaf analysis_dir resolved → geometry from real `setup_info` → Approve). ✅ Grouping live-unfold + split-flag segmentation (real stage-position jumps) + Confirm-gate. ✅ Corpus render (5-stat masthead 3·35·138·5.9h·2, amber review banner, sheet with real detector thumbnails + slot chips, §7 dock). 🐛 Found **P0 .dat-pattern (analyze indexes nothing)**, **P1 route→404 (fixed)**, **P2 scanning-no-auto-complete** — all in Inbox. *The walk did its job: caught a P0 that 1000+ green tests + the funnel preview all missed.*

### 8. P0/P1 from the live walk (do NEXT — blocks a usable real ingest)
- [x] **8a. P0 — integration-`.dat` resolution** (`365ba492`) — `config_from_db` now prefers the per-experiment pattern COLUMNS over the deprecated blob; root cause was the builtin `{name}.dat` fallback. Live-verified: 173 indices / 333 auto_peaks / 0 dat-errors on a real re-ingest; suite green.
- [x] **8b. P1 — grouping Confirm/Back → 404** (`9abb40c5`) — absolute route in `AppRoutes.tsx`.
- [x] **8c. P2 — scanning UI auto-complete** (`95ef696f` + review fix `0e136f58`) — root cause was a stale live overlay: the at-most-once `/api/events` stream drops the terminal `ingest_complete` frame on an EventSource reconnect (benign EPIPE), leaving `ingestInFlight` stuck. Fixed with one shared `effectiveIngestStatus` rule (**a TERMINAL persisted state always wins**, self-healing both scan kinds) + a scoped `useExperiment` `refetchInterval`. frontend-reviewer caught 2 P1s in the first cut (verified vs `routes_experiments.jl`: the rescan route sets the row scanning→complete just like an initial scan, so the original "analyzing-overlay-wins" branch reintroduced the bug for rescans, and the corpus surface choice keyed on persisted-`scanning` which both routes set) → fixed in `0e136f58` (drop the branch; corpus surface gates on overlay phase + `last_scanned_at`). Live-verified on a fresh DB: scan completes → grouping auto-flips to review + enabled Confirm, no reload.

## Inbox (unplaced discoveries)

_Append `- [ ] <what> · <anchor> · blocks/blocked-by <item>` here, then re-sort into the list above._

**All three item-7 live-walk findings were PLACED into §8 (do not re-place); detail kept here for provenance:**
- **(→ 8a, RESOLVED `365ba492`)** Real ingest indexed ZERO samples. First-written diagnosis pointed at `analyze_exposure!`; the actual root cause was narrower: `config_from_db` read the deprecated TOML `config` blob, so analyze got the builtin `{name}.dat` instead of the stored `{name}_tot.dat` column. Real files `agbe_S1963_tot.dat`; scan looked for `agbe_S1963.dat` → "dat file not found" ×138 → nothing indexed. **The funnel preview MASKED it** (`/api/fs/manifest` uses the pattern → 138/138). Fixed by making `config_from_db` prefer the pattern columns; live-verified 173 indices / 333 auto_peaks / 0 errors.
- **(→ 8b, RESOLVED `9abb40c5`)** Confirm-groups / Back navigated to bare `/corpus` (404) — relative `../corpus` from the top-level `/grouping` takeover. Fixed to absolute `/experiments/:id/corpus`.
- **(→ 8c, RESOLVED `95ef696f`)** Scanning UI did not auto-flip on scan completion. Root cause confirmed on the second fresh-DB walk: the at-most-once `/api/events` stream drops the terminal `ingest_complete` frame on an EventSource reconnect (the benign EPIPE), leaving `ingestInFlight` stuck `"scanning"`; the derivation preferred that overlay unconditionally. Fixed via `effectiveIngestStatus` + scoped `useExperiment` `refetchInterval`. Live-verified.

(prior: nearest-file enrichment placed as 6g)

**Second fresh-DB live walk (2026-06-22, new mini-fresh.db) — 3 NEW first-run findings (only reachable from a truly empty DB):**
- [x] **9a. P2 — `serve` fatally errors on a missing DB** (`9a0823a9`) — resolved the bootstrap chicken-and-egg: `serve` now creates+migrates an empty DB via `open_db` when absent, so a new user lands on the empty experiments home and creates one via the funnel. A loud `@warn` fires on bootstrap so a typo'd `HIMALAYA_DB_PATH` can't silently mask existing data. The `analyze`/`upgrade-grouping` commands keep their require-existing-DB guard.
- [x] **9b. P2 — onboarding welcome tour copy STALE** (`418b7e10`) — the 4-step `OnboardingFlow` tour described the retired *"Triage, Index, Series stage tabs run across the top"* model. Rewrote all 4 slides to the current flow (directory → grouping review → screen the corpus → Focus the trace → Series). No em-dashes.
- [x] **9c. P2 — `Input` not a `forwardRef`** (`d3583916`) — `Popover` clones its trigger (`Input`) and attaches a ref for Escape focus-return; a plain function component dropped it ("Function components cannot be given refs" on every funnel visit, + broken focus-return). `Input` now forwards its ref to the inner `<input>`, merged with the existing `inputRef`. Live-verified: funnel console 1 error → 0.
- [ ] **9d. P3 — picker placeholder slightly over-promises.** "Start typing and we suggest matches" implies suggestions appear on type, but `Popover` is self-managed (no `open` prop) and only opens on a trigger *click*. In normal use (click-to-focus, then type) the autocomplete works — verified live, the suggestion popover shows `mini-1p7m`. Only a programmatic fill without a focus-click misses it. Minor: either reword the hint or give `Popover` a controlled `open` prop so DirectoryPickerField can auto-open on `hasSuggestions`. *Pre-existing; not introduced by 8c/9c.*

## Done

_Move checked items here with their commit sha as they land._

(none yet)
