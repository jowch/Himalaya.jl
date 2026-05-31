# Design-system catalog (extracted from The Print mockups)

Source mockups: `docs/redesign-mockups/*.html` + `DESIGN.md`
Baseline primitives: `src/components/ui/**` (to be seeded into `src/print/ui/` in a later task)

---

## Existing primitives — confirmed in mockups

All 13 `.tsx` files found in `src/components/ui/`:
`Button.tsx`, `Card.tsx`, `Dot.tsx`, `HintText.tsx`, `IconButton.tsx`, `Kicker.tsx`,
`ModalShell.tsx`, `PeakGlyph.tsx`, `PhaseChip.tsx`, `PhaseStrip.tsx`, `ScoreBar.tsx`,
`SegmentedControl.tsx`, `Toast.tsx`

| Primitive | Appears in | Variants/states seen | Gap vs current |
|---|---|---|---|
| **Button** | All 7 mockups | ghost (most chrome), solid (ink fill — "Confirm & build", "Export figure"), accent (terracotta — "+ New series", "Drop X" cull bar), armed/active state (tool-btn.armed — "+ Peak" armed), disabled (`.build-btn:disabled` in scoping) | Current: `solid`, `accent`, `ghost`, `danger`. **Missing: armed/active toggle state** (tool-btn.armed uses accent-fill; there is no `armed` prop on Button), and **disabled-via-data-state** styling for the scoping gate button. |
| **Card** | All 7 mockups | elevated plate (trace plate, loupe-plate, scope-plate, series plate), flat panel (detector panel, reflections panel), draft card style (`.is-draft` dashed border + no shadow) in folio | `elevated` + `border` props exist. **Missing: `draft` variant** (dashed border, reduced opacity) seen on folio draft cards. |
| **SegmentedControl** | All 7 mockups | bordered (log/linear, waterfall/heatmap, comb/indexing-space, contact-sheet/loupe, sort), plain (stage tabs — Samples/Index/Series in topbar, which uses a different visual pattern) | `bordered` + `plain` variants exist. **Missing: `size="xs"` step** — mini comb-panel toggle (focus-plot: `.seg.mini button` at `4px 9px` / `10px`), smaller than the current `sm`. |
| **PhaseChip** | sample-table (status cell), focus-workspace (rail candidate name, reflections table), series-builder (trace list, right-margin of waterfall), focus-plot (assignment cart, candidates, custom-index preview), series-plot (phase reading panel) | tint (default), solid (series waterfall right margin), coexistence label ("Pn3m + Lam") | Current: `tint` + `solid`. The coexistence two-phase chip ("Pn3m + Lam") is rendered by the builder as a plain text string — **no `coexist` variant** in the primitive. |
| **PhaseStrip** | series-folio (card body), series-scoping (preview), series-builder (implicit via card rendering), series-plot (companion strip beside waterfall SVG) | horizontal (folio, scoping), vertical (series-plot companion), `form_factor` / `null` 3-state assignment | `orientation` prop + `form_factor`/`null` states already landed (Plan E). No gap beyond confirming all states appear in mockups. |
| **ModalShell** | focus-plot (custom-index builder `.scrim`/`.custom-sheet`) | dialog (centered, 600px max-width), scrim click to close, Esc close, kicker + title + body + footer layout inside shell | Current: `dialog` + `drawer` variants, `sm/md/lg` sizes. The custom-index sheet is `max-width: 600px` — `size="md"` (640px cap) covers it well. No gap. |
| **Kicker** | All 7 mockups | accent (terracotta section eyebrow — "Contact sheet", "Integration", "Folio", "New series", "Speculative"), faint (panel labels — "Phase call", "Reflections", "Ordered by", "Compose") | `accent` + `faint` tones exist. No gap. |
| **IconButton** | focus-workspace (topbar Notes badge counts — not strictly icon-only), series-builder (railtog collapse ‹/›, rail-back restore tab), focus-plot (custom-index close ×) | ghost (default), close × (dismiss) | `ghost`/`accent`/`danger` tones + `dismiss` prop exist. **The topbar "Notes" button with a monospace badge count** is not an icon-only button — it is a ghost Button with a `.badge` span appended; that badge pattern has no primitive. |
| **ScoreBar** | focus-workspace (phase-call `.pcb-bar`, candidate `.cs-bar`), focus-plot (assignment cart `.asb-bar`, candidate `.cs-bar`) | `bar` (full-width, 4px, in phase-call block), `compact` (46px × 3.5px, in candidate rows) | `bar` + `compact` sizes exist. No gap. |
| **Dot** | sample-table (`.stage-tab .dot` active indicator — 5px, terracotta accent when active), focus-workspace (same topbar dot), series-folio (same), series-scoping (`.fs-dot` confirmation state — 8px, sage/accent), series-plot (legend dots) | neutral (inactive tab dot), accent (active tab dot), success (ready state), muted (open circle outline) | `neutral`/`accent`/`success`/`muted` tones, `xs`/`sm` sizes. The tab dot is 5px (`xs`). The scoping confirmation dot is 8px — slightly larger than `sm` (8px = `sm`). No gap — `sm` = `h-2 w-2` = 8px covers it. |
| **HintText** | focus-workspace (rail note), series-builder (rail hint text), series-scoping (`.of-note`, `.foot-note`) | single italic faint style | Single appearance. No gap. |
| **PeakGlyph** | focus-plot (trace overlay: downward triangle auto ▽, diamond manual ◆, caret predicted-absent; legend swatches), series-plot (track annotation anchors, minimal dot anchors) | triangle-down (auto, filled), diamond (manual, filled), caret (predicted-absent, hollow), excluded (hollow/faint of either) | Full glyph vocabulary already implemented (`triangle-down`, `diamond`, `caret`, `ring`, opacity). No gap. |
| **Toast** | Not explicitly shown in any mockup surface — implied by the mutation queue's confirmation feedback | info / success / warning / error | No mockup surface directly shows toasts; all surfaces are static. No changes needed from mockup inspection. |

---

## Gap primitives — recurring patterns with NO current primitive

| Proposed primitive | Description | Mockups | Prop sketch |
|---|---|---|---|
| **TopBar** | Fixed 56px header shell: wordmark + stage tabs + spacer + right slot. Appears identically in all surfaces. Currently assembled inline on every page with duplicated markup. | All 7 | `rightSlot?: ReactNode; activeStage?: "samples" \| "index" \| "series"` |
| **StageTabs** | The three uppercase label tabs (Samples / Index / Series) with 5px dot indicator. Each tab carries a `Dot` whose tone flips accent on active. Currently inlined in every topbar. | All 7 | `active: StageKey; onChange: (s: StageKey) => void` |
| **FacetChip** | Pill-shaped corpus filter button (beamtime facet selector). Pill radius (999px), `plate` bg + `hair-strong` border, `ink`/`ink-faint` two-level label. Appears only on the contact-sheet topbar. | sample-table | `label: string; value: string; onClick: () => void` |
| **DataTable** | The contact-sheet's column-grid table (sample rows, exposure strip, kept count, tags, status). 5-column grid, `paper-sunk` header, hairline row separators, `plate` background, hover highlight. Rows are clickable. | sample-table | `columns: ColDef[]; rows: RowData[]; onRowClick?: (id) => void` |
| **SampleRow** | One row in the contact-sheet DataTable: specimen cell (screened-mark + name + mono id), exposure-thumbnail strip, kept-count cell, tags cell, status cell. | sample-table | `sample: SampleData; selected?: boolean; onThumbClick?, onNameClick?` |
| **ScreenedMark** | 13px circular completion badge (empty outline = unscreened, ink-filled checkmark = screened). | sample-table | `done: boolean` |
| **ThumbnailStrip** | Horizontal row of 62px detector thumbnails, each with rejection ×-mark overlay and representative-dot badge; shift-click range selection. | sample-table, focus-workspace (exposure switcher strip), sample-table loupe | `thumbs: ThumbData[]; selected?: Set<string>; onSelect?, onDoubleClick?` |
| **Thumbnail** | Single dark-window detector-image tile (SVG-backed or `<img>`), with frame-edge border, rejected opacity + terracotta ×-mark overlay, mono frame-number label, representative accent dot. | sample-table, focus-workspace, sample-table loupe | `src: string; rejected?: boolean; isRep?: boolean; selected?: boolean; frameNo?: number` |
| **RejectOverlay** | The grease-pencil ×-mark (two hand-skewed terracotta lines in SVG) rendered over rejected detector frames. | sample-table (thumb and big-frame), sample-table loupe (big-frame), focus-workspace (not shown but designed) | `size?: "sm" \| "lg"` (62px thumb vs full-frame) |
| **ProgressBar** | Accent-filled capsule progress track (screened-so-far counter). Appears in the contact-sheet header progress widget. Distinct from ScoreBar in intent and color (always accent, not phase-colored). | sample-table | `value: number; total: number; className?` |
| **MetaList** | Vertical key-value list in monospace; key=ink-faint, value=ink. Used in the loupe sidebar for exposure metadata (frame, integration time, collected, signal). | sample-table (loupe), series-scoping (implied field labels) | `entries: {key: string; value: ReactNode}[]` |
| **SignalBars** | 5-bar vertical signal-strength indicator, inline-flex, bars colored ink-soft when active / hair-strong when inactive. | sample-table (loupe metadata) | `value: number; max?: number` (0–5) |
| **VerdictCard** | Loupe verdict/status card: dot + bold state text + hint text + toggle button, `paper-sunk` fill, `hair-strong` border, 7px radius. | sample-table (loupe) | `kept: boolean; onToggle: () => void` |
| **RepresentativeBox** | Loupe representative-exposure selection card: "Representative for indexing" accent dot row + body copy + "Set as representative" button. | sample-table (loupe) | `isRep: boolean; rejected: boolean; onSet: () => void` |
| **TagList** | Horizontal wrapping list of pill-shaped free-form sample tags (999px radius, `plate` bg, `hair` border, ink-soft text). Includes a `+ tag` add button that becomes visible on row hover (or always-visible invite state when no tags). | sample-table (row and loupe), series-scoping (sample names context) | `tags: string[]; onAdd?: () => void; editable?: boolean` |
| **TagPill** | Single rounded pill tag (10.5px, 600, ink-soft). The unit of TagList. | sample-table, series-folio (`.fig-tag` metadata tags on plate header) | `children: ReactNode; onRemove?: () => void` |
| **FloatingCullBar** | Fixed bottom-center floating action bar that appears on multi-select of exposure thumbnails: count + "Drop / Restore / Clear" buttons on an ink-fill slab. | sample-table | `count: number; onDrop, onRestore, onClear: () => void` |
| **KbLegend** | Horizontal keyboard-shortcut legend row below surfaces: `<kbd>` key badge + description text in ink-faint. | sample-table (sheet footer), sample-table loupe (sidebar) | `shortcuts: {key: string; description: string}[]` |
| **KbKey** | A single `<kbd>`-styled key badge: mono, `plate` bg, `hair-strong` border + 2px bottom border for 3D effect. | sample-table, focus-workspace (implied), series-scoping (`.sh-r` count mono), focus-plot (custom-index foot note) | `children: ReactNode` |
| **PhasecallBlock** | The focus-workspace rail "Phase call" panel: per-phase block with serif phase name in phase color + score value + meta line (lattice · N peaks) + ScoreBar + series ratios. Multi-phase shows a "Coexistence · N phases" tag at top. | focus-workspace, focus-plot (assignment cart as `.assign`) | `phases: PhaseCall[]; empty?: ReactNode` |
| **CandidateRow** | Clickable candidate-indexing row: checkbox mark + phase name (mono/serif) + covers N peaks + Bonnet badge (focus-plot) + score + mini ScoreBar. Hover previews peaks in phase color; checked adds a terracotta inset ring. | focus-workspace, focus-plot | `candidate: CandidateData; checked: boolean; onToggle, onMouseEnter, onMouseLeave` |
| **BonnetBadge** | Small pill badge: "BONNET" text, accent color, accent-tinted background, 999px radius. Appears inline in the candidate row when the Gauss–Bonnet ratio flag is set. | focus-plot only | `(no props — purely decorative badge)` |
| **Stepper** | Previous/next navigation control with a centered label slot showing sample name + "sample N of M" sub-label. Used in the focus workspace and focus-plot topbar. | focus-workspace, focus-plot | `label: string; sublabel?: string; onPrev, onNext: () => void; prevDisabled?, nextDisabled?` |
| **SampleStepper** | Alias/specialization of Stepper for "sample N of M" navigation; may just be a configured `<Stepper>`. | focus-workspace, focus-plot | see Stepper |
| **ExposureSwitcher** | Row of small (30–32px) exposure thumbnail buttons in the detector panel header; selected frame has terracotta 2px ring; rep frame has small accent dot badge. | focus-workspace, focus-plot | `exposures: ExposureData[]; selected: number; onSelect: (n) => void` |
| **NotesMargin** | Third column aside: ruled quiet column, not a card. "Notes" heading + mono count + note entries (with author/age meta) + "Add a note" dashed-bottom prompt. Collapses to a topbar button below 1320px. | focus-workspace | `notes: NoteData[]; onAdd: () => void; collapsed?: boolean` |
| **NoteEntry** | Single timestamped note in the notes margin: small uppercase meta (author, age, accent author color) + body copy with inline `<code>`-style q-ref spans. | focus-workspace | `author: string; age: string; children: ReactNode` |
| **SearchInput** | Full-width search bar: magnifier icon + text input, `plate` bg, `hair-strong` border, 7px radius, ink-faint placeholder. | series-folio | `value: string; onChange: (v: string) => void; placeholder?` |
| **FilterChip** | Toggle pill button for filtering (not the facet chip — this one sits in a `.chips` row with `on` state that fills ink-on-paper). Series folio "All / Has transition / Cross-experiment". | series-folio | `label: string; active: boolean; onClick: () => void` |
| **FolioCard** | The series folio card: mini-waterfall SVG figure + card body (fig-n kicker, title, meta, PhaseStrip, card foot with provenance + edit time). Hover lifts 3px. Draft state: dashed border. | series-folio | `series: SeriesData; onClick: () => void` |
| **MiniWaterfall** | Small read-only SVG stacked waterfall (the "frozen plate" figure on folio cards). N traces, each a colored fill + line + peak ticks; baseline dashes. | series-folio, series-scoping (spark preview for candidates) | `traces: MiniTraceData[]; width?: number` |
| **Sparkline** | Tiny 76×28px single-trace SVG trace preview. Colored by dominant phase. Used in the scoping surface per-sample row and candidate row. | series-scoping | `comps: PhaseComp[]; width?: number; height?: number` |
| **ScopingPlate** | Centered worksheet plate (max 760px): kicker + editable title + autogroup summary + ordering field + samples list + candidates section + phase strip preview + confirmation footer. | series-scoping | The whole plate; internal components are AutogroupBanner, OrderingField, ScopingSampleRow, etc. |
| **AutogroupBanner** | Machine-proposal info card: terracotta star icon + body text explaining the grouping. `paper-sunk` bg, `hair` border, 7px radius. | series-scoping, series-builder rail | `children: ReactNode` |
| **OrderingField** | Re-openable field showing the ordering variable value with a chevron: `plate` bg, `hair-strong` border, 7px radius, cursor pointer. | series-scoping, series-builder rail | `value: string; onClick: () => void` |
| **ScopingSampleRow** | One row in the scoping series: grip handle + sparkline + sample name + mono ID + parsed value cell (flagged/ok states with dashed-underline affordance). | series-scoping | `sample: ScopingSample; onToggleFlag: (id) => void` |
| **GripHandle** | Vertical drag handle (two columns of dots, `⋮⋮`). Shown on hover in sortable lists. | series-scoping (sample rows), series-builder (trace list rows) | `(no props — decorative)` |
| **RailSlidersSection** | Series-builder rail section containing range slider + value display + label. Offset control + scale toggle + track toggle. | series-builder | Could just be a composed use of a `Slider` primitive. |
| **Slider** | Range input styled with accent thumb, `hair-strong` track, 3px height. | series-builder | `min, max, step, value: number; onChange: (v: number) => void; label?: string` |
| **ToggleSwitch** | Small 32×18px pill toggle (on=accent fill, off=hair-strong; white disc thumb). "Track reflections" annotation toggle in series-builder rail. | series-builder | `checked: boolean; onChange: (v: boolean) => void; label: string` |
| **TraceListRow** | Series-builder rail trace list row: grip + phase-colored dot (9×9, 2.5px radius) + sample name + dose/id + phase chip. Draggable to reorder. | series-builder | `trace: TraceData; onDragStart?` |
| **FloatingDock** | Fixed bottom-right floating control dock (shown in full-bleed mode only): label + range slider + value display. Single-control dock in the series-builder. | series-builder | `label: string; children: ReactNode` |
| **RailbackTab** | Fixed right-edge tab to restore the collapsed rail: left-rounded plate, arrow + rotated vertical "Compose" label. | series-builder | `label: string; onClick: () => void` |
| **CombPanel** | The focus-plot "Reflections — comb" panel: q-ruler with teeth per assigned phase (observed=solid stem+cap; absent=dashed+caret). Includes comb/indexing-space toggle and legend row. | focus-plot | `assignment: AssignedPhase[]; peaks: ObservedPeak[]; view: "comb" \| "resid"` |
| **CombLegend** | Three-item icon+label legend for the comb panel: predicted & observed / predicted absent / leftover peak. | focus-plot | `(no props — static, uses SVG swatches inline)` |
| **CustomIndexSheet** | Modal content: symmetry segmented control + lattice slider + number input + preview comb SVG + fit summary + Cancel/Add buttons. | focus-plot | `onAdd: (candidate) => void; onClose: () => void` |
| **RegionBracket** | Left-margin vertical bracket on the series waterfall: a colored line + rotated region label ("Pn3m", "coexistence", "Lamellar"). | series-builder, series-plot | `regions: Region[]; traces: TraceData[]` |
| **PhaseReadingPanel** | Series-plot rail "Phases present" card: per-phase row with dot + name + span (sample range) + lattice line + coexistence summary. | series-plot | `phases: PhaseReading[]` |
| **MemberList** | Series-plot rail trace member list: per-trace row with swatch dot + sample name (+N coexisting) + var value + data line. Hover highlights trace. | series-plot | `members: MemberData[]; hotIndex?: number; onHover?: (i) => void` |
| **ExportButton** | Wide plate-fill button with icon + title + subtitle layout. "Export as CSV / SVG" in series-plot rail. | series-plot | `icon?: ReactNode; title: string; subtitle?: string; onClick: () => void` |
| **EmptyState** | Centered layout: serif h2 + body copy. Shown on the folio when search/filter matches nothing; on contact-sheet when no samples; zero-exposure state. | series-folio (filter empty), implied by functional-redesign for other surfaces | `title: string; body?: ReactNode` |

---

## Tokens

Tokens defined in `src/styles.css` `@theme` that are confirmed used in the mockups: all paper/ink/hair/accent/phase/status colors, radius scale (5px), type scale (xs through display), font families (sans/serif/mono).

| Token role | Status | Note |
|---|---|---|
| `--color-unindexed` | **MISSING from @theme** | Used in focus-workspace + focus-plot + series-plot for gray "unindexed" rings/peaks/strip cells: `oklch(0.660 0.012 76)`. Not in `styles.css @theme`; referenced only by inline CSS vars in mockups. |
| `--color-manual` | **MISSING from @theme** | focus-plot defines `--manual: oklch(0.550 0.200 340)` for manual-peak magenta. Exists in `DESIGN.md` as `peak-manual` but not in `styles.css`. |
| `--text-display` step 2 sizes | Confirmed | `--text-display: 27px`, `--text-headline-lg: 26px` both in @theme. Folio uses 31–33px serif; builder uses 33px — these exceed `--text-display` (27px). **A `--text-display-lg` step (~30–33px) is absent.** |
| `--color-frame-tag` | Confirmed present | `oklch(0.82 0.01 80)` — already in @theme. |
| `--color-scrim` | Confirmed present | `oklch(0.05 0 0 / 0.65)` in @theme (mockups use 0.5–0.55). Minor discrepancy but close. |
| Phase color `--color-phase-*` | **MISSING from @theme** | Phase colors live in `src/phases.ts` (`phaseColor()`), not as `@theme` tokens. Mockups define them as `:root` CSS vars. For the greenfield `src/print/` layer, all 8 phases should be `@theme` tokens (`--color-phase-pn3m`, etc.) so Tailwind generates utilities (`text-phase-pn3m`, `bg-phase-pn3m`). |
| `--text-kicker` role class | Confirmed — `.text-kicker`, `.text-kicker-accent`, `.text-kicker-faint` in styles.css | No gap. |
