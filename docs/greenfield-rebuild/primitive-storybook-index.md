# Primitive Storybook index — Phase-1 review map

The full `src/print/ui/` primitive set, browsable in Storybook (`npm run storybook` →
http://localhost:6006, sidebar group **ui/**). Every entry renders on the warm-paper
field with live Print tokens. This is the reviewer's map for the Phase-1 primitive sweep
(plan: `docs/superpowers/plans/2026-05-31-greenfield-ui-primitive-sweep.md`).

Each primitive: TDD'd against `data-*`/roles (never class strings), passes the
`docs/greenfield-rebuild/impeccable-primitive-checklist.md` bar (A–J), and clears
`npx impeccable detect` as a set (exit 0). 24 primitives carry their own `test/print-ui/`
spec (98 tests); the 9 seeded primitives carry their pre-existing `test/ui/` specs.

## Hardened existing primitives (mockup gaps closed)

| Storybook title | Story states | From mockup |
|---|---|---|
| `ui/Button` | solid · accent · ghost · danger · **Armed** | focus-plot "+ Peak" armed tool button |
| `ui/Card` | Default · Elevated · **Draft** | folio recipe **draft** card (dashed, no Plate-Lift) |
| `ui/SegmentedControl` | Small · Medium · **MiniXs** | focus-plot mini comb/resid toggle (`xs`) |
| `ui/PhaseChip` | Default · Solid · **Coexistence** | builder coexistence chip ("Pn3m + Lam", dominant color) |

## New atomic / layout / form primitives

| Storybook title | Story states | From mockup |
|---|---|---|
| `ui/Badge` | Default · WithinButton | notes-button inline mono count |
| `ui/TopBar` | Default | fixed 56px header shell (banner landmark) |
| `ui/StageTabs` | Samples · Index · Series | workspace stage tabs (leading accent Dot) |
| `ui/FacetChip` | Default | facet dropdown-trigger pill (`aria-haspopup`) |
| `ui/KbKey` | Default · Letter | keyboard key cap (`<kbd>`, 3D bottom border) |
| `ui/KbLegend` | Default | shortcut legend row (composes KbKey) |
| `ui/TagPill` | Default · Removable | single tag pill (× → accent) |
| `ui/TagList` | WithTags · Empty · Editable | wrapping tag list + reveal-on-hover add invite |
| `ui/EmptyState` | FilterNoMatch · TitleOnly | centered teaching empty message (serif title) |
| `ui/ScreenedMark` | Screened · Unscreened | completion badge (fill + check + aria-label) |
| `ui/RejectOverlay` | OverThumbnail · OverFrame | grease-pencil terracotta ✕ over a rejected frame |
| `ui/ProgressBar` | Partial · Full · Empty | accent capsule progress track (role=progressbar) |
| `ui/BonnetBadge` | Default · InCandidateRow | Gauss–Bonnet flag pill (color-mix tint, word) |
| `ui/GripHandle` | Default · InRow | drag handle (⋮⋮, row-hover reveal) |
| `ui/SearchInput` | Empty · WithValue | search field (focus-within accent ring) |
| `ui/FilterChip` | Off · On | filter toggle pill (ink/paper inversion) |
| `ui/Slider` | Default · WithValue | offset range input (`.print-range` thumb) |
| `ui/ToggleSwitch` | Off · On | pill switch (sliding knob, accent ON track) |
| `ui/MetaList` | Default | mono key/value list (semantic `<dl>`) |
| `ui/SignalBars` | Strong · Weak · Full | 5-bar strength indicator (count + aria-label) |

## Seeded primitives (carried from the closed-look library; stories added for review)

| Storybook title | Story states |
|---|---|
| `ui/Dot` | per tone (accent/success/muted/neutral) + Labeled |
| `ui/HintText` | Default |
| `ui/IconButton` | per tone (ghost/accent/danger) + WithGlyph + Disabled |
| `ui/Kicker` | Default + per tone |
| `ui/ModalShell` | Default · Small · Large · TopAligned · Drawer (rendered open) |
| `ui/PeakGlyph` | TriangleAuto · Diamond · Caret · Hot · Excluded |
| `ui/PhaseStrip` | Default · Throughout · Coexistence · Unindexed · Small · Vertical |
| `ui/ScoreBar` | High · Mid · Low + Compact |
| `ui/Toast` | Info · Success · Warning · Error (seeded via `showToast`) |

## Notes for the reviewer

- **Second-Channel (color-blind) primitives to spot-check in grayscale:** `ScreenedMark`
  (fill + checkmark), `SignalBars` (bar count), `ProgressBar` (fill width), `PhaseChip`
  (phase name text), `RejectOverlay` (✕ shape + consumer frame dimming), `FilterChip`
  (ink/paper inversion). None rest meaning on hue alone.
- **Accent is rationed** to interaction marks only (focus rings, the armed/ON state, the
  reject ✕, the progress fill, the live slider thumb, the tag-remove/add hover) — never on
  resting chrome.
- **Touch-target exceptions** (documented in-code): `SegmentedControl` `xs`, `StageTabs`,
  `FacetChip`, `FilterChip`, `ToggleSwitch` — dense in-row/in-bar controls whose effective
  hit area is enlarged by their container row.
- **Phase-2 consumer responsibility:** `GripHandle` is `aria-hidden`; the consuming row
  must provide a keyboard reorder path.
- **Motion:** transitions are color/opacity/transform only at 120ms; a global
  `prefers-reduced-motion` reset in `styles.css` near-zeros all motion.

## Deferred to Phase 2 (not in this sweep)

`Sparkline`, `MiniWaterfall` (figure renderers needing real trace data) and all composites
(`DataTable`, `FolioCard`, `CombPanel`, `PhasecallBlock`, `CandidateRow`, `ScopingPlate`,
`NotesMargin`, `ExposureSwitcher`, `Stepper`, …) that assemble these primitives.
