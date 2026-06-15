# Greenfield Phase-4 Cutover Ledger

> Generated 2026-06-03 by a survey workflow. The evidence base for the Phase-4 cutover strategy brainstorm. Branch worktree-greenfield-ui-rebuild (unmerged).

## How to read this

This is a JOIN, not a list — each row joins a legacy capability to its greenfield replacement, tagged **carry** / **drop** / **new** / **refactor** / **gap**. A *carry* survives the cutover untouched; a *drop* is dead code to delete; a *new* is a greenfield-only piece; a *refactor* changes shape across the boundary; a *gap* is a behaviour that must survive but has no greenfield piece yet. The load-bearing rows are the **GAPS** (carry with no greenfield piece) and the **DEAD CODE** (drop) — those two sets drive the cutover plan. Where the surveys disagree (most notably whether `comparisons.jl` is dead or shared) the conflict is flagged inline rather than silently resolved.

---

## Frontend — per route surface

### Samples / contact-sheet (`/samples`)

`SamplesPage` renders the corpus contact-sheet: a screened-progress header, a 5-column grid of rows (specimen + thumbnail strip + kept count + tags + status door), a floating CullBar, and a sticky kb-legend, wrapped by the `CorpusShell` layout route + `CorpusTopbar`. The entire panel/renderer layer (`SheetTable`, `SampleTableRow`, `CullBar`, `Thumbnail/ThumbnailGallery`, cell composites, `ProgressBar`) is built in `src/print/`. The cutover is a greenfield Layer-4 page that lifts the same data flow (`useCorpusSamples` + `useCorpusExposures` + `useScreenedProgress` + three queue mutations) into a new shell composing the built Print panels.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| SamplesPage | component | refactor | pages/SamplesPage.tsx | — | Page shell not in print/pages/; owns beamtime URL, filter, screened-progress, Skeleton, kb-legend |
| CorpusShell | component | refactor | components/CorpusShell.tsx | print/ui/TopBar.tsx | Layout route wrapping /samples, /loupe, /sample; must be a Layer-4 shell |
| CorpusTopbar | component | refactor | components/CorpusTopbar.tsx | print/ui/TopBar.tsx | TopBar+StageTabs+Stepper built; router-aware wiring is the gap |
| ContactSheetRow | component | drop | components/ContactSheetRow.tsx | print/components/SampleTableRow.tsx | Stateful row → presentational; select/keyboard/mutation move to page |
| ExposureThumb (inline) | component | drop | ContactSheetRow.tsx:ExposureThumb | print/components/Thumbnail.tsx | rep dot/reject/frame-badge covered; dbl-click→loupe + shift-click are page gaps |
| CullBar (legacy) | component | drop | components/CullBar.tsx | print/components/CullBar.tsx | Greenfield adds Restore; aria-hidden/tabIndex when hidden |
| SampleStatusChip | component | drop | components/SampleStatusChip.tsx | print/components/StatusCell.tsx | Greenfield StatusCell has no door Link — page must wrap |
| RejectXMark | component | drop | components/RejectXMark.tsx | print/ui/RejectOverlay.tsx | grease-pencil mark primitive |
| DetectorImage (legacy) | component | carry | components/DetectorImage.tsx | print/detector/DetectorImage.tsx | Wire only the print/ version post-cutover |
| Kicker / HintText (legacy ui) | component | drop | components/ui/* | print/ui/* | tone variants exist |
| SheetTable | component | new | — | print/components/SheetTable.tsx | Children-slotting table container, SAMPLE_TABLE_COLS shared |
| SampleTableRow | component | new | — | print/components/SampleTableRow.tsx | Presentational 5-col; selectedExposureIds Set seam |
| SpecCell/KeptCell/StatusCell | component | new | — | print/components/SpecCell.tsx | StatusCell lacks the focus door Link |
| ThumbnailGallery + Thumbnail | component | new | — | print/components/ThumbnailGallery.tsx | onSelect is (id)=>void; no shift/dbl-click/loupe |
| ProgressBar | component | new | — | print/ui/ProgressBar.tsx | replaces inline width div |
| useCorpusSamples / useCorpusExposures / useScreenedProgress / useExperiments | hook | carry | queries.ts | — | useCorpusExposures O(1)-observer bulk hook is load-bearing |
| useSetExposureStatus / useAddCorpusSampleTag / useRemoveCorpusSampleTag | hook | carry | queries.ts | — | Carry hook; call site moves to page |
| GET /api/samples · /:id/exposures · /experiments | endpoint | carry | api.ts | — | No backend change |
| GET /api/events (SSE) | sse | carry | App.tsx | — | Global; contact sheet heals via cache invalidation |
| Mutation queue | queue | carry | lib/queue/ | — | Global; page calls same hooks |
| isSampleScreened / sampleDisplayName | hook | carry | lib/sample/* | — | Carry into page data mapping |
| boneyard contact-sheet skeleton | capability | carry | bones/contact-sheet.bones.json | — | May need recapture for new SheetTable layout |
| Shift-click range-select | capability | gap | ContactSheetRow:handleSelect | — | Needs onSelectExtend or page-level range logic |
| Dbl-click / Enter/Space → loupe | capability | gap | ContactSheetRow:openLoupe | — | No onOpenLoupe on ThumbnailGallery |
| X / Esc cull shortcuts | capability | gap | ContactSheetRow:useEffect keydown | — | Re-mount listener at page against page-owned selection |
| Restore in CullBar | capability | new | — | print/components/CullBar.tsx | New batch-restore loop |
| Tag inline form | capability | gap | ContactSheetRow:tagFormOpen | print/ui/TagList.tsx | Mutation wiring same; form UX moves into TagList |
| corpus→focus door (StatusCell) | capability | gap | ContactSheetRow:Link | print/components/StatusCell.tsx | M-1 door — page must add Link/href |
| Beamtime ?beamtime= filter | capability | gap | SamplesPage:beamtime | — | Read at page, write at topbar |
| View-seg switch | capability | gap | CorpusTopbar:view-seg | — | Topbar Layer-4 wiring |
| SeriesCommitConflictModal (in App) | component | drop | components/SeriesCommitConflictModal.tsx | — | **CONFLICT: this survey marks it drop ("409 flow cancelled"); app-shell & data-plane surveys mark it carry (series_commit conflict still live). See Risks.** |

### Loupe (`/samples/loupe/:sampleId`)

Single-sample inspection: full detector image + filmstrip left (`LoupeFrame`), metadata/verdict/representative/tags sidebar right (`LoupeSidebar`). Queries samples, exposures, peaks; two queue mutations (status, select) + two tag mutations. Greenfield `BigFrame` + `LoupeSidePanel` + sub-components are built; the Layer-4 page is the gap, plus the Exposure→GalleryExposure src adapter, the "Open in Index stage" link, defaultExposureId logic, keyboard handler, and beamtime back-nav.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| LoupePage | route | gap | pages/LoupePage.tsx | — | Biggest gap; owns activeId, defaultExposureId, keyboard, back-nav, Skeleton, chrome |
| LoupeFrame | component | drop | components/LoupeFrame.tsx | print/components/BigFrame.tsx | Page computes caption + maps to GalleryExposure |
| LoupeSidebar | component | drop | components/LoupeSidebar.tsx | print/components/LoupeSidePanel.tsx | Presentational; Index-stage link missing |
| DetectorImage (legacy) | component | drop | components/DetectorImage.tsx | print/detector/DetectorImage.tsx | Greenfield takes raw src; page builds URL |
| ThumbnailGallery (legacy) | component | drop | components/ThumbnailGallery.tsx | print/components/ThumbnailGallery.tsx | Map Exposure[]→GalleryExposure[] |
| RejectXMark (legacy) | component | drop | components/RejectXMark.tsx | print/ui/RejectOverlay.tsx | Already wired inside BigFrame |
| LoupeTagsEditor (inline) | component | drop | LoupeSidebar:LoupeTagsEditor | print/ui/TagList.tsx | State encapsulated in TagList |
| useCorpusSamples / useExposures / usePeaks | hook | carry | queries.ts | — | usePeaks drives signal level |
| useSetExposureStatus / useSelectExposure / useAddCorpusSampleTag / useRemoveCorpusSampleTag | hook | carry | queries.ts | — | Tag.id adapter: onRemove=(t)=>mutate(t.id) |
| GET /api/samples · /:id/exposures · /:id/image · /:id/peaks; PATCH status · select; POST/DELETE tags | endpoint | carry | api.ts | — | No backend change |
| SSE (via queue) | sse | gap | — | — | No direct SSE; inherits queue replay |
| useQueueMutation | queue | carry | lib/queue/ | — | All four mutations queue-backed |
| boneyard Skeleton | component | gap | LoupePage.tsx | — | No loupe.bones.json yet |
| defaultExposureId logic | capability | new | LoupePage.tsx:20-26 | — | rep→first-accepted→first; reproduce in page |
| Open in Index stage link | capability | gap | LoupeSidebar.tsx:328 | — | RepresentativeBox has no Link; inject via children/sibling |
| beamtime back-nav | capability | new | LoupePage.tsx:137 | — | preserve ?beamtime= on back |
| Loupe keyboard handler | capability | new | LoupePage.tsx:146 | — | ←/→ flip, X drop, R rep, Esc back |
| SignalMeter / level derivation | capability | new | LoupeSidebar.tsx:43 | print/ui/SignalBars.tsx | Pass SignalBars node into MetaList entry |

### Focus workspace (`/sample/:sampleId`)

Five panels around Zustand `activeSampleId`/`activeExposureId`/`hoveredQ`: trace hero, detector panel, comb, assignment rail, notes margin. Cutover: swap each legacy panel for its built greenfield counterpart (`TracePlate`, `DetectorPanel`, `CombsPanel`, `AssignmentRail`+`AssignmentCart`+`CandidateList`, `CustomIndexModal`), write a Layer-4 page owning shared q-hover state. Three behaviours live only in legacy with no greenfield piece: `StaleIndicesBanner`, `SpeculativeBuilder`, `FocusNotesMargin`.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| FocusWorkspacePage | route | refactor | pages/FocusWorkspacePage.tsx | — | New print/pages/FocusPage.tsx owning scale/xDomain/hoveredQ/combView |
| FocusWorkspaceLayout | component | drop | components/FocusWorkspaceLayout.tsx | — | Grid + notes drawer collapse into page |
| useSyncActiveSampleFromRoute | hook | carry | hooks/useSyncActiveSampleFromRoute.ts | same | URL→Zustand shim |
| PlotCard | component | drop | components/PlotCard.tsx | print/components/TracePlate.tsx | auto-fit/armed/export lift to page; build TraceModel + PlotPeak[] |
| TraceViewer | component | drop | components/TraceViewer.tsx | print/plot/TracePlot.tsx | Observable Plot → d3; losingPeakIds dim needs inversion |
| FocusPlotHeader | component | refactor | components/FocusPlotHeader.tsx | print/components/PlateHeader.tsx | logic → page builds props |
| FocusDetectorPanel | component | drop | components/FocusDetectorPanel.tsx | print/components/DetectorPanel.tsx | ring derivation + switcher → page; switcher = ThumbnailGallery in tools slot |
| CombPanel | component | drop | components/CombPanel.tsx | print/components/CombsPanel.tsx | ghost-preview row injected as flagged CombSeries |
| IndicesCard / PhasePanel | component | drop | components/PhasePanel.tsx | print/components/AssignmentRail.tsx | Decomposes to AssignmentCart+CandidateList+CustomIndexModal |
| AssignmentCart / CandidateList/Row | component | new | — | print/components/* | Built Batch 6/8 |
| StaleIndicesBanner | component | gap | components/StaleIndicesBanner.tsx | — | hash compare + debounce + pending-ops suppression; no greenfield port |
| SpeculativeBuilder | component | gap | components/SpeculativeBuilder.tsx | — | anchor+ratio dialog; useSpeculativeSnap + queue gate; no port |
| FocusNotesMargin | component | gap | components/FocusNotesMargin.tsx | — | q-ref highlight + focus-gated textarea + drawer; cache from useSamples(expId) |
| CustomIndexModal (legacy) | component | drop | components/CustomIndexModal.tsx | print/components/CustomIndexModal.tsx | Greenfield presentational; page owns sym/val/snap |
| FigureExportControls | component | drop | components/FigureExportControls.tsx | print/components/ExportButton.tsx | + useFigureExport; traceAdapter carries |
| Stepper | component | new | — | print/components/Stepper.tsx | prev/next sample nav; page wires router |
| useTrace/usePeaks/useIndices/useAssignment/useExposures | hook | carry | queries.ts | same | Shared cache keys |
| useCorpusSamples/useSamples/useExperiment | hook | carry | queries.ts | same | notes MUST read useSamples(expId) |
| useAddPeak/useRemovePeak/useSetPeakExcluded | hook | carry | queries.ts | same | TracePlate interaction handlers |
| useAddAssignmentPhase/useRemoveAssignmentPhase/useDeleteIndex/useCommitCustomIndex/useReanalyzeExposure/useSpeculativeSnap/useCreateSpeculative/useSelectExposure/useUpdateSample | hook | carry | queries.ts | same | All queue-backed |
| useAutoPickExposure / deriveActiveIndices | hook | carry | hooks/ · lib/assignment.ts | same | Move call to page; pure fn |
| trace/peaks/indices/assignment/reanalyze/speculative_snap endpoints | endpoint | carry | api.ts | same | No backend change |
| SSE /api/events | sse | carry | App.tsx | same | Consumer not owner |
| q-link hover (hoveredQ/previewIndexId) | capability | carry | state.ts | same | TracePlot onHoverQ emission deferred (see gap) |
| losingPeakIds dim | capability | gap | PlotCard:losingPeakIds | — | TracePlot has highlightPeakIds (inverse); extend or invert in page |
| notesDrawerOpen / toggleNotesDrawer | capability | carry | state.ts | same | ModalShell drawer at <xl |

### Series folio (`/series`)

Corpus masonry of saved series cards. Fetches the listing once (`useSeriesList`) + per-card `useSeries(id)` for the mini-waterfall/phase strip. Greenfield `SeriesCard`/`Gallery`/`FolioHeader`/`CardFigure` are complete at L3. Cutover replaces `SeriesFolioPage` + `SeriesFolioCard`, bridges API data to `WaterfallRow`, preserves filter/sort/empty-state/SSE. The sharpest issue: `CardFigure` needs real trace data per member, which `useSeries` does not return.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| SeriesFolioPage | route | refactor | pages/SeriesFolioPage.tsx | print/pages (gap) | Owns search/filter/sort/skeleton/empty; replace inline appearance |
| SeriesFolioCard | component | drop | components/SeriesFolioCard.tsx | print/components/SeriesCard.tsx | Presentational; page owns useSeries(id) + toWaterfallRows; stale-dot has no slot |
| SeriesMiniWaterfall | component | drop | components/SeriesMiniWaterfall.tsx | print/waterfall/CardFigure.tsx | Synthetic → real TracePlot; needs trace data |
| buildMiniWaterfall / folioFigure.ts | component | drop | lib/series/folioFigure.ts | print/waterfall/waterfallModel.ts | toWaterfallRows needs tracesById; buildPhaseStrip reusable |
| folioFilter | component | carry | lib/series/folioFilter.ts | same | Pure, unit-tested |
| useSeriesList / useSeries | hook | carry | queries.ts | same | per-card useSeries moves to page scope |
| GET /api/series · /:id | endpoint | carry | api.ts | same | No shape change |
| SSE series invalidation | sse | carry | applyRemoteToCache.ts:267 | same | created/updated/deleted/committed/pinned arms |
| Skeleton / FOLIO_FIXTURE | component | refactor | SeriesFolioPage:FOLIO_FIXTURE | print/pages (gap) | New bone fixture vs print/ui/Card masonry |
| Card / SegmentedControl / Kicker / PhaseStrip (legacy ui) | component | drop | components/ui/* | print/ui/* | Drop-in API parity |
| FolioHeader / Gallery / SeriesCard / CardFigure | component | new | — | print/components/* | FolioHeader has no actions slot |
| SearchInput / FilterChip | component | new | — | print/ui/* | Replace raw inputs/buttons |
| Button accent (+New series) | component | refactor | SeriesFolioPage:94 | print/ui/Button.tsx | Inline appearance is guard violation |
| Empty-state CTAs | component | refactor | SeriesFolioPage:154 | print/ui/EmptyState.tsx | corpus-empty + no-match |

### Series scoping (`/series/new`)

Machine-proposes/human-confirms worksheet that writes (key,value) `sample_tags` with source='scoping' then navigates to the folio. The whole greenfield presentational set (`ScopePlate`, `ScopeSampleRow`, `ScopeCandidateRow`, `Sparkline`, `FlagButton`, `AutoGroup`, `Field`) is built + Storybook-verified. Remaining work: the Layer-4 page + a thin greenfield confirm modal.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| SeriesScopingPage | route | refactor | pages/SeriesScopingPage.tsx:60 | print/pages/SeriesScopingPage.tsx | Port seed/undo/sort/rekey/confirm/nav-gate; presentational shell built |
| ScopingRow | component | drop | components/ScopingRow.tsx:21 | print/components/ScopeSampleRow.tsx | grip restored; FlagButton replaces value cell |
| ScopingValueCell | component | drop | components/ScopingValueCell.tsx:21 | print/ui/FlagButton.tsx | inline edit DROPPED — flag toggle only (intentional) |
| ScopingSparkline | component | drop | components/ScopingSparkline.tsx:16 | print/plot/Sparkline.tsx | self-contained renderer |
| ScopingAutogroupCard | component | drop | components/ScopingAutogroupCard.tsx:12 | print/components/AutoGroup.tsx | grouping prose as ReactNode |
| ScopingOrderField | component | drop | components/ScopingOrderField.tsx:15 | print/ui/Field.tsx | static read-out → real dropdown (upgrade) |
| ScopingLooseMatches | component | drop | components/ScopingLooseMatches.tsx:26 | print/components/ScopeCandidateRow.tsx | show-more must move to page |
| ScopingFoot | component | drop | components/ScopingFoot.tsx:17 | print/components/ScopePlate.tsx | inlined into plate footer |
| ScopingConfirmModal | component | gap | components/ScopingConfirmModal.tsx:17 | — | No greenfield modal; ModalShell ready, ~20 lines |
| ScopePlate / ScopeSampleRow / ScopeCandidateRow | component | new | — | print/components/* | Batch 12 |
| useDragReorder | hook | new | — | print/components/useDragReorder.ts | New manual-override reorder; kb a11y not yet wired |
| useCorpusSampleTags / useCorpusPickerSamples / useScopeSeries / useMemberTraces / useMemberIndices | hook | carry | queries.ts | same | scopeSeriesMutator no-op optimistic |
| GET /api/sample-tags · /picker-samples; POST /samples/tags/batch | endpoint | carry | api.ts | same | N add_tag frames under one client_op_id |
| SSE add_tag fan-out | sse | carry | applyRemoteToCache.ts:309 | same | corpusSampleTags + corpusPickerSamples invalidation |
| proposeOrdering / parseSortKey / dominantPhase / humanizeKey | hook | carry | lib/scoping/* | same | Pure, framework-agnostic |
| lib/plot/sparkline.ts | component | drop | lib/plot/sparkline.ts:15 | print/plot/Sparkline.tsx | retire after wiring (verify no other consumer) |
| Skeleton (boneyard) | component | gap | SeriesScopingPage.tsx:262 | — | No scoping bone yet |
| Empty-state CTA / error banner / pendingBuildRef gate | capability | carry | SeriesScopingPage.tsx | — | ScopePlate has no emptyState/error slot — branch above it |

### Series builder (`/series/:id`)

Combined read/edit surface; draft-gated edit mode via Zustand `seriesDraft` driving `SeriesRecipeEditor`. Greenfield `SeriesPlate`/`WaterfallChart`/`BuilderRail`/`ReadingPanel`/`MemberList`/`ExportButton`/`RailBack`/`Dock` are built. Two hard gaps: (1) no greenfield recipe-editor composite, (2) grouping-mode control + vertical phase-strip companion have no print/ equivalent. Heatmap + cross-trace tracking are intentionally dropped.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| SeriesBuilderPage | route | refactor | pages/SeriesBuilderPage.tsx | print/pages (gap) | Wire draft + 5 hydration hooks + save/commit |
| MultiTracePlot | component | drop | components/MultiTracePlot.tsx | print/components/SeriesPlate.tsx + waterfall/WaterfallChart.tsx | d3 declarative; heatmap out-of-scope; tracking deferred |
| offsetToBandFraction | component | drop | MultiTracePlot:offsetToBandFraction | WaterfallChart offsetScale | semantics equivalent, formula differs — calibrate |
| SeriesBuilderRail | component | drop | components/SeriesBuilderRail.tsx | print/components/BuilderRail.tsx | drops Representation + Track-reflections; recipe editor gap |
| AutogroupCard | component | drop | components/AutogroupCard.tsx | print/components/AutoGroup.tsx (compose) | — |
| RepresentationToggle | component | drop | components/RepresentationToggle.tsx | — | Heatmap out-of-scope; drop cleanly |
| OffsetSlider / ScaleToggle | component | drop | components/* | print/ui/Slider · SegmentedControl | inlined in rail/plate |
| OffsetDock | component | drop | components/OffsetDock.tsx | print/components/Dock.tsx | condition simplifies to 'collapsed' |
| Rail collapse restore | component | drop | SeriesBuilderRail (rail-restore) | print/components/RailBack.tsx | vertical tab |
| SeriesReadingPanel | component | drop | components/SeriesReadingPanel.tsx | print/components/ReadingPanel.tsx | SeriesReading→ReadingRow[] adapter needed |
| SeriesMemberRow | component | drop | components/SeriesMemberRow.tsx | print/components/SeriesMemberRow.tsx | via MemberList |
| MemberMetaGutter | component | drop | components/MemberMetaGutter.tsx | print/components/MemberList.tsx | band-resize + band-aligned layout NOT replicated |
| BandResizeDivider / ActiveBandContext | component | drop | components/* | — | Compare-era; no greenfield, dropped |
| SeriesPhaseStrip | component | gap | components/SeriesPhaseStrip.tsx | print/ui/PhaseStrip.tsx (primitive only) | segmentsFromMembers pure-portable; no composite |
| GroupingModeToggle | component | gap | components/GroupingModeToggle.tsx | — | bySample/byPhase/distinct; decision: keep via ToolBar slot or drop |
| AnnotationToggles | component | gap | components/AnnotationToggles.tsx | — | peak ticks/labels; thread Zustand flags into WaterfallChart rows |
| SeriesRecipeEditor | component | gap | components/SeriesRecipeEditor.tsx | — | **Largest behavioral gap**; title/desc/add-sample/order/Save/Commit wired to queue |
| FigureExportControls | component | drop | components/FigureExportControls.tsx | print/components/ExportButton.tsx + useFigureExport | — |
| buildMultiTraceExportSpec / resolveDisplayLabels / seriesReading / buildSeriesSaveBody/CommitBody | component | carry | lib/* | same shared lib | adapters carry directly |
| useSeries / useMemberTraces(+Loading) / useMemberExposures / useMemberSamples / useCorpusSamples | hook | carry | queries.ts | same | feeds WaterfallRow adapter |
| useSaveSeries / useCommitSeriesPlate | hook | carry | queries.ts | same | commit 409 → ConflictError via bridge |
| GET /api/series/:id; POST /series; PATCH /series/:id; POST /:id/commit | endpoint | carry | api.ts | same | expected_content_hash discipline preserved |
| seriesDraft Zustand slice | hook | carry | state.ts | same | sessionStorage persistence; draft.id===series.id coercion |
| SSE series_save / series_commit | sse | carry | replayCoordinator.ts | same | inherited from queue |
| SeriesCommitConflictModal + ConflictModalShell | component | drop | components/* | — | **CONFLICT: this survey says RETIRED (LWW, ⛔); app-shell/data-plane say CARRY. See Risks.** |
| Track-reflections / Heatmap representation | capability | drop | CrossTraceTracking* · MemberHeatmapLayer | — | deferred/out-of-scope; no placeholder control |
| Fig-caption / fig-kicker tag row | component | refactor | SeriesBuilderPage | print/components/SeriesPlate · PlateHeader | assemble text into foot/kicker slots |
| Empty-state / Skeleton / showPeakTicks·Labels | component/hook | carry | SeriesBuilderPage · state.ts | print/ui/EmptyState · bones | fixture shape needs update |

---

## Frontend — app shell, router, SSE, data plane, mutation queue

### App shell + router + SSE lifecycle

The legacy shell roots at `main.tsx` (StrictMode → ErrorBoundary → QueryClientProvider → BrowserRouter → App). `App.tsx` mounts four singletons: EventSource(`/api/events`, "curation"), `attachPersistence`, `attachConflictBridge` (series_commit only), one-shot `rehydrate`. `AppRoutes` defines one `CorpusShell` layout route hosting seven children + redirect/stale handling. The greenfield has TopBar/StageTabs/Wordmark primitives but no wired router, no App.tsx, no SSE mount — Layer 4 is entirely ⬜.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| main.tsx (providers, QueryClient, boneyard) | component | carry | main.tsx | — | staleTime 30s/retry 1/no-focus-refetch is a product decision; oklch literal must stay |
| ErrorBoundary | component | carry | ErrorBoundary.tsx | — | inline style exempt from guard |
| App.tsx (4 effects + layout) | component | refactor | App.tsx | — | Rebuild as greenfield App.tsx; swap mounted components as they land |
| EventSource('/api/events') | sse | carry | App.tsx:34 | — | single ES, 'curation' listener only |
| AppRoutes (full table) | route | refactor | components/AppRoutes.tsx | — | 6 pages need greenfield L4; redirects carry structurally |
| CorpusShell | component | refactor | components/CorpusShell.tsx | — | TopBar + Outlet; data-testid='corpus-shell' tested |
| CorpusTopbar | component | refactor | components/CorpusTopbar.tsx | print/ui/TopBar.tsx | 305-line; wire StageTabs onChange→navigate, derive active from pathname |
| StageTabs / TopBar / Wordmark | component | new | — | print/ui/* | Built Batch 1; not router-wired |
| NavModal | component | carry | components/NavModal.tsx | — | app-level; uses legacy ModalShell — import swap |
| OnboardingFlow | component | carry | components/OnboardingFlow.tsx | — | username gate; calls api.listUsers/createUser directly |
| SeriesCommitConflictModal | component | carry | components/SeriesCommitConflictModal.tsx | — | **CARRY here (series 409 live). CONFLICT w/ Samples & Series-builder surveys (drop).** |
| InfrastructureBanner | component | carry | components/InfrastructureBanner.tsx | — | token classes only; not ported |
| ToastContainer | component | carry | components/ui/Toast.tsx | print/ui/Toast.tsx | verify print Toast exports a container host |
| handleRemoteEvent / replayCoordinator | hook | carry | lib/queue/replayCoordinator.ts | — | most load-bearing queue module; pure lib |
| attachPersistence / attachConflictBridge | hook | carry | lib/queue/* | — | SCHEMA_VERSION=4; series_commit-only bridge |
| useAppState (Zustand) | hook | carry | state.ts | — | persist v5; migrate strips activePage/theme |
| useGlobalShortcuts / useFocusTrap / useStateFromUrl / useSyncActiveSampleFromRoute / useAutoPickExposure / useCurrentUserId | hook | carry | hooks/* | — | preserve hoisting + pathname-aware step shortcuts |
| queryKeys + all query/mutation hooks | hook | carry | queries.ts | — | keys are the SSE-invalidation contract; comparison hooks dormant |
| api.ts | endpoint | carry | api.ts | — | not UI-specific |
| StaleUrlPage / ResolvingFallback / IndexSlugRedirect | component | carry | components/* | — | catch-all machinery; one-frame race guard load-bearing |
| postSampleMessage/postComparisonMessage mutators | hook | drop | lib/queue/mutators/trivial.ts | — | chat retired; drop unless data plane reinstated |

### Data plane — queries.ts / api.ts / lib/queue

The mutation engine (`useQueueMutation` + `handleRemoteEvent` + `applyRemoteToCache` + `replayCoordinator` + `deferred` + `persistence`) survives intact. ~Half the hooks (core SAXS + series CRUD + scoping + corpus) carry. The retired-feature hooks (comparisons, comparison pins/messages, sample messages) are dead and must be dropped. `ConflictError` survives for the series_commit path. `seriesPins` queryKey + SSE arm exist but no hook is exported.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| useQueueMutation / handleRemoteEvent / applyRemoteToCache / deferred / persistence / mutatorRegistry / queue hooks | hook | carry | lib/queue/* | — | engine intact; comparison_* registry arms removed when mutators drop |
| queryKeys | hook | carry | queries.ts:51 | — | comparison keys stay until SSE arms pruned |
| Core read hooks (experiments, samples, corpus, exposures, trace, members, peaks, indices, assignment) | hook | carry | queries.ts | — | all greenfield pages consume |
| Core mutation hooks (peak add/remove/exclude, reanalyze, assignment add/remove/state, custom-index, speculative, deleteIndex, updateSample, tags, status, select) | hook | carry | queries.ts | — | all queue-backed |
| Series hooks (list/get/scopeSeries/save/commit/delete) + corpus tags/picker | hook | carry | queries.ts | — | folio/builder/scoping |
| usePickerSamples / useSampleTags / useRecentlyPickedExposures | hook | carry **or** drop | queries.ts | — | **CONFLICT: data-plane survey = carry (future picker); dead-code sweep = drop (no live consumer). Experiment-scoped variants are the drop set.** |
| useSampleMessages / usePostSampleMessage | hook | drop | queries.ts:559 | — | chat retired |
| useComparisons/useComparison/Forks/Messages/PostMessage/Save/Delete/Pins/Pin/Unpin | hook | drop | queries.ts:699-942 | — | Compare retired; routes redirect |
| saveComparisonMutator / deleteComparisonMutator / postSample·ComparisonMessageMutator | hook | drop/refactor | lib/queue/mutators/* | — | registered for rehydration+SSE replay; remove registry arms together |
| contentHash.ts (canonicalJson/contentHash/comparisonHash8) | hook | drop | lib/comparison/contentHash.ts | — | zero prod callers; Julia canonical_json must relocate first |
| applyRemoteToCache comparison SSE arms | hook | drop | applyRemoteToCache.ts:181-266 | — | unreachable; drop with mutators (NOT before backend stops emitting) |
| ConflictError + conflictBridge + SeriesCommitConflictModal | component | carry | api.ts:628 · lib/queue/conflictBridge.ts | — | series_commit only |
| ConflictModalShell | component | refactor/gap | components/ConflictModalShell.tsx | — | not ported to greenfield; refactor onto print ModalShell |
| useFigureExport / useDragReorder | hook | new | — | print/components/* | greenfield-only |
| queryKeys.seriesPins (no hook) | hook | gap | queries.ts:102 | — | key+SSE arm exist; useSeriesPins/usePin/useUnpin missing |
| useSeriesForksOf | hook | gap | api.forksOfSeries:1013 | — | api fn exists; no hook |
| CorpusSample.screened / .phase | endpoint | gap | api.ts:49 | print/components/SampleTableRow.tsx | backend #162 not wired; derived client-side |
| lib/comparison/coloring·draft·normalization·yBands·labels·dragThreshold | hook | carry | lib/comparison/* | — | consumed by surviving Series components |
| lib/sample/screened.ts | hook | carry | lib/sample/screened.ts | — | derive until #162 lands |

---

## Backend

### API route surface

Routes register in `server.jl` via `register_routes!()`. Comparison HTTP routes are entirely absent (routes_comparisons.jl deleted; comparisons.jl helper remains). Sample/series message routes exist but have zero greenfield consumers. The 409 gate survives only on POST /api/series/:id/commit. The export endpoint uses legacy index_groups. Phase-4: drop ~8 dead endpoints, refactor 2, add 2 net-new.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| GET /api/health · /api/events (SSE) | endpoint/sse | carry | server.jl:74,84 | same | core |
| GET/POST /api/users; /:username/actions | endpoint | carry | routes_users.jl | same | onboarding + audit surface |
| Experiment routes (list/get/PATCH-stub/analyze/reingest) | endpoint | carry | routes_experiments.jl | same | PATCH always-400 stub |
| Samples routes (list/corpus/PATCH/tags/batch/:id) | endpoint | carry | routes_samples.jl | same | scoping batch contract tested |
| Exposures routes (list/image/status/select/tags/analyze/:id) | endpoint | carry | routes_exposures.jl | same | select LWW intentional |
| Trace / peaks / indices / assignment / speculative / custom-index routes | endpoint | carry | routes_trace·peaks·analysis.jl | same | custom-index 2-frame ordering load-bearing |
| GET /api/experiments/:id/export | endpoint | refactor | routes_export.jl:91 | same | reads index_groups → must JOIN assignment_members (state='indexed') |
| Series routes (list/get/forks/POST/PATCH/commit/delete) | endpoint | carry | routes_series.jl | same | commit 409 RETAINED for series |
| GET/POST /api/series/:id/messages | endpoint | drop | routes_series.jl:291 | — | series chat retired |
| GET /api/users/me/series-pins; POST/DELETE /:id/pin | endpoint | gap | routes_series.jl:358 | — | routes IMPLEMENTED; frontend hooks missing |
| GET/POST /api/samples/:id/messages | endpoint | drop | routes_messages.jl | — | sample chat retired |
| /api/comparisons/* (all) | endpoint | drop | (no handler file) | — | HTTP surface already absent; orphaned frontend hooks |
| GET /api/users/:id/recently-picked-exposures | endpoint | refactor | routes_picker.jl:33 | same | queries comparison_members → port to series_members |
| GET /api/sample-tags · /picker-samples (corpus) | endpoint | carry | routes_picker.jl:77,99 | same | used by scoping |
| GET /api/experiments/:eid/sample-tags · /picker-samples | endpoint | carry **or** drop | routes_picker.jl:61,93 | same | **CONFLICT: API survey = carry (builder picker); dead-code sweep = drop (no live caller).** |
| GET /api/resolve | endpoint | carry | routes_resolve.jl:38 | same | permalink resolution |

### Schema + pipeline + data model

Three layers: (1) core pipeline tables, (2) assignment model (assignments live; index_groups legacy write-only), (3) collaboration/series. Cutover keeps pipeline+assignment+series; drops legacy comparison/message plumbing once migration relay is no longer needed.

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| experiments / samples / sample_tags / exposures / exposure_sources·tags | table | carry | db.jl | — | fully live |
| auto_peaks / peak_curations | table | carry | db.jl:100,110 | — | trace-plot engine load-bearing |
| indices / index_peaks | table | carry | db.jl:75,88 | — | inputs_hash drives StaleIndicesBanner |
| assignments / assignment_members | table | carry | db.jl:142 | — | backs confirmed_index in snapshot |
| index_groups / index_group_members | table | refactor | db.jl:124 | — | legacy write target + migration relay; export must move off them |
| series / series_members / series_samples | table | carry | db.jl:759 | series.jl | confirmed_index sources assignment_members (D-10) |
| series_messages | table | drop | db.jl:824 | — | chat retired; DROP deferred (FK CASCADE) |
| series_pins | table | carry | db.jl:838 | — | pin routes live |
| sample_messages | table | drop | db.jl:154 | — | routes_messages.jl retired |
| comparisons / comparison_members / comparison_messages | table | refactor | db.jl:488 | comparisons.jl, series.jl | **permanent replay machinery — do NOT drop**; rebuild_views re-folds events |
| comparison_pins | table | drop | db.jl:569 | — | migrated to series_pins; dispatcher branches stay |
| user_actions | table | carry | db.jl:165 | — | event log; undoes_event_id reserved for undo/redo |
| idempotent_responses / schema_migrations | table | carry | db.jl:199,852 | idempotency.jl | sentinel-gated migrations |
| analyze_exposure! / persist_analysis! | component | carry | pipeline.jl | — | hash-guarded skip paths; seed_assignment_if_absent! |
| migrate_assignments! / migrate_comparisons_to_series! | migration | carry | db.jl:1069,978 | — | sentinel-gated; pre-Plan-A state not log-derivable |
| insert_speculative/custom_index! / Gauss-Bonnet predictor / rebuild_views_from_log! | component | carry | speculative.jl · routes_analysis.jl · events.jl | — | — |
| edit-tracking / undo-redo / versioning tables | table | new | — | — | not yet defined; undoes_event_id is the hook |

### Commit contract (events / SSE / idempotency / content-hash gate)

`apply_event!` (atomic append + view-update, InTransaction variant), `broadcast_event!` (post-commit SSE fan-out carrying client_id/client_op_id/post_state), `with_idempotency` (per-op-id dedup). The only 409 gate is conditional, in POST /api/series/:id/commit (routes_series.jl:225-244).

| item | kind | class | legacy ref | greenfield ref | notes |
|---|---|---|---|---|---|
| apply_event! dispatcher (all view kinds) | event | carry | events.jl:42 | — | InTransaction variant is the route form |
| peak_added/excluded/unexcluded/removed | event | carry | events.jl:311 | — | post_state envelope; greenfield trace editor |
| assignment_add/remove/set_state | event | carry | events.jl:369 | — | greenfield assignment rail |
| index_confirmed/unconfirmed | event | carry | events.jl:363 | — | retired but MUST remain as no-op guards |
| speculative_created/deleted | event | carry | events.jl:416 | — | custom-index + speculative |
| analyze_run | event | carry | events.jl:422 | — | M0.4 SSE suppression when both skips true |
| update_sample/add_tag/remove_tag/set_exposure_status/select_exposure | event | carry | events.jl:404 | — | M2.1 trivial routes; all live |
| post_message (sample/series) | event | drop | events.jl:408 | — | keep exhaustiveness branch; routes droppable |
| comparison_created/submitted/deleted/pinned/unpinned | event | drop | events.jl:427 | — | **keep branches for replay; comparisons.jl helpers still shared** |
| series_created/recipe_updated/deleted/plate_committed/pinned/unpinned | event | carry | events.jl:469 | — | plate_committed is the only post_state series kind |
| broadcast_event! + SSE stream | sse | carry | events.jl:1027 · server.jl:84 | — | per-subscriber Channel(64) + 15s heartbeat |
| with_idempotency + gc | capability | carry | idempotency.jl:83 | — | Stripe-style; I2 atomic event+cache write |
| post-commit broadcast queue | capability | carry | events.jl:217 | — | task-local; flush on commit, clear on rollback/4xx |
| expected_content_hash / 409 gate | refactor | refactor | routes_series.jl:225 | — | CONDITIONAL; relax = delete lines 225-227 + 236-244 |
| hash_trace_file / hash_peak_set | capability | carry | hash.jl:10 | — | drives skip predicates; distinct from content_hash |
| edit-tracking → undo/redo → versioning | capability | new | — | — | hook-in: new 'revert_to_event' kind; user_actions.undoes_event_id is the FK slot |
| routes_messages.jl | endpoint | drop | routes_messages.jl:24 | — | de-register + delete |
| comparisons.jl business-logic file | component | refactor | comparisons.jl | series.jl, routes_picker.jl | shared helpers must relocate first |
| rebuild_views_from_log! | capability | carry | events.jl:1075 | — | needs a branch + test for every new kind |

---

## Consolidated GAPS (greenfield/backend pieces still missing)

Ranked by cutover criticality, deduped across all areas.

1. **[frontend] Six Layer-4 page shells do not exist** — SamplesPage, LoupePage, FocusPage, SeriesFolioPage, SeriesScopingPage, SeriesBuilderPage. All L0-L3 composites are built; each page wires queries+state+navigation+Skeleton. This is the bulk of the cutover.
2. **[frontend] Greenfield App.tsx + AppRoutes + CorpusShell do not exist** — no SSE mount, queue persistence, conflict bridge, rehydrate, or router wiring in the greenfield tree.
3. **[frontend] SeriesRecipeEditor greenfield composite** — the largest single behavioral gap for the builder: title/desc/add-sample/order/Save/Commit wired to useSaveSeries + useCommitSeriesPlate; BuilderRail only exposes a presentational traces slot.
4. **[frontend] StaleIndicesBanner, SpeculativeBuilder, FocusNotesMargin** — three Focus behaviours with no print/ port (hash-stale debounce; anchor+ratio speculative dialog; q-ref notes + focus-gated textarea + drawer).
5. **[backend→frontend] series-pins hooks missing** — GET /api/users/me/series-pins + POST/DELETE /:id/pin are implemented backend; no useSeriesPins/usePinSeries/useUnpinSeries, no folio pin affordance.
6. **[backend] CardFigure trace-data source** — folio cards need real WaterfallRow trace data; useSeries returns only members. Needs GET /api/series/:id/traces (or embedded traces map) or cards render beads-only. **Highest-priority new backend route implied by greenfield.**
7. **[frontend] q-link / TracePlot onHoverQ emission** — deferred Layer-1 addition; trace-peak hover → ring → tooth circuit incomplete.
8. **[frontend] StatusCell corpus→focus door + ThumbnailGallery dbl-click/shift-click/loupe + X/Esc cull shortcuts** — contact-sheet interaction gaps; loupe unreachable from sheet without them.
9. **[frontend] losingPeakIds dim** — TracePlot has highlightPeakIds (inverse); extend or invert in page.
10. **[frontend] GroupingModeToggle / AnnotationToggles / SeriesPhaseStrip companion** — builder controls with no print/ composite; decision needed on keep-vs-drop for grouping mode.
11. **[frontend] greenfield ScopingConfirmModal + ConflictModalShell port** — both small; ModalShell ready.
12. **[frontend] Exposure→GalleryExposure src adapter + WaterfallRow adapter** — page-level pure helpers; URL-pattern must match legacy for cache coherence; WaterfallRow adapter replicates buildMemberMarks coloring/peak/band logic.
13. **[frontend] boneyard skeleton fixtures/bones** — contact-sheet recapture, loupe.bones.json, scoping bone, folio bone, builder fixture all need (re)capture for new layouts.
14. **[frontend] InfrastructureBanner, NavModal, OnboardingFlow, StaleUrlPage, ResolvingFallback, IndexSlugRedirect, ToastContainer** — not ported to greenfield (mostly import swaps).
15. **[backend] recently_used_exposures series-side query** — currently queries comparison_members.created_by; hollows out as comparison data ages.
16. **[backend] export endpoint** — reads index_groups (pre-D10); silently wrong for assigned exposures.
17. **[backend] edit-tracking/undo-redo/versioning schema** — no tables; no revert_to_event kind; no GET /api/exposures/:id/history; no inverse-dispatch logic.
18. **[backend] comparisons.jl shared-helper relocation** — series.jl/events.jl depend on canonical_json, compute_member_snapshot, is_member_stale, recently_used_exposures, picker_samples, comparison_now_iso, etc.; not yet moved.
19. **[frontend] CorpusSample.screened (#162) / .phase** — derived client-side; backend fields not wired.
20. **[backend] no bulk exposures endpoint** — useCorpusExposures fans out N requests (bounded by useQueries; not a blocker).
21. **[frontend] useSeriesForksOf hook** — api.forksOfSeries exists; no hook.
22. **[frontend] canonicalJson has no prod callers** — if cross-language hash parity matters later, move to a dedicated utilities module.

---

## DEAD CODE to delete (the drop set)

### Frontend — clean drops
- Legacy contact-sheet components: ContactSheetRow, ExposureThumb, CullBar (legacy), SampleStatusChip, RejectXMark, legacy DetectorImage/ThumbnailGallery, legacy ui Kicker/HintText/Card/SegmentedControl/PhaseStrip.
- Legacy Loupe: LoupeFrame, LoupeSidebar, LoupeTagsEditor.
- Legacy Focus: FocusWorkspaceLayout, PlotCard, TraceViewer, FocusDetectorPanel, CombPanel, IndicesCard, PhasePanel, legacy CustomIndexModal, FigureExportControls.
- Legacy Series folio/builder: SeriesFolioCard, SeriesMiniWaterfall, buildMiniWaterfall/folioFigure.ts, MultiTracePlot, offsetToBandFraction, SeriesBuilderRail, AutogroupCard, RepresentationToggle, OffsetSlider, ScaleToggle, OffsetDock, rail-restore, SeriesReadingPanel, legacy SeriesMemberRow, MemberMetaGutter, BandResizeDivider, ActiveBandContext.
- Legacy scoping: ScopingRow, ScopingValueCell, ScopingSparkline, ScopingAutogroupCard, ScopingOrderField, ScopingLooseMatches, ScopingFoot, lib/plot/sparkline.ts.
- Retired-feature hooks: useSampleMessages, usePostSampleMessage, useComparisons/useComparison/Forks/Messages/PostMessage/Save/Delete/Pins/Pin/Unpin.
- Drop set tests: queue/saveComparison.test.tsx, queue/deleteComparison.test.tsx, queue/applyRemoteToCache.compare.test.ts, queries-messages.test.tsx.
- applyRemoteToCache comparison SSE arms + comparison queryKeys (after mutator cleanup).

### Frontend — shared, EXTRACT-FIRST (do not delete naively)
- **postSampleMessage/postComparisonMessage/saveComparison/deleteComparison mutators** — registered for sessionStorage rehydration + SSE replay; remove registry arms together and bump SCHEMA_VERSION so stale blobs are ignored. The shared `post_message` OpKind case must keep its surviving arm.
- **contentHash.ts** — has zero prod callers but the cross-language parity test pairs with Julia canonical_json; move canonicalJson to a utilities module + repoint the test, or delete atomically with the Julia fixture.

### Backend — clean drops
- routes_messages.jl + register_messages_routes!() in server.jl + include() in HimalayaUI.jl + test/test_routes_messages.jl.
- GET/POST /api/series/:id/messages route registrations in routes_series.jl.

### Backend — shared, EXTRACT-FIRST / KEEP-FOR-REPLAY (do not delete naively)
- **comparisons.jl** — NOT dead. Shared helpers (canonical_json, comparison_now_iso, compute_member_snapshot, is_member_stale, recently_used_exposures, picker_samples, _topk_phases, _count_distinct_phases, user_id_for_event, _ngc_for_phase) are consumed by series.jl + routes_picker.jl + events.jl. compute_content_hash + member_ids_for_comparison are used only by the comparison dispatcher arms. Relocate shared helpers to a new helpers file BEFORE deleting; the comparison-specific arms stay for replay.
- **comparison* tables (comparisons, comparison_members, comparison_messages, comparison_pins)** — KEEP. rebuild_views_from_log! re-folds historical comparison events. Dropping them breaks replay on any DB with comparison history.
- **comparison event dispatcher branches in events.jl** — KEEP for replay exhaustiveness.
- **series_messages / sample_messages tables** — drop deferred; FK CASCADE / SET NULL requires a PRAGMA foreign_keys=OFF rebuild migration.

> **Explicit conflict — SeriesCommitConflictModal + ConflictModalShell:** the Samples and Series-builder surveys classify these **drop** (asserting the 409/conflict UI is cancelled and series commit is LWW). The App-shell, Data-plane, and both backend commit-contract surveys classify them **carry** (the series_commit 409 gate is explicitly RETAINED; only the comparison If-Match flow was cancelled). The backend confirms the 409 gate is still live in routes_series.jl:225-244. **Treat as carry/retain until a deliberate decision relaxes the gate** (see next section) — do not drop the modal before the gate is relaxed and the frontend stops handling 409.

> **Explicit conflict — usePickerSamples / useSampleTags / useRecentlyPickedExposures + experiment-scoped picker routes:** the data-plane and API-route surveys mark these **carry** (future series-builder picker); the dead-code sweep marks them **drop** (no live consumer today). Corpus-wide siblings (GET /api/sample-tags, GET /api/picker-samples) are unambiguously **carry** (used by scoping). Resolve by deciding whether the greenfield builder gets an experiment-scoped picker + "recently used" section.

---

## Backend refactors & NET-NEW features

### Refactors
- **Relax expected_content_hash / 409 gate (DECIDED)** — the only live optimistic-concurrency gate, conditional, in POST /api/series/:id/commit (routes_series.jl:225-244). Relax = delete the expected_hash extraction (225-227) + the 409 response branch (236-244) + the current_series_content_hash read at 237. KEEP compute_series_content_hash / current_series_content_hash (still write content_hash for series_plate_committed post_state + future fork/stale checks). Must be coordinated with the frontend ceasing to send/handle the field. **Cross-survey caveat:** several frontend surveys still describe SeriesCommitConflictModal as live — this relaxation is the deliberate act that retires it; do it atomically with the modal removal.
- **build_export()** — replace index_groups/index_group_members SELECT with a JOIN on assignment_members WHERE state='indexed'; rename/update :active_group_kind → :assignment_state.
- **recently_used_exposures()** — port from comparison_members.created_by to series_members.created_by (or a transition union).
- **comparisons.jl helper relocation** — move shared helpers to a neutral helpers file so comparisons.jl can eventually be deleted without breaking series.jl (include-order coupling is load-bearing).
- **CorpusTopbar / StageTabs** — wire onChange→navigate, derive active from pathname (presentational primitive → router-aware shell).
- **ConflictModalShell / SeriesRecipeEditor / FigureExportControls** — refactor onto print primitives.

### Net-new features
- **GET /api/series/:id/traces (or embedded traces map)** — highest-priority new route: a single exposure_id→Trace map per series so CardFigure renders real curves without N×M fan-out across the folio masonry.
- **series-pins frontend wiring** — useSeriesPins/usePinSeries/useUnpinSeries hooks over already-implemented routes; folio pin affordance.
- **Edit-tracking → undo/redo → versioning subsystem (DECIDED, replaces cancelled 409 UI; Layer 4)** — concrete hook-in detail from the commit-contract survey: (1) user_actions is already an append-only ordered log with undoes_event_id (already set on peak_unexcluded — the undo-pointer schema exists); (2) apply_event! InTransaction + with_idempotency give dedup; (3) broadcast carries actor+ts. MISSING: (a) a new event kind `revert_to_event {target_event_id}` whose dispatcher re-folds the log to that point for the affected entity (the first concrete hook-in point); (b) per-kind inverse-dispatch logic (non-trivial for multi-step kinds like series_plate_committed); (c) a version-snapshot concept/table (content_hash on series is a partial form covering only the committed plate, not peaks/assignments); (d) GET /api/exposures/:id/history returning the user_actions log; (e) frontend hooks presenting an undo stack from recent user_actions rows.
- **Optional GET /api/samples/bulk-exposures** — eliminate the useCorpusExposures fan-out (deferred; bounded by useQueries).

---

## Risks & cutover-ordering constraints

Deduped across all surveys.

- **Atomic page swaps** — every route must be swapped atomically in AppRoutes (import legacy ↔ greenfield). Partial swap (e.g. replacing only TracePlate) is unsafe because greenfield panels expect page-owned state the legacy layout does not provide.
- **SSE-after-commit ordering** — broadcast must always fire AFTER the SQLite transaction commits (_flush_post_commit_broadcasts! at the end of with_idempotency). Reordering lets subscribers see uncommitted state.
- **useCorpusExposures O(1)-observer pattern is load-bearing** — do NOT regress to per-row useExposures; the fan-out caused the ERR_INSUFFICIENT_RESOURCES bug on the 139-sample corpus.
- **SeriesCommitConflictModal cross-survey conflict** — the series_commit 409 gate is RETAINED (backend-confirmed). Do not drop the modal until the gate is deliberately relaxed AND the frontend stops handling 409. Removing the modal mount from App.tsx while legacy SeriesBuilderPage is still wired risks unhandled ConflictError during the transition.
- **comparisons.jl include-order** — Julia loads comparisons.jl before series.jl; deleting it before moving shared helpers produces UndefVarError at module load. Relocate first.
- **comparison* tables must not be dropped** — rebuild_views_from_log! re-folds historical comparison events; dropping breaks replay on existing DBs. Same for keeping comparison dispatcher branches.
- **rebuild_views_from_log! exhaustiveness** — every new event kind (undo/redo/versioning) needs a dispatcher branch + a round-trip test or disaster recovery silently breaks.
- **Cache-coherence for Focus notes** — notes MUST read useSamples(experimentId), not useCorpusSamples; breaking this silently reverts notes after save.
- **TracePlot engine swap (Observable Plot → d3)** — deep entanglement; peak hit-testing, scroll-zoom, label placement differences are user-visible. losingPeakIds inversion risks dimming the wrong peaks.
- **defaultExposureId ordering (rep→first-accepted→first)** — load-bearing for which frame opens; mis-implementation is a subtle regression not caught by render tests.
- **Exposure→GalleryExposure URL adapter** — must exactly match the legacy ?v=/thumb pattern for cache coherence.
- **SCHEMA_VERSION / persist version bumps** — bump sessionStorage SCHEMA_VERSION when adding/removing OpKinds (stale comparison/message ops dropped with a toast); bump Zustand persist version + migrate when adding persisted fields (v5 migration strips activePage/theme).
- **StageTabs wiring** — if the L4 shell calls setActiveSample instead of navigate, the , / . step shortcuts break on /sample/:id (they switch on pathname.startsWith).
- **PageBody one-frame ResolvingFallback** — moving useStateFromUrl out of PageBody reintroduces #181 (typo'd URLs bounce to contact sheet instead of "Page not found").
- **Active-experiment gate on samples query** — preserve `enabled: experimentId !== undefined` or GET /api/samples?experiment=0 fires on cold mount.
- **CardFigure N+1 / trace-fetch** — resolve the per-series traces route before folio cutover or render beads-only; otherwise O(N×M) fetches on page load.
- **offsetToBandFraction vs offsetScale** — different formulas; calibrate so default offset≈1.2 looks identical; needs a regression test.
- **Export silently wrong** — until build_export() moves to assignment_members, exports differ from the UI for any non-auto assignment.
- **recently_used_exposures hollows out** — silent UX regression (not a crash) as comparison_members ages; fix before any "recently used" section ships.
- **Boneyard skeleton ordering** — bones.json must be (re)captured + committed before a greenfield page ships or the skeleton falls back to plain text.
- **OP_LOCKS / lock-ordering** — OP_LOCKS_MU must never be held across _DB_WRITE_LOCK; long-lived undo op_ids may keep stale locks past the 1h GC TTL.
- **Tag.id optimistic-add** — TagList onRemove receives a Tag object; if an optimistic-add lacks id the remove handler gets undefined — test the add→remove path.
- **Keyboard input/textarea guard** — verify TagList's editor inputs report INPUT/TEXTAREA so loupe/focus shortcuts don't clobber tag editing.
