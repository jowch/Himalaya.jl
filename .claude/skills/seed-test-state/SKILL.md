---
name: seed-test-state
description: Navigate the running Vite preview to a specific experiment/sample/exposure/page by writing Zustand state to localStorage and reloading. Use when manually testing the UI on a known data path (e.g. AgBe trace, cubic sample, Inspect page on exposure 12). Skips the click-through-the-nav-modal dance.
disable-model-invocation: true
---

# seed-test-state

Quickly point the running preview at a specific app state — useful when iterating on a UI change that requires being on a specific experiment/sample/exposure/page.

## Args (in any order)

- `experiment=<id>` — experiment id (default: leave unchanged)
- `sample=<id>` — sample id (default: leave unchanged)
- `exposure=<id>` — exposure id, or `auto` to clear (default: clear, so PlotCard auto-picks)
- `page=index|inspect|compare` — which tab (default: leave unchanged)

If no args, shows current state and a list of available experiments/samples/exposures from the live backend.

## Procedure

1. **Discover the running preview server.**
   ```
   mcp__Claude_Preview__preview_list
   ```
   If no server, surface the launch.json name (`himalaya-frontend`) and tell the user to start it via `mcp__Claude_Preview__preview_start name=himalaya-frontend`.

2. **List available IDs.** Hit the live backend on `:8080` so the user can pick:
   ```bash
   curl -s http://localhost:8080/api/experiments
   curl -s http://localhost:8080/api/experiments/<id>/samples
   curl -s http://localhost:8080/api/samples/<id>/exposures
   ```
   Show id + name pairs.

3. **Write the state.** The Zustand persisted state lives in `localStorage` under the key `himalaya-ui:state` with shape `{ state: { ... }, version: 3 }`.
   ```js
   (() => {
     const raw = localStorage.getItem('himalaya-ui:state') || '{}';
     const parsed = JSON.parse(raw);
     parsed.state = parsed.state || {};
     // Apply only the args the user passed:
     if (typeof EXPERIMENT === 'number') parsed.state.activeExperimentId = EXPERIMENT;
     if (typeof SAMPLE === 'number')     parsed.state.activeSampleId     = SAMPLE;
     if (EXPOSURE === 'auto')            parsed.state.activeExposureId   = undefined;
     else if (typeof EXPOSURE === 'number') parsed.state.activeExposureId = EXPOSURE;
     if (PAGE)                           parsed.state.activePage         = PAGE;
     // Required defaults so we don't trip the onboarding overlay:
     parsed.state.username     = parsed.state.username     ?? 'test-user';
     parsed.state.tutorialSeen = true;
     parsed.state.theme        = parsed.state.theme        ?? 'dark';
     parsed.version = 3;
     localStorage.setItem('himalaya-ui:state', JSON.stringify(parsed));
     location.reload();
     return 'reloaded';
   })()
   ```
   Send via `mcp__Claude_Preview__preview_eval`.

4. **Wait ~2s** for the reload, then `mcp__Claude_Preview__preview_screenshot` so the user can confirm the seeded state.

## Examples

- `/seed-test-state experiment=1 sample=12 page=index` → AgBe trace on Index page
- `/seed-test-state page=inspect sample=12` → AgBe images on Inspect page
- `/seed-test-state exposure=auto` → reset to PlotCard's auto-picked first exposure

## Why user-only (`disable-model-invocation: true`)

Writing state to a live preview is a side-effecting operation that the user should always be the one initiating. Claude shouldn't unilaterally jump the preview to a different state mid-task.

## Gotchas

- **Version mismatch:** Zustand persist key is `himalaya-ui:state`, version `3`. If `state.ts` bumps the version, this skill breaks until updated.
- **Onboarding overlay:** if `username` isn't set, the app shows the onboarding modal and ignores other state. Skill always seeds a default username.
- **Backend down:** if `:8080` doesn't respond, exposure lists won't load — but the page state seed still works, the user just sees empty data.
