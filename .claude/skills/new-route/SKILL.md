---
name: new-route
description: Scaffold a new REST resource in HimalayaUI — Julia route handler, json.jl helpers if needed, api.ts fetcher, and queries.ts hook. Read existing routes as live reference so the skill stays accurate after refactors. Usage: /new-route <resource> [--parent <parent>]
---

# new-route

Scaffolds the four-file pattern for a new REST resource in HimalayaUI. Rather than relying on static examples, the skill reads current source files as its reference — so it stays accurate even after the codebase is restructured.

## Args

- `<resource>` — the new resource name in snake_case (e.g. `tags`, `audit_entries`)
- `--parent <parent>` — parent resource if nested (e.g. `--parent samples` → `/api/samples/{id}/tags`)
- `--readonly` — omit mutation handlers and mutation hooks (GET only)

## Procedure

### 1. Read the current reference implementations

Before writing any code, read these files to understand the current conventions:

```
packages/HimalayaUI/src/routes_peaks.jl      ← canonical CRUD route file
packages/HimalayaUI/src/json.jl              ← row serialization helpers
packages/HimalayaUI/src/server.jl            ← route registration wiring
packages/HimalayaUI/src/actions.jl           ← log_action! signature
packages/HimalayaUI/frontend/src/api.ts      ← typed fetcher pattern
packages/HimalayaUI/frontend/src/queries.ts  ← queryKeys + useQuery/useMutation pattern
```

The `routes_peaks.jl` file is the best all-around reference: it has GET (list), POST (create), PATCH (update), and DELETE, all with correct SQLite/Oxygen/logging patterns. If the resource you're adding is simpler, use whichever existing route most closely matches.

### 2. Identify what handlers are needed

Decide which HTTP verbs apply to the new resource, following REST conventions:
- `GET /api/<parents>/{id}/<resource>` — list
- `POST /api/<parents>/{id}/<resource>` — create
- `GET /api/<resource>/{id}` — get one (if needed)
- `PATCH /api/<resource>/{id}` — update
- `DELETE /api/<resource>/{id}` — delete

### 3. Write `packages/HimalayaUI/src/routes_<resource>.jl`

Following exactly the patterns you read in step 1:

- **Function name**: `register_<resource>_routes!()`
- **Imports**: copy the `using` line from the reference — `HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite`
- **Reads**: `Tables.rowtable(DBInterface.execute(db, sql, params))` — never access raw cursor values, they vanish after the query closes
- **Inserts**: capture `res = DBInterface.execute(...)` then `id = Int(DBInterface.lastrowid(res))` — `lastrowid` takes the result, not the db
- **Mutations**: call `log_action!(db, req; action = "...", entity_type = "...", entity_id = ...)` after every write
- **Responses**: `HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(...))`
  - 200 for GET/PATCH
  - 201 for POST (include the newly-created row)
  - 204 for DELETE (no body)
  - 404 with `Dict(:error => "... not found")` for missing rows
- **JSON helpers**: use `rows_to_json` (list) and `row_to_json` (single row) from `json.jl`; pass `bool_keys` for any columns stored as 0/1 integers that should be booleans in JSON
- **FK on users**: if the new table has a `created_by` column referencing `users(id)`, the schema DDL must include `ON DELETE SET NULL` — not at call sites

### 4. Register in `server.jl`

In `register_routes!()`, add `register_<resource>_routes!()` alongside the existing calls. Match the ordering (grouped by resource family). Also add the include line at the top of `HimalayaUI.jl` if it isn't there.

### 5. Add JSON helpers to `json.jl` if needed

Only if the new resource needs custom serialization beyond what `rows_to_json`/`row_to_json` already provide. Check `json.jl` first — the helpers are generic and often sufficient.

### 6. Write the TypeScript fetcher in `frontend/src/api.ts`

Following the pattern from the file you read in step 1:

- Export named async functions: `export const list<Resource>s = ...`, `export const create<Resource> = ...`
- Use the `request<T>()` helper — it handles JSON parsing, error wrapping, and auth headers
- Reads get no `AuthOpts`; mutations get `opts?: AuthOpts` as the last parameter and pass it to `request()`
- Define an interface for the resource shape above the functions

### 7. Write the TanStack Query hook in `frontend/src/queries.ts`

Following the pattern from the file you read in step 1:

- Add a key to `queryKeys` at the top: `<resource>s: (parentId: number) => [..., parentId, '<resource>s'] as const`
- Export `use<Resource>s(parentId)` returning `useQuery({ queryKey: queryKeys.<resource>s(parentId), queryFn: ... })`
- For mutations, export `use<Verb><Resource>()` returning `useMutation({ mutationFn: ..., onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.<resource>s(parentId) }) })`
- Use `authOpts(username)` (already in `queries.ts`) for mutation auth — never pass `{ username: undefined }`

### 8. Write tests

- **Backend**: add a `@testset` in `packages/HimalayaUI/test/` — follow the pattern in the nearest existing test file. Use a temp DB (`open_db(joinpath(tmp, "test.db"))`) and `start_test_server!` / `stop_test_server!`.
- **Frontend**: add `test/<Resource>.test.tsx` (or extend an existing test file) — mock the fetcher at the module boundary, not at the fetch level.

### 9. Verify

```bash
# Backend
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'

# Frontend types
cd packages/HimalayaUI/frontend && npx tsc --noEmit -p tsconfig.build.json

# Frontend tests
npx vitest run test/<Resource>.test.tsx
```

## Gotchas

- **`Tables.rowtable` is mandatory for reads.** Raw `DBInterface.execute` rows silently lose their values after the query closes. Always materialize with `Tables.rowtable(...)` before accessing fields.
- **`lastrowid` takes the result, not the db.** `DBInterface.lastrowid(res)`, not `DBInterface.lastrowid(db)`.
- **Stdlib deps must be in `Project.toml`.** If the new route file uses a stdlib not already listed (`Printf`, `Sockets`, etc.), add it to `[deps]` in `packages/HimalayaUI/Project.toml`.
- **`exactOptionalPropertyTypes: true` in TS.** The `authOpts(username)` helper in `queries.ts` returns `{}` or `{ username }` — never `{ username: undefined }`. Use it; don't inline `{ username: someVar }`.
- **Mutations must invalidate scoped query keys.** `queryKeys.<resource>s(parentId)` — not the root `queryKeys.experiments` or similar broad key.
- **Phase serialization is unrelated but easy to mix up.** This pattern is for REST resources, not for Julia `Phase` types. If you're adding a route that returns phase data, use `string(nameof(P))` not `string(P)` — the latter returns `"Himalaya.Pn3m"` which breaks SQLite roundtrips.
