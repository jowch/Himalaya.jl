# Changelog

## Unreleased

### Breaking changes

- **`samples.label` removed; replaced by `name` (stable identifier, was column 3) and
  `display_name` (editable label, was column 2). Manifest column meanings swapped.**
  Existing experiments need `himalaya migrate-toml <experiment-dir>` to upgrade
  their `experiment.toml`. Existing DBs auto-migrate on first `open_db` after
  deploy. Issue #88.
- **`/api/export` CSV header changed: `sample_label` → `sample_display_name`.**
  Downstream pipelines parsing this CSV need to update their column names.
- **`PATCH /api/samples/:id` no longer accepts `name`; use `display_name`.**
  `samples.name` is now the stable identifier and is set only at ingest time.
- **`PATCH /api/experiments/:id` no longer accepts any field.** The route is
  retained as defensive surface for future fields. Experiment renames must go
  through `experiment.toml` + reingest.
- **First boot after migration purges `idempotent_responses`.** In-flight
  `client_op_id` retries from before the deploy will get fresh executions
  rather than cached pre-rename response bodies.
