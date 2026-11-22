CREATE TABLE users (
    id         INTEGER PRIMARY KEY,
    username   TEXT UNIQUE NOT NULL,
    first_name TEXT,
    last_name  TEXT
);
CREATE TABLE experiments (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT,
    path            TEXT NOT NULL,
    data_dir        TEXT NOT NULL,
    analysis_dir    TEXT NOT NULL,
    manifest_path   TEXT,
    config          TEXT,
    experiment_type TEXT,
    image_pattern        TEXT,
    metadata_pattern     TEXT,
    integration_pattern  TEXT,
    energy_kev      REAL,
    flight_path_m   REAL,
    created_at      DATETIME DEFAULT CURRENT_TIMESTAMP,
    energy_kev_source        TEXT DEFAULT 'default',
    flight_path_m_source     TEXT DEFAULT 'default',
    beam_center_x            REAL,
    beam_center_x_source     TEXT DEFAULT 'default',
    beam_center_y            REAL,
    beam_center_y_source     TEXT DEFAULT 'default',
    pixel_size_um            REAL,
    pixel_size_um_source     TEXT DEFAULT 'default',
    q_units                  TEXT,
    q_units_source           TEXT DEFAULT 'default',
    last_scanned_at          TEXT,
    ingest_status            TEXT DEFAULT 'idle',
    ingest_error             TEXT,
    last_scan_tier           TEXT DEFAULT 'fast',
    consecutive_empty_ticks  INTEGER NOT NULL DEFAULT 0
, description TEXT);
CREATE TABLE loads (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    load_index    INTEGER NOT NULL,
    session_id    INTEGER,
    start_time    TEXT,
    end_time      TEXT,
    frame_count   INTEGER NOT NULL DEFAULT 0,
    note          TEXT
);
CREATE TABLE samples (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER REFERENCES experiments(id),
    name  TEXT,
    notes         TEXT,
    load_id          INTEGER REFERENCES loads(id),
    slot_index       INTEGER,
    grouping_source  TEXT DEFAULT 'auto_position',
    name_source      TEXT DEFAULT 'auto'
, merged_into_id INTEGER REFERENCES samples(id));
CREATE TABLE exposures (
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id            INTEGER REFERENCES samples(id),
    filename             TEXT,
    kind                 TEXT DEFAULT 'file',
    selected             BOOLEAN DEFAULT FALSE,
    status               TEXT CHECK (status IN ('accepted', 'rejected')),
    image_path           TEXT,
    trace_hash           TEXT,
    analysis_inputs_hash TEXT,
    experiment_id        INTEGER REFERENCES experiments(id) ON DELETE CASCADE,
    prp_path             TEXT,
    timestamp            TEXT,
    exposure_time        REAL,
    horizontal_position  REAL,
    scan_id              INTEGER,
    frame_no             INTEGER,
    load_id              INTEGER REFERENCES loads(id),
    content_fingerprint  TEXT,
    transmitted_counts   REAL,
    incident_monitor     REAL,
    counting_time        REAL,
    monitors_read        INTEGER NOT NULL DEFAULT 0,
    analysis_version     TEXT,
    analysis_error       TEXT,
    analysis_trace_stat  TEXT
);
CREATE TABLE exposure_sources (
    averaged_exposure_id INTEGER REFERENCES exposures(id),
    source_exposure_id   INTEGER REFERENCES exposures(id),
    role                 TEXT DEFAULT 'signal',
    PRIMARY KEY (averaged_exposure_id, source_exposure_id)
);
CREATE TABLE indices (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    phase       TEXT NOT NULL,
    basis       REAL NOT NULL,
    score       REAL,
    r_squared   REAL,
    lattice_d   REAL,
    status      TEXT DEFAULT 'candidate',
    kind        TEXT NOT NULL DEFAULT 'auto',
    inputs_hash TEXT,
    -- 1 = `basis`/`lattice_d` are a lattice the USER chose (the custom-index
    -- modal), not a fit. `_persist_analysis_inner!` re-resolves a locked
    -- index's peaks like any other speculative but skips the least-squares
    -- refit, so no reanalysis can move the lattice out from under the user.
    -- This is what makes insert_custom_index!'s scan-derived intents safe:
    -- they resolve WHICH peaks are claimed, never WHERE the comb sits.
    basis_locked INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE index_peaks (
    index_id       INTEGER REFERENCES indices(id),
    peak_id        INTEGER NOT NULL,
    peak_kind      TEXT NOT NULL DEFAULT 'auto'
                   CHECK (peak_kind IN ('auto', 'curation')),
    ratio_position INTEGER,
    residual       REAL,
    PRIMARY KEY (index_id, peak_id, peak_kind)
);
CREATE TABLE auto_peaks (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id     INTEGER REFERENCES exposures(id),
    q               REAL NOT NULL,
    intensity       REAL,
    prominence      REAL,
    sharpness       REAL,
    findpeaks_index INTEGER
);
CREATE TABLE peak_curations (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL CHECK (kind IN ('exclude', 'add')),
    q           REAL NOT NULL,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE index_groups (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    kind        TEXT NOT NULL DEFAULT 'auto',
    active      BOOLEAN DEFAULT FALSE,
    created_by  INTEGER REFERENCES users(id) ON DELETE SET NULL,
    created_at  DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE index_group_members (
    group_id  INTEGER REFERENCES index_groups(id),
    index_id  INTEGER REFERENCES indices(id),
    PRIMARY KEY (group_id, index_id)
);
CREATE TABLE assignments (
    exposure_id INTEGER PRIMARY KEY REFERENCES exposures(id),
    state       TEXT NOT NULL DEFAULT 'indexed'
                CHECK (state IN ('indexed', 'form_factor', 'null'))
);
CREATE TABLE assignment_members (
    exposure_id INTEGER NOT NULL REFERENCES exposures(id),
    index_id    INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
    PRIMARY KEY (exposure_id, index_id)
);
CREATE TABLE sample_calls (
    sample_id INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
    state     TEXT NOT NULL DEFAULT 'indexed'
              CHECK (state IN ('indexed', 'form_factor', 'null'))
);
CREATE TABLE user_actions (
    id              INTEGER PRIMARY KEY,
    user_id         INTEGER REFERENCES users(id) ON DELETE SET NULL,
    timestamp       DATETIME DEFAULT CURRENT_TIMESTAMP,
    action          TEXT,
    entity_type     TEXT,
    entity_id       INTEGER,
    note            TEXT,
    payload         TEXT,
    undoes_event_id INTEGER REFERENCES user_actions(id),
    client_id       TEXT,
    client_op_id    TEXT
);
CREATE TABLE idempotent_responses (
    client_op_id  TEXT PRIMARY KEY,
    status_code   INTEGER NOT NULL,
    body          TEXT NOT NULL,
    created_at    DATETIME DEFAULT CURRENT_TIMESTAMP
);
CREATE TABLE schema_migrations (
    name        TEXT PRIMARY KEY,
    applied_at  TEXT
);
CREATE TABLE speculative_peak_intents (
    index_id       INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
    ratio_position INTEGER NOT NULL,
    q              REAL    NOT NULL,
    PRIMARY KEY (index_id, ratio_position)
);
CREATE TABLE sample_bridges (
    sample_id INTEGER PRIMARY KEY REFERENCES samples(id) ON DELETE CASCADE,
    group_id  INTEGER NOT NULL
);
CREATE TABLE sample_labels (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id INTEGER NOT NULL REFERENCES samples(id) ON DELETE CASCADE,
    label     TEXT NOT NULL,
    state     TEXT NOT NULL CHECK (state IN ('suggested', 'confirmed', 'removed')),
    source    TEXT NOT NULL DEFAULT 'user'
);
CREATE TABLE label_colors (
    label TEXT NOT NULL,
    color TEXT NOT NULL
);
CREATE INDEX idx_auto_peaks_exposure
    ON auto_peaks(exposure_id);
CREATE INDEX idx_peak_curations_exposure
    ON peak_curations(exposure_id);
CREATE UNIQUE INDEX idx_one_custom_group_per_exposure
    ON index_groups(exposure_id) WHERE kind = 'custom';
CREATE INDEX idx_events_by_exposure
    ON user_actions(entity_type, entity_id, id);
CREATE INDEX idx_idempotent_responses_created
    ON idempotent_responses(created_at);
CREATE INDEX idx_events_by_client_op_id ON user_actions(client_op_id) WHERE client_op_id IS NOT NULL;
CREATE UNIQUE INDEX idx_events_unique_op
ON user_actions(client_op_id, action, entity_id)
WHERE client_op_id IS NOT NULL;
CREATE INDEX loads_experiment_idx ON loads(experiment_id);
CREATE UNIQUE INDEX exposures_unique_filename ON exposures(experiment_id, filename);
CREATE INDEX exposures_experiment_idx ON exposures(experiment_id);
CREATE INDEX exposures_load_idx ON exposures(load_id);
CREATE INDEX samples_experiment_idx ON samples(experiment_id);
CREATE INDEX samples_load_slot_idx ON samples(load_id, slot_index);
CREATE INDEX idx_sample_bridges_group ON sample_bridges(group_id);
CREATE UNIQUE INDEX idx_sample_labels_unique ON sample_labels(sample_id, label COLLATE NOCASE);
CREATE UNIQUE INDEX idx_label_colors_label ON label_colors(label COLLATE NOCASE);
INSERT INTO schema_migrations (name, applied_at) VALUES ('assignments_v1', '2026-11-21T23:56:51.647Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('loads_table_v1', '2026-11-21T23:56:51.648Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('exposures_experiment_id_v1', '2026-11-21T23:56:51.950Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('experiments_geometry_v1', '2026-11-21T23:56:51.950Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('samples_name_collapse_v1', '2026-11-21T23:56:51.962Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('samples_merged_into_v1', '2026-11-21T23:56:51.962Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('experiments_config_cols_v1', '2026-11-21T23:56:51.963Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('load_sessions_backfill_v1', '2026-11-21T23:56:52.061Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('drop_scan_signature_v1', '2026-11-21T23:56:52.061Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('exposures_beam_monitors_v1', '2026-11-21T23:56:52.061Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('exposures_ingest_state_v1', '2026-11-21T23:56:52.061Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('experiments_ingest_error_v1', '2026-11-21T23:56:52.062Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('speculative_peak_durability', '2026-11-21T23:56:52.120Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('hex_sqrt11_removal_v1', '2026-11-21T23:56:52.147Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('index_peaks_ratio_position_int_v1', '2026-11-21T23:56:52.192Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('sample_calls_v1', '2026-11-21T23:56:52.371Z');
INSERT INTO schema_migrations (name, applied_at) VALUES ('label_suggestions_v1', '2026-11-21T23:56:52.408Z');
