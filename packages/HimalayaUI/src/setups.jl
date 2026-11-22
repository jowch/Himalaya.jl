using DBInterface, Tables, SQLite

const SETUP_CONFIG_COLUMNS = (
    "data_dir", "analysis_dir",
    "image_pattern", "metadata_pattern", "integration_pattern",
    "energy_kev", "energy_kev_source", "flight_path_m", "flight_path_m_source",
    "beam_center_x", "beam_center_x_source", "beam_center_y", "beam_center_y_source",
    "pixel_size_um", "pixel_size_um_source", "q_units", "q_units_source",
)

const MIGRATION_SETUP_IDS = "setup_ids_v1"

const RETIRED_EXPERIMENT_COLUMNS = (
    SETUP_CONFIG_COLUMNS...,
    "manifest_path", "config", "experiment_type", "last_scan_tier", "consecutive_empty_ticks",
)

function setup_root_dir(data_dir::AbstractString)::String
    d = String(rstrip(String(data_dir), '/'))
    basename(d) == "data" ? dirname(d) : d
end

setup_folder(root_dir::AbstractString) = String(basename(rstrip(String(root_dir), '/')))

function create_setup!(db::SQLite.DB;
        experiment_id::Integer,
        data_dir::AbstractString,
        analysis_dir::AbstractString,
        root_dir::AbstractString = setup_root_dir(data_dir),
        folder::AbstractString = setup_folder(root_dir),
        display_name = nothing,
        position::Integer = 0,
        included::Bool = true,
        image_pattern = nothing, metadata_pattern = nothing, integration_pattern = nothing,
        energy_kev = nothing, energy_kev_source = "default",
        flight_path_m = nothing, flight_path_m_source = "default",
        beam_center_x = nothing, beam_center_x_source = "default",
        beam_center_y = nothing, beam_center_y_source = "default",
        pixel_size_um = nothing, pixel_size_um_source = "default",
        q_units = nothing, q_units_source = "default")
    res = DBInterface.execute(db, """
        INSERT INTO setups
            (experiment_id, folder, display_name, position, included, root_dir, data_dir, analysis_dir,
             image_pattern, metadata_pattern, integration_pattern,
             energy_kev, energy_kev_source, flight_path_m, flight_path_m_source,
             beam_center_x, beam_center_x_source, beam_center_y, beam_center_y_source,
             pixel_size_um, pixel_size_um_source, q_units, q_units_source)
        VALUES (?,?,?,?,?,?,?,?, ?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?)
    """, [experiment_id, folder, display_name, position, Int(included), String(root_dir),
          String(data_dir), String(analysis_dir),
          image_pattern, metadata_pattern, integration_pattern,
          energy_kev, energy_kev_source, flight_path_m, flight_path_m_source,
          beam_center_x, beam_center_x_source, beam_center_y, beam_center_y_source,
          pixel_size_um, pixel_size_um_source, q_units, q_units_source])
    Int(DBInterface.lastrowid(res))
end

setups_of(db::SQLite.DB, experiment_id::Integer; included_only::Bool = false) =
    Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM setups WHERE experiment_id = ?" *
        (included_only ? " AND included = 1" : "") * " ORDER BY position, id",
        [Int(experiment_id)]))

function setup_row(db::SQLite.DB, setup_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM setups WHERE id = ?", [Int(setup_id)]))
    isempty(rows) ? nothing : first(rows)
end

function default_setup_id(db::SQLite.DB, experiment_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM setups WHERE experiment_id = ? ORDER BY position, id LIMIT 1", [Int(experiment_id)]))
    isempty(rows) ? nothing : Int(first(rows).id)
end

function exposure_setup(db::SQLite.DB, exposure_id::Integer)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT st.* FROM exposures e JOIN setups st ON st.id = e.setup_id
         WHERE e.id = ?""", [Int(exposure_id)]))
    isempty(rows) ? nothing : first(rows)
end

setup_display_name(row) = coalesce(row.display_name, row.folder)

_same_path(a, b) = !ismissing(a) && rstrip(String(a), '/') == rstrip(String(b), '/')

function _copy_experiments_into_setups!(db::SQLite.DB)
    exp_cols = cols_of(db, "experiments")
    "data_dir" in exp_cols || return nothing
    shared = [c for c in SETUP_CONFIG_COLUMNS if c in exp_cols]
    sel(c) = endswith(c, "_source") ? "COALESCE(x.$c, 'default') AS $c" : "x.$c"
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT x.id, x.path, $(join((sel(c) for c in shared), ", ")) FROM experiments x
         WHERE x.data_dir IS NOT NULL
           AND NOT EXISTS (SELECT 1 FROM setups st WHERE st.experiment_id = x.id)
         ORDER BY x.id"""))
    for r in rows
        root = setup_root_dir(String(r.data_dir))
        vals = [getproperty(r, Symbol(c)) for c in shared]
        DBInterface.execute(db, """
            INSERT INTO setups (experiment_id, folder, position, included, root_dir, $(join(shared, ", ")))
            VALUES (?, ?, 0, 1, ?, $(join(fill("?", length(shared)), ", ")))""",
            [Int(r.id), setup_folder(root), root, vals...])
        _same_path(r.path, r.data_dir) &&
            DBInterface.execute(db, "UPDATE experiments SET path = ? WHERE id = ?", [root, Int(r.id)])
    end
    nothing
end

function _backfill_setup_ids!(db::SQLite.DB)
    for t in ("loads", "samples", "exposures")
        _table_exists(db, t) || continue
        cols = cols_of(db, t)
        ("setup_id" in cols && "experiment_id" in cols) || continue
        DBInterface.execute(db, """
            UPDATE $t SET setup_id = (
                SELECT st.id FROM setups st WHERE st.experiment_id = $t.experiment_id
                 ORDER BY st.position, st.id LIMIT 1)
             WHERE setup_id IS NULL AND experiment_id IS NOT NULL""")
    end
    nothing
end

function _add_setup_id_columns!(db::SQLite.DB)
    for t in ("loads", "samples", "exposures")
        _table_exists(db, t) || continue
        "setup_id" in cols_of(db, t) ||
            DBInterface.execute(db, "ALTER TABLE $t ADD COLUMN setup_id INTEGER REFERENCES setups(id)")
    end
    nothing
end

function _setup_keyed_indexes!(db::SQLite.DB)
    if _table_exists(db, "exposures") && "setup_id" in cols_of(db, "exposures")
        idx = Tables.rowtable(DBInterface.execute(db,
            "SELECT sql FROM sqlite_master WHERE type = 'index' AND name = 'exposures_unique_filename'"))
        if isempty(idx) || !occursin("setup_id", String(coalesce(first(idx).sql, "")))
            DBInterface.execute(db, "DROP INDEX IF EXISTS exposures_unique_filename")
            DBInterface.execute(db,
                "CREATE UNIQUE INDEX exposures_unique_filename ON exposures(setup_id, filename)")
        end
        DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS exposures_setup_idx ON exposures(setup_id)")
    end
    if _table_exists(db, "loads") && "setup_id" in cols_of(db, "loads")
        DBInterface.execute(db,
            "CREATE INDEX IF NOT EXISTS loads_setup_idx ON loads(setup_id, load_index)")
    end
    if _table_exists(db, "samples") && "setup_id" in cols_of(db, "samples")
        DBInterface.execute(db, "CREATE INDEX IF NOT EXISTS samples_setup_idx ON samples(setup_id)")
    end
    nothing
end

function _drop_retired_experiment_columns!(db::SQLite.DB)
    existing = cols_of(db, "experiments")
    for c in RETIRED_EXPERIMENT_COLUMNS
        c in existing && DBInterface.execute(db, "ALTER TABLE experiments DROP COLUMN $c")
    end
    nothing
end

function migrate_setups!(db::SQLite.DB; finish::Bool = true)
    _table_exists(db, "experiments") || return nothing
    SQLite.transaction(db) do
        _copy_experiments_into_setups!(db)
        if finish
            _add_setup_id_columns!(db)
            if !_migrated(db, MIGRATION_SETUP_IDS)
                _backfill_setup_ids!(db)
                _record_migration!(db, MIGRATION_SETUP_IDS)
            end
            _setup_keyed_indexes!(db)
            _drop_retired_experiment_columns!(db)
        end
    end
    nothing
end
