using SQLite, DBInterface, Tables
import HimalayaUI

"""
    build_fixture(path, analysis_dir) -> NamedTuple

Create a schema-correct DB at `path` (via HimalayaUI.open_db) populated with a
deterministic sample: 3 auto peaks (q = 0.10, 0.1414, 0.1732), one `exclude`
curation on the middle peak, one `add` curation at q = 0.20, one candidate
`Pn3m` index supported by two auto peaks, and a confirmed assignment
(`assignments.state = 'indexed'` + an `assignment_members` row) for `exposure_id`.
Also creates a second exposure (`exposure2_id`) with an assignment member row
whose `assignments.state = 'form_factor'` (NOT `'indexed'`) — used to prove
`confirmed_indices` gates on `assignments.state = 'indexed'`, not merely on
`assignment_members` row presence.
Direct INSERTs are fine here — a reader only cares that the view tables are populated.
"""
function build_fixture(path::AbstractString, analysis_dir::AbstractString)
    db = HimalayaUI.open_db(path)  # creates schema + migrations
    try
        experiment_id = HimalayaUI.create_experiment!(db; name="fixture",
            path=analysis_dir, data_dir=analysis_dir, analysis_dir=analysis_dir)
        # HimalayaUI.create_sample! takes only `name` (no `display_name`): the
        # post-redesign migration `migrate_samples_name_collapse!` unconditionally
        # renames samples.display_name -> samples.name and drops the old `name`
        # column on every fresh DB, so there is no separate display_name column
        # to populate after open_db. `name` IS the display string here.
        sample_id = HimalayaUI.create_sample!(db; experiment_id=experiment_id,
            name="Sample 1")
        # HimalayaUI.create_exposure! requires experiment_id (no default).
        exposure_id = HimalayaUI.create_exposure!(db; experiment_id=experiment_id,
            sample_id=sample_id, filename="s1", kind="file", status="accepted")

        auto_peak_ids = Int[]
        for (q, sharp, fidx) in [(0.10, 1.0, 10), (0.1414, 0.8, 20), (0.1732, 0.6, 30)]
            r = DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness, findpeaks_index) VALUES (?,?,?,?,?,?)",
                [exposure_id, q, 100.0, 5.0, sharp, fidx])
            push!(auto_peak_ids, Int(DBInterface.lastrowid(r)))
        end

        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
            [exposure_id, 0.1414])
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)",
            [exposure_id, 0.20])

        ri = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, kind, status) VALUES (?, 'Pn3m', ?, ?, 'auto', 'candidate')",
            [exposure_id, 0.10, 0.9])
        index_id = Int(DBInterface.lastrowid(ri))
        # supporting peaks: ratio positions 1 and 2 -> first two auto peaks
        for (pos, pid) in [(1, auto_peak_ids[1]), (2, auto_peak_ids[3])]
            DBInterface.execute(db,
                "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position) VALUES (?, ?, 'auto', ?)",
                [index_id, pid, pos])
        end

        DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [exposure_id])
        DBInterface.execute(db,
            "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [exposure_id, index_id])

        # Second exposure: an assignment_members row exists, but its assignments
        # row is state='form_factor' (NOT 'indexed'). confirmed_indices must return
        # EMPTY here — this pins the filter to assignments.state='indexed', not
        # merely the presence of an assignment_members row.
        exposure2_id = HimalayaUI.create_exposure!(db; experiment_id=experiment_id,
            sample_id=sample_id, filename="s2", kind="file", status="accepted")
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness, findpeaks_index) VALUES (?,?,?,?,?,?)",
            [exposure2_id, 0.11, 100.0, 5.0, 1.0, 12])
        ri2 = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, kind, status) VALUES (?, 'Im3m', ?, ?, 'auto', 'candidate')",
            [exposure2_id, 0.11, 0.5])
        index2_id = Int(DBInterface.lastrowid(ri2))
        DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')", [exposure2_id])
        DBInterface.execute(db,
            "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [exposure2_id, index2_id])

        return (; experiment_id, sample_id, exposure_id, exposure2_id,
                  index_id, auto_peak_ids)
    finally
        close(db)
    end
end
