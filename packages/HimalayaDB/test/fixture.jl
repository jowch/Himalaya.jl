using SQLite, DBInterface, Tables
import HimalayaUI

"""
    build_fixture(path, analysis_dir) -> NamedTuple

Create a schema-correct DB at `path` (via HimalayaUI.open_db) populated with a
deterministic sample: 3 auto peaks (q = 0.10, 0.1414, 0.1732), one `exclude`
curation on the middle peak, one `add` curation at q = 0.20, one candidate
`Pn3m` index supported by two auto peaks, and a confirmed custom index group.
Also creates a second, UNCURATED exposure (`exposure2_id`): an `active` auto
group with a member but no custom group — used to prove `confirmed_indices`
filters on `kind='custom'`, not `active=1`.
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

        rg = DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'custom', 1)",
            [exposure_id])
        group_id = Int(DBInterface.lastrowid(rg))
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group_id, index_id])

        # Second, UNCURATED exposure: an active auto group with a member, but NO
        # custom group. confirmed_indices must return EMPTY here (the curator
        # committed nothing) — this pins the filter to kind='custom', not active=1.
        exposure2_id = HimalayaUI.create_exposure!(db; experiment_id=experiment_id,
            sample_id=sample_id, filename="s2", kind="file", status="accepted")
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, intensity, prominence, sharpness, findpeaks_index) VALUES (?,?,?,?,?,?)",
            [exposure2_id, 0.11, 100.0, 5.0, 1.0, 12])
        ri2 = DBInterface.execute(db,
            "INSERT INTO indices (exposure_id, phase, basis, score, kind, status) VALUES (?, 'Im3m', ?, ?, 'auto', 'candidate')",
            [exposure2_id, 0.11, 0.5])
        index2_id = Int(DBInterface.lastrowid(ri2))
        rg2 = DBInterface.execute(db,
            "INSERT INTO index_groups (exposure_id, kind, active) VALUES (?, 'auto', 1)",
            [exposure2_id])
        group2_id = Int(DBInterface.lastrowid(rg2))
        DBInterface.execute(db,
            "INSERT INTO index_group_members (group_id, index_id) VALUES (?, ?)",
            [group2_id, index2_id])

        return (; experiment_id, sample_id, exposure_id, exposure2_id,
                  index_id, group_id, auto_peak_ids)
    finally
        close(db)
    end
end
