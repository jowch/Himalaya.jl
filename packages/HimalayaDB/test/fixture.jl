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
The experiment's one setup carries geometry and a `{name}_tot.dat` integration
pattern. A second sample (`sample2_id`, a buffer) is sample 1's subtraction
background and shares a bridge with it; sample 1 has an `indexed` sample call
and confirmed, suggested and removed labels. A third sample (`merged_id`) is
merged into sample 1, carries a label, is in the same bridge and is sample 2's
background; every reader must leave it out, as HimalayaUI does.
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

        # Near-boundary q (not an exact match to the auto peak's 0.1414): within
        # the auto peak's MAX(1e-6, ABS(q)*0.001) ≈ 0.0001414 tolerance band, so
        # this exercises curated_peaks' tolerance formula rather than equality.
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)",
            [exposure_id, 0.14149])
        radd = DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)", [exposure_id, 0.20])
        add_curation_id = Int(DBInterface.lastrowid(radd))

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
        # third supporting peak, backed by the add-curation (not an auto peak):
        # exercises reconstruct_index's peak_curations join and its sharpness=0
        # fidelity gap (see reconstruct.jl docstring).
        DBInterface.execute(db,
            "INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position) VALUES (?, ?, 'curation', ?)",
            [index_id, add_curation_id, 3])

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

        setup_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM setups WHERE experiment_id = ?", [experiment_id]))).id)
        DBInterface.execute(db,
            "UPDATE setups SET energy_kev = 11.0, flight_path_m = 1.1, integration_pattern = '{name}_tot.dat' WHERE id = ?",
            [setup_id])

        sample2_id = HimalayaUI.create_sample!(db; experiment_id=experiment_id, name="Buffer 1")
        DBInterface.execute(db,
            "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'indexed')", [sample_id])
        for (sid, label, state, source) in [(sample2_id, "Buffer", "confirmed", "user"),
                                            (sample_id, "Calibrant", "suggested", "system"),
                                            (sample_id, "Empty capillary", "removed", "user"),
                                            (sample_id, "Hot", "confirmed", "user")]
            DBInterface.execute(db,
                "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, ?, ?, ?)",
                [sid, label, state, source])
        end
        DBInterface.execute(db,
            "INSERT INTO sample_subtractions (sample_id, ref_sample_id) VALUES (?, ?)", [sample_id, sample2_id])
        for sid in (sample_id, sample2_id)
            DBInterface.execute(db,
                "INSERT INTO sample_bridges (sample_id, group_id) VALUES (?, 1)", [sid])
        end

        merged_id = HimalayaUI.create_sample!(db; experiment_id=experiment_id, name="Merged away")
        DBInterface.execute(db, "UPDATE samples SET merged_into_id = ? WHERE id = ?", [sample_id, merged_id])
        DBInterface.execute(db,
            "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'Ghost', 'confirmed', 'user')", [merged_id])
        DBInterface.execute(db, "INSERT INTO sample_bridges (sample_id, group_id) VALUES (?, 1)", [merged_id])
        DBInterface.execute(db,
            "INSERT INTO sample_subtractions (sample_id, ref_sample_id) VALUES (?, ?)", [sample2_id, merged_id])

        return (; experiment_id, setup_id, sample_id, sample2_id, merged_id, exposure_id, exposure2_id,
                  index_id, index2_id, auto_peak_ids)
    finally
        close(db)
    end
end
