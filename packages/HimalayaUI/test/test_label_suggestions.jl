using Test, SQLite, DBInterface, Tables
using HimalayaUI

function _seed_named_samples(db, names)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    ids = [HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
               slot_index=i, name=n) for (i, n) in enumerate(names)]
    (exp_id, ids)
end

_rows(db) = Set((Int(r.sample_id), String(r.label), String(r.state), String(r.source)) for r in
    Tables.rowtable(DBInterface.execute(db, "SELECT sample_id, label, state, source FROM sample_labels")))

@testset "label suggestions" begin
    @testset "suggests a role only when the sample name says so" begin
        db = open_prepared_clone(mktempdir())
        names = ["buffer_blank_S2411", "HEPES_buf", "PBS", "blank", "Tris-2",
                 "AgBeh_cal", "silver behenate", "LaB6",
                 "empty_capillary", "EmptyCap3",
                 "HEPES_S01", "PBS_2", "buf_3",
                 "MO_H2O_S2401", "vesicles_ULV_S2412", "JC001", "tristearin", "pbsx",
                 "MO_HEPES_S01", "DOPC_in_PBS", "lipid_10mM_Tris"]
        (exp_id, ids) = _seed_named_samples(db, [HimalayaUI._auto_name(n, 3, i) for (i, n) in enumerate(names)])
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        got = Dict(names[findfirst(==(sid), ids)] => label for (sid, label, _, _) in _rows(db))
        @test got == Dict(
            "buffer_blank_S2411" => "Buffer", "HEPES_buf" => "Buffer", "PBS" => "Buffer",
            "blank" => "Buffer", "Tris-2" => "Buffer",
            "HEPES_S01" => "Buffer", "PBS_2" => "Buffer", "buf_3" => "Buffer",
            "AgBeh_cal" => "Calibrant", "silver behenate" => "Calibrant", "LaB6" => "Calibrant",
            "empty_capillary" => "Empty capillary", "EmptyCap3" => "Empty capillary")
        @test all(state == "suggested" && source == "ingest" for (_, _, state, source) in _rows(db))
    end

    @testset "a long digit run in a solvent sample name matches quickly and never throws" begin
        for n in (16, 22, 40, 200)
            @test HimalayaUI.suggested_role("PBS_" * "1"^n * "x") === nothing
            @test HimalayaUI.suggested_role(HimalayaUI._auto_name("PBS_" * "1"^n, 1, 1)) == "Buffer"
        end
    end

    @testset "never touches a sample that already carries a role, in any state" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_named_samples(db, ["buffer_A", "buffer_B", "buffer_C"])
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'Calibrant', 'confirmed', 'user')", [a])
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'buffer', 'removed', 'user')", [b])
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'retake', 'confirmed', 'user')", [c])
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        @test _rows(db) == Set([
            (a, "Calibrant", "confirmed", "user"),
            (b, "buffer", "removed", "user"),
            (c, "retake", "confirmed", "user"),
            (c, "Buffer", "suggested", "ingest"),
        ])
    end

    @testset "is idempotent" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_named_samples(db, ["buffer_A", "buffer_B"])
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        @test _rows(db) == Set([(a, "Buffer", "suggested", "ingest"), (b, "Buffer", "suggested", "ingest")])
    end

    @testset "an ingest pass suggests roles for the experiment" begin
        root = mktempdir()
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_named_samples(db, ["buffer_blank"])
        mkpath(joinpath(root, "data")); mkpath(joinpath(root, "analysis"))
        DBInterface.execute(db, "UPDATE setups SET data_dir = ?, analysis_dir = ? WHERE experiment_id = ?",
            [joinpath(root, "data"), joinpath(root, "analysis"), exp_id])
        @test HimalayaUI.run_ingest!(db, exp_id) === :ok
        HimalayaUI.wait_prewarm(exp_id)
        @test _rows(db) == Set([(a, "Buffer", "suggested", "ingest")])
    end

    @testset "a failed suggestion step leaves the ingest pass complete" begin
        root = mktempdir()
        db = open_prepared_clone(mktempdir())
        (exp_id, _) = _seed_named_samples(db, ["buffer_blank"])
        mkpath(joinpath(root, "data")); mkpath(joinpath(root, "analysis"))
        DBInterface.execute(db, "UPDATE setups SET data_dir = ?, analysis_dir = ? WHERE experiment_id = ?",
            [joinpath(root, "data"), joinpath(root, "analysis"), exp_id])
        DBInterface.execute(db, "DROP TABLE sample_labels")
        @test (@test_logs (:warn, "label suggestions failed") match_mode=:any HimalayaUI.run_ingest!(db, exp_id)) === :ok
        HimalayaUI.wait_prewarm(exp_id)
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT ingest_status FROM experiments WHERE id = ?", [exp_id]))).ingest_status == "complete"
    end

    @testset "opening a database suggests roles once for experiments ingested before suggestions existed" begin
        db = open_prepared_clone(mktempdir())
        (_, (a, b)) = _seed_named_samples(db, ["buffer_A", "JC001"])
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_LABEL_SUGGESTIONS])
        HimalayaUI.migrate_schema!(db)
        @test _rows(db) == Set([(a, "Buffer", "suggested", "ingest")])
        DBInterface.execute(db, "DELETE FROM sample_labels")
        HimalayaUI.migrate_schema!(db)
        @test isempty(_rows(db))
    end

    @testset "a failed suggestion step on open never blocks the database from opening" begin
        db = open_prepared_clone(mktempdir())
        (_, (a,)) = _seed_named_samples(db, ["buffer_A"])
        (_, (b,)) = _seed_named_samples(db, ["PBS"])
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_LABEL_SUGGESTIONS])
        DBInterface.execute(db, """CREATE TRIGGER fail_one BEFORE INSERT ON sample_labels
            WHEN NEW.sample_id = $a BEGIN SELECT RAISE(ABORT, 'boom'); END""")
        @test_logs (:warn, "label suggestions failed") match_mode=:any HimalayaUI.migrate_schema!(db)
        @test _rows(db) == Set([(b, "Buffer", "suggested", "ingest")])
    end
    @testset "suggests the flight path when every frame of a sample agrees" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d, e, f)) = _seed_named_samples(db, ["MO_A", "MO_B", "MO_C", "MO_D", "MO_E", "MO_F"])
        frame(sid, m) = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = sid,
            filename = "f$(sid)_$(rand(UInt32))", pipe_length_m = m)
        frame(a, 1.7095); frame(a, 1.7095)
        frame(b, 0.7); frame(b, 1.7)
        frame(c, nothing)
        frame(d, 1.7); frame(d, nothing)
        frame(e, 1.7)
        frame(f, 1.7)
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, '1.70 m', 'removed', 'user')", [e])
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, '1 M', 'confirmed', 'user')", [f])
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        HimalayaUI.suggest_labels_for_experiment!(db, exp_id)
        @test _rows(db) == Set([
            (a, "1.7 m", "suggested", "ingest"),
            (d, "1.7 m", "suggested", "ingest"),
            (e, "1.70 m", "removed", "user"),
            (f, "1 M", "confirmed", "user"),
            (f, "1.7 m", "suggested", "ingest"),
        ])
        @test HimalayaUI.flight_path_label(0.5) == "0.5 m"
    end

end
