using Test, SQLite, DBInterface, Tables, JSON3
using HimalayaUI: create_schema!, migrate_schema!, migrate_comparisons_to_series!,
                  MIGRATION_COMPARISONS_TO_SERIES

@testset "migrate_comparisons_to_series!" begin
    @testset "empty corpus: sentinel written, no series rows, idempotent" begin
        db = SQLite.DB()
        create_schema!(db)
        migrate_schema!(db)   # runs the migration once as part of the sequence

        # Sentinel row present after the full migrate_schema! run.
        sent = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES]))
        @test length(sent) == 1

        # No comparisons → no series.
        n_series = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM series"))[1].n
        @test n_series == 0

        # Calling the migration again is a gated no-op (does not throw, no dup sentinel).
        migrate_comparisons_to_series!(db)
        sent2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES]))[1].n
        @test sent2 == 1
        close(db)
    end

    @testset "single comparison → committed series with recipe + plate" begin
        db = SQLite.DB()
        create_schema!(db)
        migrate_schema!(db)   # creates all tables incl. series* + sentinel

        # Seed minimal experiment/sample/exposure graph so members reference real
        # exposures with real sample_ids (recipe derivation needs exposures.sample_id).
        DBInterface.execute(db, """INSERT INTO experiments
            (id, name, path, data_dir, analysis_dir)
            VALUES (1, 'exp', '/x', '/x/data', '/x/analysis')""")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id) VALUES (10, 1), (11, 1)")
        DBInterface.execute(db,
            "INSERT INTO exposures (id, sample_id) VALUES (100, 10), (101, 11)")

        # A comparison with a known created_by/created_at and two members.
        DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (5, 'alice')")
        DBInterface.execute(db, """INSERT INTO comparisons
            (id, title, description, content_hash, created_by, created_at, updated_at)
            VALUES (1, 'Cmp A', 'desc', 'oldhash', 5, '2026-01-01T00:00:00.000Z',
                    '2026-01-02T00:00:00.000Z')""")
        snap = JSON3.write(Dict(:effective_peaks => [], :confirmed_index => nothing,
                                :analysis_inputs_hash => nothing))
        # `band_height`/`y_offset`/`normalization` are NOT NULL DEFAULT in the
        # comparison_members schema (db.jl: REAL/REAL/TEXT with defaults
        # 1.0/0/'none'), so the plate-payload helper's Float64(...)/String(...)
        # coercions always get a real value even when the fixture omits them. We
        # set them EXPLICITLY here so a future schema change that drops the
        # defaults surfaces as a fixture failure, not a misdiagnosed migration bug.
        DBInterface.execute(db, """INSERT INTO comparison_members
            (comparison_id, exposure_id, display_order, band_height, y_offset,
             normalization, snapshot, created_at)
            VALUES (1, 100, 0, 1.0, 0.0, 'none', ?, '2026-01-01T00:00:00.000Z'),
                   (1, 101, 1, 1.0, 0.0, 'none', ?, '2026-01-01T00:00:00.000Z')""",
            [snap, snap])

        # IMPORTANT: clear the sentinel so the migration re-runs against the now-
        # populated corpus (migrate_schema! already wrote it on the empty DB).
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES])

        migrate_comparisons_to_series!(db)

        # Exactly one series, committed, content_hash set (plate folded).
        srows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, title, description, state, content_hash, created_by, created_at FROM series"))
        @test length(srows) == 1
        s = srows[1]
        @test String(s.title) == "Cmp A"
        @test String(s.state) == "committed"
        @test !ismissing(s.content_hash)            # plate committed ⇒ hash set
        @test Int(s.created_by) == 5                # carried via user_actions.user_id
        # NOTE: do NOT assert s.created_at == the comparison's created_at — the
        # dispatcher stamps the series row's created_at with comparison_now_iso()
        # (migration time). The HISTORICAL timestamp lives on the event row's
        # `timestamp` column (asserted below), not the series view row.

        sid = Int(s.id)

        # Two members, ids minted fresh, display_order preserved.
        mrows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, display_order FROM series_members WHERE series_id = ? ORDER BY display_order", [sid]))
        @test length(mrows) == 2
        @test Int(mrows[1].exposure_id) == 100
        @test Int(mrows[2].exposure_id) == 101

        # Recipe: two distinct samples (10, 11) at positions 0, 1.
        rrows = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, position FROM series_samples WHERE series_id = ? ORDER BY position", [sid]))
        @test [Int(r.sample_id) for r in rrows] == [10, 11]
        @test [Int(r.position) for r in rrows] == [0, 1]

        # Two synthesized user_actions rows, entity_type='series', NULL client ids,
        # user_id carried, in created/committed order.
        ev = Tables.rowtable(DBInterface.execute(db,
            """SELECT action, user_id, client_op_id, client_id, timestamp FROM user_actions
               WHERE entity_type = 'series' AND entity_id = ? ORDER BY id""", [sid]))
        @test [String(e.action) for e in ev] == ["series_created", "series_plate_committed"]
        @test all(e -> Int(e.user_id) == 5, ev)
        @test all(e -> ismissing(e.client_op_id), ev)
        @test all(e -> ismissing(e.client_id), ev)
        # Historical timestamp carried onto the EVENT row (not created_at; column
        # is named `timestamp`).
        @test all(e -> String(e.timestamp) == "2026-01-01T00:00:00.000Z", ev)

        # No idempotent_responses rows written (migration is not a client op).
        nresp = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM idempotent_responses"))[1].n
        @test nresp == 0

        close(db)
    end

    @testset "messages and pins copy forward" begin
        db = SQLite.DB()
        create_schema!(db)
        migrate_schema!(db)

        DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (5, 'alice'), (6, 'bob')")
        DBInterface.execute(db, """INSERT INTO comparisons
            (id, title, created_by, created_at, updated_at)
            VALUES (1, 'Cmp A', 5, '2026-01-01T00:00:00.000Z', '2026-01-01T00:00:00.000Z')""")
        DBInterface.execute(db, """INSERT INTO comparison_messages
            (comparison_id, author_id, body, created_at)
            VALUES (1, 5, 'hello', '2026-01-03T00:00:00.000Z'),
                   (1, 6, 'world', '2026-01-04T00:00:00.000Z')""")
        DBInterface.execute(db, """INSERT INTO comparison_pins (user_id, comparison_id, pinned_at)
            VALUES (5, 1, '2026-01-05T00:00:00.000Z')""")

        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)

        sid = Int(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM series"))[1].id)

        msgs = Tables.rowtable(DBInterface.execute(db,
            "SELECT author_id, body, created_at FROM series_messages WHERE series_id = ? ORDER BY created_at", [sid]))
        @test [String(m.body) for m in msgs] == ["hello", "world"]
        @test String(msgs[1].created_at) == "2026-01-03T00:00:00.000Z"   # carried, not now()
        @test Int(msgs[2].author_id) == 6

        pins = Tables.rowtable(DBInterface.execute(db,
            "SELECT user_id, pinned_at FROM series_pins WHERE series_id = ?", [sid]))
        @test length(pins) == 1
        @test Int(pins[1].user_id) == 5
        @test String(pins[1].pinned_at) == "2026-01-05T00:00:00.000Z"    # carried
        close(db)
    end

    @testset "migrated series survives a rebuild_views_from_log! round-trip" begin
        db = SQLite.DB()
        create_schema!(db)
        migrate_schema!(db)
        DBInterface.execute(db, """INSERT INTO experiments
            (id, name, path, data_dir, analysis_dir)
            VALUES (1, 'exp', '/x', '/x/data', '/x/analysis')""")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id) VALUES (10, 1)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id) VALUES (100, 10)")
        DBInterface.execute(db, """INSERT INTO comparisons
            (id, title, created_by, created_at, updated_at)
            VALUES (1, 'Cmp A', NULL, '2026-01-01T00:00:00.000Z', '2026-01-01T00:00:00.000Z')""")
        snap = JSON3.write(Dict(:effective_peaks => [], :confirmed_index => nothing,
                                :analysis_inputs_hash => nothing))
        DBInterface.execute(db, """INSERT INTO comparison_members
            (comparison_id, exposure_id, display_order, band_height, y_offset,
             normalization, snapshot, created_at)
            VALUES (1, 100, 0, 1.0, 0.0, 'none', ?, '2026-01-01T00:00:00.000Z')""", [snap])
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)

        sid = Int(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM series"))[1].id)
        hash_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT content_hash FROM series WHERE id = ?", [sid]))[1].content_hash
        members_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, display_order FROM series_members WHERE series_id = ? ORDER BY display_order", [sid]))
        recipe_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, position FROM series_samples WHERE series_id = ? ORDER BY position", [sid]))

        # Empty the view rows, then re-fold the synthesized log (exercises the
        # MIGRATED path, not only native series — master plan §6.3).
        DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])
        DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
        DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
        HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")

        hash_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT content_hash FROM series WHERE id = ?", [sid]))[1].content_hash
        members_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, display_order FROM series_members WHERE series_id = ? ORDER BY display_order", [sid]))
        recipe_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, position FROM series_samples WHERE series_id = ? ORDER BY position", [sid]))

        @test String(hash_after) == String(hash_before)
        @test [Int(m.exposure_id) for m in members_after] == [Int(m.exposure_id) for m in members_before]
        @test [Int(r.sample_id) for r in recipe_after] == [Int(r.sample_id) for r in recipe_before]
        # NOTE: do NOT assert series_messages/series_pins survive this round-trip.
        # They are NOT event-sourced (copied by raw INSERT in Task 3); the
        # DELETE FROM series CASCADE drops them and rebuild_views_from_log!
        # (which folds only series_created/series_plate_committed) does not
        # restore them. That is correct and matches native series.
        close(db)
    end

    @testset "second migrate_schema! run is a gated no-op (no duplicate series)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "himalaya.db")
            db = HimalayaUI.open_db(path)
            DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (5, 'alice')")
            DBInterface.execute(db, """INSERT INTO comparisons
                (id, title, created_by, created_at, updated_at)
                VALUES (1, 'Cmp A', 5, '2026-01-01T00:00:00.000Z', '2026-01-01T00:00:00.000Z')""")
            # Simulate the real cutover: the comparison was created AFTER the
            # empty-DB migration already wrote the sentinel. Clear it so the
            # first "real" copy runs, then re-open to prove idempotency.
            DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
                [MIGRATION_COMPARISONS_TO_SERIES])
            HimalayaUI.migrate_comparisons_to_series!(db)
            n1 = Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM series"))[1].n
            close(db)

            # Re-open → migrate_schema! runs migrate_comparisons_to_series! again;
            # sentinel is now present, so it must be a no-op (no second series).
            db2 = HimalayaUI.open_db(path)
            n2 = Tables.rowtable(DBInterface.execute(db2, "SELECT COUNT(*) AS n FROM series"))[1].n
            @test n1 == 1
            @test n2 == 1   # no duplicate copy on re-open
            close(db2)
        end
    end

    @testset "orphan member (NULL exposure_id) contributes no recipe row but is copied to plate" begin
        db = SQLite.DB(); create_schema!(db); migrate_schema!(db)
        DBInterface.execute(db, """INSERT INTO experiments
            (id, name, path, data_dir, analysis_dir)
            VALUES (1, 'e', '/x', '/x/data', '/x/analysis')""")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id) VALUES (10,1)")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id) VALUES (100,10)")
        DBInterface.execute(db, """INSERT INTO comparisons (id,title,created_at,updated_at)
            VALUES (1,'C','2026-01-01T00:00:00.000Z','2026-01-01T00:00:00.000Z')""")
        snap = JSON3.write(Dict(:effective_peaks=>[], :confirmed_index=>nothing, :analysis_inputs_hash=>nothing))
        # member 0 has a real exposure; member 1 is an orphan (NULL exposure_id).
        DBInterface.execute(db, """INSERT INTO comparison_members
            (comparison_id, exposure_id, display_order, band_height, y_offset, normalization, snapshot, created_at)
            VALUES (1, 100, 0, 1.0, 0.0, 'none', ?, '2026-01-01T00:00:00.000Z'),
                   (1, NULL, 1, 1.0, 0.0, 'none', ?, '2026-01-01T00:00:00.000Z')""", [snap, snap])
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)
        sid = Int(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM series"))[1].id)
        # Plate has both members (orphan included).
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM series_members WHERE series_id = ?", [sid]))[1].n == 2
        # Recipe has one entry (only the non-orphan exposure's sample).
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM series_samples WHERE series_id = ?", [sid]))[1].n == 1
        close(db)
    end

    @testset "forked comparison → migrated series has NULL fork lineage (P1: id-space mismatch)" begin
        db = SQLite.DB(); create_schema!(db); migrate_schema!(db)
        # Parent + child comparison; the child forks the parent (comparisons id-space).
        DBInterface.execute(db, """INSERT INTO comparisons (id,title,created_at,updated_at)
            VALUES (1,'Parent','2026-01-01T00:00:00.000Z','2026-01-01T00:00:00.000Z')""")
        DBInterface.execute(db, """INSERT INTO comparisons
            (id,title,forked_from_id,forked_at_hash,created_at,updated_at)
            VALUES (2,'Child',1,'parenthash','2026-01-02T00:00:00.000Z','2026-01-02T00:00:00.000Z')""")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)
        # Both comparisons copied; the forked child's series row NULLs its lineage
        # (comparisons.forked_from_id is a comparisons-id-space ref; copying it
        # straight through would dangle against series(id)). forked_at_hash NULL too.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT title, forked_from_id, forked_at_hash FROM series ORDER BY id"))
        @test length(rows) == 2
        @test all(r -> ismissing(r.forked_from_id), rows)
        @test all(r -> ismissing(r.forked_at_hash), rows)
        # And the copy leaves no dangling FK: foreign_key_check is clean.
        fk = Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check(series)"))
        @test isempty(fk)
        close(db)
    end

    @testset "comparison with zero members → empty recipe + empty plate, committed" begin
        db = SQLite.DB(); create_schema!(db); migrate_schema!(db)
        DBInterface.execute(db, """INSERT INTO comparisons (id,title,created_at,updated_at)
            VALUES (1,'Empty','2026-01-01T00:00:00.000Z','2026-01-01T00:00:00.000Z')""")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)
        sid = Int(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM series"))[1].id)
        @test Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM series_members WHERE series_id=?", [sid]))[1].n == 0
        @test Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM series_samples WHERE series_id=?", [sid]))[1].n == 0
        # A committed plate with zero members still gets state='committed' + a hash.
        st = Tables.rowtable(DBInterface.execute(db, "SELECT state FROM series WHERE id=?", [sid]))[1].state
        @test String(st) == "committed"
        close(db)
    end

    @testset "multiple comparisons copy in id order, ids independent" begin
        db = SQLite.DB(); create_schema!(db); migrate_schema!(db)
        for cid in (1, 2, 3)
            DBInterface.execute(db, """INSERT INTO comparisons (id,title,created_at,updated_at)
                VALUES (?, ?, '2026-01-01T00:00:00.000Z','2026-01-01T00:00:00.000Z')""",
                [cid, "C$cid"])
        end
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [MIGRATION_COMPARISONS_TO_SERIES])
        migrate_comparisons_to_series!(db)
        titles = [String(r.title) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT title FROM series ORDER BY id"))]
        @test titles == ["C1", "C2", "C3"]
        close(db)
    end
end
