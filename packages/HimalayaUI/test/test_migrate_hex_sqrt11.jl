using Test, SQLite, DBInterface, Tables
using Himalaya

# #304: √11 left phaseratios(Hexagonal), renumbering every position past 5.
# Production data does NOT exercise the renumber (all intents sit at positions
# 1–3), so these fixtures are the only thing that ever drives it.

function _hex_fixture()
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="H1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")
    (; db, experiment_id = exp_id, exposure_id = e_id)
end

# open_db already ran the one-shot migration on the fresh DB; drop the sentinel
# so it can be driven again against hand-seeded rows.
_rearm!(db) = DBInterface.execute(db,
    "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_HEX_SQRT11])

function _mkindex!(db, exposure_id, phase)
    res = DBInterface.execute(db, """
        INSERT INTO indices (exposure_id, phase, basis, status, kind)
        VALUES (?, ?, 0.1, 'candidate', 'speculative')""", [exposure_id, phase])
    Int(DBInterface.lastrowid(res))
end

_intents(db, ix_id) = [(Int(r.ratio_position), round(Float64(r.q), digits = 4))
    for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT ratio_position, q FROM speculative_peak_intents WHERE index_id = ? ORDER BY ratio_position",
        [ix_id]))]

function _seed_intents!(db, ix_id, positions)
    for p in positions
        DBInterface.execute(db,
            "INSERT INTO speculative_peak_intents (index_id, ratio_position, q) VALUES (?,?,?)",
            [ix_id, p, 0.01 * p])
    end
end

@testset "hex √11 migration: position 6 orphaned, tail renumbered, q carried" begin
    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    # 7 AND 8 on the same index is the case the +1000 two-pass exists for: a
    # naive single decrement hits 8→7 while 7 still exists and aborts on the
    # (index_id, ratio_position) PK — which would roll back the sentinel and
    # make every subsequent open_db throw.
    _seed_intents!(fx.db, hex, [1, 5, 6, 7, 8, 9])

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)

    # 6 dropped (it claimed a √11 that cannot exist); 7/8/9 → 6/7/8 carrying
    # their q, so each intent keeps the radicand it was created for.
    @test _intents(fx.db, hex) == [(1, 0.01), (5, 0.05), (6, 0.07), (7, 0.08), (8, 0.09)]
end

@testset "hex √11 migration: other phases untouched" begin
    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    sq  = _mkindex!(fx.db, fx.exposure_id, "Square")
    _seed_intents!(fx.db, hex, [6, 7])
    _seed_intents!(fx.db, sq,  [6, 7])

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)

    @test _intents(fx.db, hex) == [(6, 0.07)]
    @test _intents(fx.db, sq)  == [(6, 0.06), (7, 0.07)]   # Square is 12 entries, unchanged
end

@testset "hex √11 migration: a qualified phase spelling is not skipped" begin
    # Nothing in the schema constrains `indices.phase`, and resolve_phase strips
    # a `Himalaya.` prefix because the qualified form is a known hazard. A row
    # spelled that way must not be silently passed over while the sentinel burns.
    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Himalaya.Hexagonal")
    _seed_intents!(fx.db, hex, [6, 7])

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)

    @test _intents(fx.db, hex) == [(6, 0.07)]
end

@testset "hex √11 migration: index_peaks deleted for Hexagonal only" begin
    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    sq  = _mkindex!(fx.db, fx.exposure_id, "Square")
    pid = Int(DBInterface.lastrowid(DBInterface.execute(fx.db,
        "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.1, 1.0)",
        [fx.exposure_id])))
    for ix in (hex, sq)
        DBInterface.execute(fx.db, """
            INSERT INTO index_peaks (index_id, peak_id, peak_kind, ratio_position, residual)
            VALUES (?, ?, 'auto', 1, 0.0)""", [ix, pid])
    end

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)

    n(ix) = first(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT count(*) AS n FROM index_peaks WHERE index_id = ?", [ix]))).n
    @test n(hex) == 0
    @test n(sq)  == 1
end

@testset "hex √11 migration: reopens the memoization gate, skipping zero-index exposures" begin
    fx = _hex_fixture()
    # Exposure WITH an index — must be invalidated so analyze re-indexes it.
    _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'abc' WHERE id = ?", [fx.exposure_id])
    # Exposure with NO index — provably unaffected by the series change, so its
    # hash must survive (this is what keeps the reset off 3148 prod exposures).
    bare = HimalayaUI.create_exposure!(fx.db;
        experiment_id=fx.experiment_id,
        sample_id=first(Tables.rowtable(DBInterface.execute(fx.db,
            "SELECT sample_id FROM exposures WHERE id = ?", [fx.exposure_id]))).sample_id,
        filename="y")
    DBInterface.execute(fx.db,
        "UPDATE exposures SET analysis_inputs_hash = 'def' WHERE id = ?", [bare])

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)

    h(id) = first(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT analysis_inputs_hash AS x FROM exposures WHERE id = ?", [id]))).x
    @test ismissing(h(fx.exposure_id))
    @test h(bare) == "def"
end

@testset "hex √11 migration: sentinel makes it one-shot" begin
    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    _seed_intents!(fx.db, hex, [7, 8])

    _rearm!(fx.db)
    HimalayaUI.migrate_hex_sqrt11!(fx.db)
    once = _intents(fx.db, hex)
    @test once == [(6, 0.07), (7, 0.08)]

    # Second call without re-arming must be a no-op. The renumber is NOT
    # idempotent — a second pass would shift these to 5 and 6.
    HimalayaUI.migrate_hex_sqrt11!(fx.db)
    @test _intents(fx.db, hex) == once
end

@testset "hex √11 migration: defers (no sentinel) against a core that still has √11" begin
    # The guard that stops a stale manifest from renumbering intents into a
    # scheme the loaded core doesn't use. Can't load an old Himalaya here, so
    # drive the predicate the guard reads.
    @test !(11 in round.(Int, Himalaya.phaseratios(Himalaya.Hexagonal) .^ 2))

    fx  = _hex_fixture()
    hex = _mkindex!(fx.db, fx.exposure_id, "Hexagonal")
    _seed_intents!(fx.db, hex, [7])
    _rearm!(fx.db)

    # With the real (fixed) core the migration proceeds and records the sentinel.
    HimalayaUI.migrate_hex_sqrt11!(fx.db)
    @test _intents(fx.db, hex) == [(6, 0.07)]
    @test !isempty(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT 1 FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_HEX_SQRT11])))

    # The deferral contract that matters: bailing must NOT write the sentinel,
    # so a corrected environment still applies the migration later.
    _rearm!(fx.db)
    @test isempty(Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT 1 FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_HEX_SQRT11])))
end
