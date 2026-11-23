using Test, Random, SQLite, DBInterface, Tables
using HimalayaUI

# Lorentzian peaks on a power-law background with Poisson-like noise.
function _pl_trace(rng, peaks; n = 700, amp = 400.0)
    q = collect(range(0.02, 0.5; length = n))
    I = @. 3000 * (q / 0.02)^-1.8 + 40
    for (q0, a, w) in peaks
        @. I += a * amp * (w^2 / ((q - q0)^2 + w^2))
    end
    σ = sqrt.(I)
    q, I .+ σ .* randn(rng, n), σ
end

function _pl_write_dat(path, q, I, σ)
    open(path, "w") do io
        for i in eachindex(q); println(io, q[i], " ", I[i], " ", σ[i]); end
    end
end

function _pl_db(tmp)
    adir = joinpath(tmp, "analysis"); mkpath(adir)
    rng = MersenneTwister(3)
    q, I, σ = _pl_trace(rng, [(0.1, 3.0, 0.002), (0.2, 3.0, 0.002), (0.3, 3.0, 0.002)])
    _pl_write_dat(joinpath(adir, "ex_tot.dat"), q, I, σ)
    _pl_write_dat(joinpath(adir, "ff_tot.dat"), _pl_trace(rng, [])...)
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path = tmp, data_dir = joinpath(tmp, "data"), analysis_dir = adir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    ff_sample = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "F")
    e_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "ex_tot")
    ff_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = ff_sample, filename = "ff_tot")
    DBInterface.execute(db, "INSERT OR REPLACE INTO sample_calls (sample_id, state) VALUES (?, 'form_factor')", [ff_sample])
    HimalayaUI.analyze_exposure!(db, e_id, adir)
    auto = [r.q for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))]
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)", [e_id, auto[1]])
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)", [e_id, 0.4])
    HimalayaUI.record_review!(db, e_id; direct = true)
    (; db, path = db.file, exp_id, s_id, ff_sample, e_id, ff_id, auto, adir)
end

@testset "labeled_exposures" begin
    tmp = mktempdir()
    b = _pl_db(tmp)
    @test length(b.auto) >= 2
    gone = HimalayaUI.create_exposure!(b.db; experiment_id = b.exp_id, sample_id = b.s_id, filename = "gone_tot")
    orphan = HimalayaUI.create_exposure!(b.db; experiment_id = b.exp_id, sample_id = b.s_id, filename = "orph_tot")
    cp(joinpath(b.adir, "ex_tot.dat"), joinpath(b.adir, "orph_tot.dat"))
    DBInterface.execute(b.db, "UPDATE exposures SET sample_id = NULL WHERE id = ?", [orphan])
    for id in (gone, orphan)
        DBInterface.execute(b.db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.4)", [id])
        HimalayaUI.record_review!(b.db, id; direct = true)
    end

    log = IOBuffer()
    got = HimalayaUI.labeled_exposures(b.db; io = log)
    @test sort([e.id for e in got]) == sort([b.e_id, b.ff_id, orphan])
    msg = String(take!(log))
    @test occursin("skipped 1 exposures", msg) && occursin("trace file not found", msg)
    e = got[findfirst(e -> e.id == b.e_id, got)]
    @test e.labelled && !e.form_factor && e.group == b.exp_id
    @test e.excluded == [b.auto[1]]
    @test e.adds == [0.4]
    @test sort(e.kept) == b.auto[2:end]
    ff = got[findfirst(e -> e.id == b.ff_id, got)]
    @test ff.form_factor && !ff.labelled && isempty(ff.adds)
    @test !any(e -> e.id == b.ff_id, HimalayaUI.labeled_exposures(b.db; include_form_factor = false, io = devnull))

    log = IOBuffer()
    by_sample = HimalayaUI.labeled_exposures(b.db; group_by = :sample, io = log)
    @test sort([e.id for e in by_sample]) == sort([b.e_id, b.ff_id])
    @test Dict(e.id => e.group for e in by_sample) == Dict(b.e_id => b.s_id, b.ff_id => b.ff_sample)
    msg = String(take!(log))
    @test occursin("skipped 2 exposures", msg) && occursin("no experiment or sample id", msg)
    @test_throws ArgumentError HimalayaUI.labeled_exposures(b.db; group_by = :beamtime)
    DBInterface.close!(b.db)
end

@testset "labeled_exposures keeps an exclude that matches no auto peak" begin
    b = _pl_db(mktempdir())
    DBInterface.execute(b.db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', 0.45)", [b.e_id])
    e = only(filter(e -> e.id == b.e_id, HimalayaUI.labeled_exposures(b.db; io = devnull)))
    @test sort(e.excluded) == sort([b.auto[1], 0.45])
    @test sort(e.kept) == b.auto[2:end]
    DBInterface.close!(b.db)
end

@testset "labeled_exposures works on a read-only handle and propagates unexpected errors" begin
    b = _pl_db(mktempdir())
    DBInterface.close!(b.db)
    ro = SQLite.DB("file:" * b.path * "?mode=ro")
    DBInterface.execute(ro, "PRAGMA query_only = ON")
    @test length(HimalayaUI.labeled_exposures(ro; io = devnull)) == 2
    DBInterface.close!(ro)

    b2 = _pl_db(mktempdir())
    DBInterface.execute(b2.db, "DROP TABLE peak_curations")
    @test_throws SQLite.SQLiteException HimalayaUI.labeled_exposures(b2.db; io = devnull)
    DBInterface.close!(b2.db)
end

_pl_set_peaks!(db, e_id, qs, found_by) = begin
    DBInterface.execute(db, "DELETE FROM auto_peaks WHERE exposure_id = ?", [e_id])
    for (q, fb) in zip(qs, found_by)
        DBInterface.execute(db, "INSERT INTO auto_peaks (exposure_id, q, found_by) VALUES (?, ?, ?)", [e_id, q, fb])
    end
end

@testset "labeled_exposures reads the latest review snapshot, not live auto peaks" begin
    b = _pl_db(mktempdir())
    DBInterface.execute(b.db, "DELETE FROM peak_curations")
    DBInterface.execute(b.db, "DELETE FROM peak_reviews")
    DBInterface.execute(b.db, "UPDATE exposures SET reviewed_at = NULL")
    @test isempty(HimalayaUI.labeled_exposures(b.db; include_form_factor = false, io = devnull))

    _pl_set_peaks!(b.db, b.e_id, [0.1, 0.2, 0.3], ["hybrid", "both", "model"])
    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "set_exposure_status",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:status => "accepted"))
    pick() = only(filter(e -> e.id == b.e_id, HimalayaUI.labeled_exposures(b.db; io = devnull)))
    @test isempty(HimalayaUI.labeled_exposures(b.db; include_form_factor = false, io = devnull))
    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "peak_excluded",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:q => 0.9))
    e = pick()
    @test e.labelled && isempty(e.adds) && e.excluded == [0.9]
    @test e.kept == [0.1, 0.2] && e.kept_model == [0.3]

    _pl_set_peaks!(b.db, b.e_id, [0.5], ["hybrid"])
    @test pick().kept == [0.1, 0.2]

    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "peak_excluded",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:q => 0.1))
    e = pick()
    @test e.kept == [0.5] && isempty(e.kept_model) && sort(e.excluded) == [0.1, 0.9]

    _pl_set_peaks!(b.db, b.e_id, [0.1, 0.3], ["hybrid", "model"])
    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "peak_excluded",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:q => 0.3))
    e = pick()
    @test e.kept == Float64[] && e.kept_model == Float64[] && sort(e.excluded) == [0.1, 0.3, 0.9]
    DBInterface.close!(b.db)
end

@testset "labeled_exposures names the missing migration" begin
    b = _pl_db(mktempdir())
    DBInterface.execute(b.db, "DROP TABLE peak_reviews")
    @test_throws ErrorException HimalayaUI.labeled_exposures(b.db; io = devnull)
    DBInterface.close!(b.db)
end

@testset "labeled_exposures skips exposures pinned only indirectly" begin
    b = _pl_db(mktempdir())
    DBInterface.execute(b.db, "UPDATE peak_reviews SET direct = 0")
    @test [e.id for e in HimalayaUI.labeled_exposures(b.db; io = devnull)] == [b.ff_id]

    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "sample_call_set_state",
        entity_type = "sample", entity_id = b.s_id, payload = Dict(:state => "indexed"))
    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "set_exposure_status",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:status => "accepted"))
    @test [e.id for e in HimalayaUI.labeled_exposures(b.db; io = devnull)] == [b.ff_id]

    HimalayaUI.apply_event!(b.db, HimalayaUI._system_request(); kind = "peak_added",
        entity_type = "exposure", entity_id = b.e_id, payload = Dict(:q => 0.4))
    @test sort([e.id for e in HimalayaUI.labeled_exposures(b.db; io = devnull)]) == sort([b.e_id, b.ff_id])
    DBInterface.close!(b.db)
end

@testset "labeled_exposures include_unlabelled adds exposures nobody curated, without labels" begin
    b = _pl_db(mktempdir())
    other = HimalayaUI.create_exposure!(b.db; experiment_id = b.exp_id, sample_id = b.s_id, filename = "other_tot")
    cp(joinpath(b.adir, "ex_tot.dat"), joinpath(b.adir, "other_tot.dat"))
    DBInterface.execute(b.db, "INSERT INTO auto_peaks (exposure_id, q) VALUES (?, 0.1)", [other])
    DBInterface.execute(b.db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.3)", [other])
    @test !any(e -> e.id == other, HimalayaUI.labeled_exposures(b.db; io = devnull))
    got = HimalayaUI.labeled_exposures(b.db; include_unlabelled = true, io = devnull)
    @test sort([e.id for e in got]) == sort([b.e_id, b.ff_id, other])
    u = only(filter(e -> e.id == other, got))
    @test !u.labelled && !u.form_factor
    @test isempty(u.adds) && isempty(u.kept) && isempty(u.excluded) && isempty(u.kept_model)
    @test only(filter(e -> e.id == b.e_id, got)).adds == [0.4]
    DBInterface.close!(b.db)
end
