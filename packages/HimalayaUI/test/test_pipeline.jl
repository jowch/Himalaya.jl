using Test
using Himalaya: indexpeaks, Index, peaks, score
using HimalayaUI: auto_group

@testset "auto_group" begin
    qs     = [0.1000, 0.1414, 0.2000]
    sharps = [1.0, 0.8, 0.6]

    candidates = indexpeaks(qs, sharps)

    group = auto_group(candidates)

    if !isempty(candidates)
        @test !isempty(group)
        peak_sets = [Set(peaks(idx)) for idx in group]
        for i in eachindex(peak_sets), j in eachindex(peak_sets)
            i == j && continue
            @test isempty(intersect(peak_sets[i], peak_sets[j]))
        end
    end
end

using HimalayaUI: create_schema!, create_experiment!, create_sample!,
                  create_exposure!, persist_analysis!, get_peaks_for_exposure,
                  get_indices_for_exposure, get_groups_for_exposure,
                  load_dat, auto_group
using Himalaya: findpeaks, indexpeaks
using SQLite

@testset "persist_analysis!" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id  = create_experiment!(db; path="/tmp", data_dir="/tmp/data",
                                     analysis_dir="/tmp/analysis")
    s_id    = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id    = create_exposure!(db; sample_id=s_id, filename="example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    peaks_result  = findpeaks(q, I, σ)
    candidates    = indexpeaks(peaks_result.q, peaks_result.sharpness)
    group_indices = auto_group(candidates)

    persist_analysis!(db, e_id, q, I, peaks_result, candidates, group_indices)

    stored_peaks   = get_peaks_for_exposure(db, e_id)
    stored_indices = get_indices_for_exposure(db, e_id)
    stored_groups  = get_groups_for_exposure(db, e_id)

    @test length(stored_peaks) == length(peaks_result.q)
    @test length(stored_indices) == length(candidates)
    @test length(stored_groups) == 1
    @test stored_groups[1].kind   == "auto"
    @test stored_groups[1].active == 1
end

using HimalayaUI: init_experiment!, analyze_exposure!, open_db, get_experiment

@testset "init_experiment!" begin
    tmp = mktempdir()
    data_dir     = joinpath(tmp, "data")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(data_dir)
    mkpath(analysis_dir)

    db = open_db(tmp)
    exp_id = init_experiment!(db;
        name         = "TestExp",
        path         = tmp,
        data_dir     = data_dir,
        analysis_dir = analysis_dir)

    @test exp_id == 1
    exp = get_experiment(db, exp_id)
    @test exp.name == "TestExp"
end

using HimalayaUI: ensure_custom_group!
using DBInterface
using Tables

@testset "persist_analysis! preserves custom-group members across reanalysis" begin
    # Regression: when a peak is added/removed/excluded the frontend triggers
    # an auto-reanalysis. That re-runs persist_analysis! which deletes all
    # `indices` rows and re-inserts with fresh PKs. If we don't translate
    # custom-group memberships from old PK → new PK by semantic identity
    # (phase + basis), the active set ends up empty.
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(tmp)
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    # Promote the auto group to a custom group (mirrors the UI flow when a
    # user clicks an Active-set member).
    custom_id, _ = ensure_custom_group!(db, e_id)

    # Snapshot the custom group's members by semantic identity before reanalysis.
    snap_before = Tables.rowtable(DBInterface.execute(db, """
        SELECT i.phase, i.basis FROM index_group_members m
        JOIN indices i ON i.id = m.index_id
        WHERE m.group_id = ?
        """, [custom_id]))
    @test !isempty(snap_before)
    snap_phases_before = sort(String[String(r.phase) for r in snap_before])

    # Mutate peaks (simulate the user adding a manual peak), then re-analyze.
    DBInterface.execute(db,
        "INSERT INTO peaks (exposure_id, q, intensity, prominence, sharpness, source, excluded)
         VALUES (?, ?, ?, ?, ?, 'manual', 0)",
        [e_id, 0.123, 1.0, 0.5, 0.5])
    analyze_exposure!(db, e_id, analysis_dir)

    # The custom group's members should still be populated, with at least one
    # phase from the original snapshot still present (the data is stable
    # enough that the same phase indexings should re-emerge).
    snap_after = Tables.rowtable(DBInterface.execute(db, """
        SELECT i.phase FROM index_group_members m
        JOIN indices i ON i.id = m.index_id
        WHERE m.group_id = ?
        """, [custom_id]))
    @test !isempty(snap_after)
    snap_phases_after = sort(String[String(r.phase) for r in snap_after])
    # At least one of the originally-active phases survives.
    @test !isempty(intersect(snap_phases_before, snap_phases_after))

    # The custom group should still be `active`; the auto group should be demoted.
    g_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, kind, active FROM index_groups WHERE exposure_id = ?", [e_id]))
    custom_row = first(filter(r -> String(r.kind) == "custom", g_rows))
    auto_row   = first(filter(r -> String(r.kind) == "auto",   g_rows))
    @test Int(custom_row.active) == 1
    @test Int(auto_row.active)   == 0
end

@testset "analyze_exposure! integration" begin
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)

    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(src, joinpath(analysis_dir, "example_tot.dat"))

    db     = open_db(tmp)
    exp_id = init_experiment!(db; path=tmp,
                                   data_dir=joinpath(tmp, "data"),
                                   analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, label="D1", name="UX1")
    e_id   = create_exposure!(db; sample_id=s_id, filename="example_tot")

    analyze_exposure!(db, e_id, analysis_dir)

    @test length(get_peaks_for_exposure(db, e_id))   > 0
    @test length(get_indices_for_exposure(db, e_id)) > 0
    @test length(get_groups_for_exposure(db, e_id))  == 1
end
