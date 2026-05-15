using Test
using SQLite, DBInterface, Tables
using Himalaya
using HimalayaUI: create_schema!, create_experiment!, create_sample!,
                  create_exposure!, effective_peaks, diff_update_auto_peaks!,
                  load_dat

@testset "effective_peaks: sharps_full passthrough matches recompute (no adds)" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id = create_experiment!(db; path = "/tmp", data_dir = "/tmp/d",
                                     analysis_dir = "/tmp/a")
    s_id   = create_sample!(db; experiment_id = exp_id, name = "S",
                                 display_name = "S")
    e_id   = create_exposure!(db; sample_id = s_id, filename = "example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    pr       = Himalaya.findpeaks(q, I, σ)
    diff_update_auto_peaks!(db, e_id, pr, I)

    # No add curations: sharps_full path must not run SG, but result is identical.
    eff_recompute  = effective_peaks(db, e_id, q, I)
    eff_passthrough = effective_peaks(db, e_id, q, I; sharps_full = pr.sharpness_full)
    @test eff_recompute.q         == eff_passthrough.q
    @test eff_recompute.sharpness == eff_passthrough.sharpness
    @test eff_recompute.peak_id   == eff_passthrough.peak_id
    @test eff_recompute.peak_kind == eff_passthrough.peak_kind
end

@testset "effective_peaks: sharps_full passthrough matches recompute (with adds)" begin
    db = SQLite.DB()
    create_schema!(db)
    exp_id = create_experiment!(db; path = "/tmp", data_dir = "/tmp/d",
                                     analysis_dir = "/tmp/a")
    s_id   = create_sample!(db; experiment_id = exp_id, name = "S",
                                 display_name = "S")
    e_id   = create_exposure!(db; sample_id = s_id, filename = "example_tot.dat")

    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    pr       = Himalaya.findpeaks(q, I, σ)
    diff_update_auto_peaks!(db, e_id, pr, I)

    # Add a curation BETWEEN two grid samples so searchsortedfirst's
    # nearest-neighbor selection is exercised (an exact grid-point value
    # would never trigger the i_hi-1 vs i_hi pick).
    k = length(q) ÷ 2
    mid_q = (q[k] + q[k + 1]) / 2  # off-grid; falls strictly between samples
    DBInterface.execute(db,
        "INSERT INTO peak_curations (exposure_id, q, kind) VALUES (?, ?, 'add')",
        [e_id, mid_q])

    eff_recompute   = effective_peaks(db, e_id, q, I)
    eff_passthrough = effective_peaks(db, e_id, q, I; sharps_full = pr.sharpness_full)
    @test eff_recompute.q         == eff_passthrough.q
    @test eff_recompute.sharpness == eff_passthrough.sharpness
    @test eff_recompute.peak_id   == eff_passthrough.peak_id
    @test eff_recompute.peak_kind == eff_passthrough.peak_kind

    # Pin nearest-neighbor selection against argmin oracle.
    sharps_full = pr.sharpness_full
    oracle_sharp = sharps_full[argmin(abs.(q .- mid_q))]
    add_pos = findfirst(==(mid_q), eff_passthrough.q)
    @test add_pos !== nothing
    @test eff_passthrough.sharpness[add_pos] == oracle_sharp
end

@testset "findpeaks return includes sharpness_full" begin
    dat_path = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    q, I, σ  = load_dat(dat_path)
    pr       = Himalaya.findpeaks(q, I, σ)
    @test hasproperty(pr, :sharpness_full)
    @test length(pr.sharpness_full) == length(I)
end
