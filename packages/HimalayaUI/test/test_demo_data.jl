# Regression floor for the synthetic demo experiment (scripts/make_demo_data.jl,
# docs/demo.md). The demo's "known answers" sit close to findpeaks/indexpeaks
# floors, so a change to peak finding, scoring or `phaseratios` can silently change
# what the demo walkthroughs show. This runs the generator (traces only, ~1 MB) and
# checks the funnel's layout detection, the sample grouping, and that the pipeline's
# auto_group call matches TRUTH's `auto` column often enough across seeds.
using Test, HimalayaUI, Himalaya, Statistics, Random

module DemoGen
    include(joinpath(@__DIR__, "..", "scripts", "make_demo_data.jl"))
end

function _demo_call(dat)
    q, I, σ = HimalayaUI.load_dat(dat)
    pr = Himalaya.findpeaks(q, I, σ)
    g = HimalayaUI.auto_group(Himalaya.indexpeaks(pr.q, pr.sharpness))
    isempty(g) ? "none" : join(sort([string(nameof(Himalaya.phase(x))) for x in g]), "+")
end

@testset "demo data: layout, grouping and phase-call floors" begin
    nexp = sum(s[2] for s in DemoGen.SAMPLES)
    hits = 0; total = 0; spurious = 0
    for seed in (7, 1, 3)
        root = mktempdir()
        redirect_stdout(devnull) do
            DemoGen.main([root, "--no-images", "--seed", string(seed)])
        end

        # Layout detection: the same resolver the ingest funnel calls.
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.integration_pattern == "{name}_tot.dat"
        @test r.image_pattern == "{name}_0_001.tif"
        @test r.setup_file !== nothing && !r.setup_ambiguous

        # Grouping: 2 holders (loads), 12 slots (samples), every scan an exposure.
        metas = HimalayaUI.scan_directory(r.data_dir, r.analysis_dir;
            tif_pattern = r.image_pattern, prp_pattern = r.metadata_pattern,
            dat_pattern = r.integration_pattern)
        @test length(metas) == nexp
        g = HimalayaUI.group_into_samples(metas)
        @test length(g.loads) == 2
        @test all(l -> l.flag == :ok, g.loads)
        @test sum(length(l.samples) for l in g.loads) == length(DemoGen.SAMPLES)

        # Phase calls vs TRUTH's `auto` column. A dropped scan that is not a flare
        # holds no sample, so it must call "none". A flare's call is not pinned
        # (the streak can hide or fake peaks), but it must not invent a phase.
        rows = Dict(begin
                f = split(l, ","; limit = 10)
                f[1] => (bases = split(f[2]), drop = Set(split(f[8])), flare = Set(split(f[9])))
            end for l in readlines(joinpath(root, "TRUTH.csv"))[2:end])
        @test sum(length(r.drop) for r in values(rows)) == sum(length(v) for v in values(DemoGen.BAD_SCANS))
        @test sum(length(r.flare) for r in values(rows)) == count(==(:flare), Iterators.flatten(values.(values(DemoGen.BAD_SCANS))))
        tot = r.analysis_dir
        for (label, _, comps, auto, _) in DemoGen.SAMPLES
            truth = Set(string(nameof(c[1])) for c in comps)
            for base in rows[label].bases
                flare = base in rows[label].flare
                bad = base in rows[label].drop && !flare
                call = _demo_call(joinpath(tot, base * "_tot.dat"))
                if !flare
                    total += 1
                    hits += call == (bad ? "none" : auto)
                end
                call == "none" || (spurious += count(p -> bad || p ∉ truth, split(call, "+")))
            end
        end
    end
    # Floors, not exact counts (CLAUDE.md): 3 seeds × 24 unflared exposures.
    # Measured 2026-11-21 over seeds 1-12 with the beamline-like noise and
    # backgrounds: all 288 unflared calls matched and the 84 empty or bad scans all
    # called "none". Of the 48 flares, 24 kept the sample's call and 24 lost a
    # phase or every phase; none invented one. Raise deliberately.
    @test hits / total >= 0.95
    @test spurious <= 1
end

# Realism floors for what the demo serves. Before 2026-11-21 the demo was far
# from the beamline: frames had a median of 3 counts per pixel, traces had 14%
# point-to-point scatter at high q, and the low-q rise reached 26000× the high-q
# floor. Those made images ~2/3 the size of real ones and peak finding see noise
# real data never has. Targets are aggregate figures from SSRL 1M frames: about
# 100 counts per pixel far from the beam, ~0.35% scatter at high q, a low-q rise
# of ~100× the floor, a beamstop shadow ~10 px in radius and no saturated pixels.
@testset "demo data: frames and traces look like the beamline's" begin
    root = mktempdir()
    redirect_stdout(devnull) do
        DemoGen.main([root, "--no-images", "--seed", "7"])
    end
    rows = Dict(begin f = split(l, ","; limit = 10); f[1] => (bases = split(f[2]), drop = Set(split(f[8]))) end
                for l in readlines(joinpath(root, "TRUTH.csv"))[2:end])
    tot = joinpath(root, "analysis", "automatic_analysis", "tot_files")
    scatter = Float64[]; rise = Float64[]
    for (label, r) in rows, base in r.bases
        base in r.drop && continue
        q, I, _ = HimalayaUI.load_dat(joinpath(tot, base * "_tot.dat"))
        q[end] > 0.38 || continue
        at(x) = I[argmin(abs.(q .- x))]
        push!(scatter, median(abs(I[k] - (I[k-1] + I[k+1]) / 2) / I[k] for k in 2:length(I)-1 if q[k] > 0.3))
        push!(rise, at(0.006) / at(0.38))
    end
    @test median(scatter) <= 0.01
    @test median(rise) <= 500

    qm, dx, dy = DemoGen.q_map()
    model, _ = DemoGen.model_trace(Xoshiro(1), "MO_H2O", DemoGen.SAMPLES[1][3])
    counts = DemoGen.make_frame(Xoshiro(2), model, qm, dx, dy, DemoGen.gap_mask(), 1.0)
    live = counts[counts .>= 0]
    @test median(live) >= 50
    @test count(==(0), live) / length(live) <= 0.03
    @test count(>=(DemoGen.PILATUS_MAX_COUNTS), counts) == 0
    ring = [counts[i] for i in eachindex(counts) if 12 <= hypot(dx[i], dy[i]) <= 18 && counts[i] >= 0]
    @test median(ring) > 0
end

# `--large` writes several realistic-scale experiments (docs/demo.md "Large mode").
# Floors, not exact counts: the mix is drawn from aggregate targets.
@testset "demo data: --large mode" begin
    root = mktempdir()
    redirect_stdout(devnull) do
        DemoGen.main([root, "--large", "--no-images", "--experiments", "2", "--samples", "60", "--seed", "3"])
    end
    exps = sort(filter(d -> isdir(joinpath(root, d)), readdir(root)))
    @test length(exps) == 2
    nframes = Int[]; nophase = 0; nsamples = 0; phases = String[]
    for e in exps
        dir = joinpath(root, e)
        r = HimalayaUI.resolve_experiment_layout(dir)
        @test r.integration_pattern == "{name}_tot.dat"
        @test r.setup_file !== nothing && !r.setup_ambiguous
        metas = HimalayaUI.scan_directory(r.data_dir, r.analysis_dir;
            tif_pattern = r.image_pattern, prp_pattern = r.metadata_pattern,
            dat_pattern = r.integration_pattern)
        g = HimalayaUI.group_into_samples(metas)
        @test all(l -> l.flag == :ok, g.loads)
        rows = [split(l, ","; limit = 10) for l in readlines(joinpath(dir, "TRUTH.csv"))[2:end]]
        found = sum(length(l.samples) for l in g.loads)
        @test found >= 0.9 * length(rows)
        @test length(g.loads) >= 2
        for f in rows
            push!(nframes, length(split(f[2])))
            nsamples += 1
            f[3] == "none" ? (nophase += 1) : append!(phases, split(f[3], "+"))
        end
    end
    @test median(nframes) == 4
    @test 0.25 <= nophase / nsamples <= 0.55
    @test "Square" in phases
    # A finished run leaves a marker; re-running with the same seed is allowed and
    # rewrites the same tree, even after a run was cut short mid-experiment.
    @test isfile(joinpath(root, DemoGen.LARGE_MARKER))
    before = Dict(e => read(joinpath(root, e, "TRUTH.csv")) for e in exps)
    rm(joinpath(root, exps[2], "TRUTH.csv")); rm(joinpath(root, DemoGen.LARGE_MARKER))
    redirect_stdout(devnull) do
        DemoGen.main([root, "--large", "--no-images", "--experiments", "2", "--samples", "60", "--seed", "3"])
    end
    @test all(e -> read(joinpath(root, e, "TRUTH.csv")) == before[e], exps)
    @test isfile(joinpath(root, DemoGen.LARGE_MARKER))

    # The two modes never write into each other's tree, or into anything else.
    curated = mktempdir(); write(joinpath(curated, "TRUTH.csv"), "")
    @test_throws ErrorException DemoGen.check_target(curated, true)
    @test_throws ErrorException DemoGen.check_target(root, false)
    @test_throws ErrorException DemoGen.check_target(let d = mktempdir(); mkdir(joinpath(d, "x")); d end, true)

    # Square is new to the generator: a clear Square trace must index as Square.
    calls = map(1:5) do s
        rng = Xoshiro(s)
        model, _ = DemoGen.model_trace(rng, "sq", DemoGen.Component[(Himalaya.Square, 80.0, 3.0, 2.5e-3)])
        I, σ = DemoGen.measure(rng, model)
        pr = Himalaya.findpeaks(DemoGen.Q, I, σ)
        [string(nameof(Himalaya.phase(x))) for x in HimalayaUI.auto_group(Himalaya.indexpeaks(pr.q, pr.sharpness))]
    end
    @test count(c -> "Square" in c, calls) >= 4
end
