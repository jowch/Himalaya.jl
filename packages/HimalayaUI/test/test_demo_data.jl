# Regression floor for the synthetic demo experiment (scripts/make_demo_data.jl,
# docs/demo.md). The demo's "known answers" sit close to findpeaks/indexpeaks
# floors, so a change to peak finding, scoring or `phaseratios` can silently change
# what the demo walkthroughs show. This runs the generator (traces only, ~1 MB) and
# checks the funnel's layout detection, the sample grouping, and that the pipeline's
# auto_group call matches TRUTH's `auto` column often enough across seeds.
using Test, HimalayaUI, Himalaya

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

        # Phase calls vs TRUTH's `auto` column.
        tot = r.analysis_dir
        for (label, _, comps, auto, _) in DemoGen.SAMPLES
            truth = Set(string(nameof(c[1])) for c in comps)
            for f in filter(f -> startswith(f, label * "_S"), readdir(tot))
                call = _demo_call(joinpath(tot, f))
                total += 1
                hits += call == auto
                call == "none" || (spurious += count(p -> p ∉ truth, split(call, "+")))
            end
        end
    end
    # Floors, not exact counts (CLAUDE.md): 3 seeds × 26 exposures. Measured
    # 2026-11-19 over seeds 1-12: every call matched except MO_AMP_10 on 3 of 36
    # scans (an extra Pn3m lattice claims HII's first peak). Raise deliberately.
    @test hits / total >= 0.95
    @test spurious <= 1
end
