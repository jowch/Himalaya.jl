using Test, Dates, SQLite, DBInterface, Tables, Statistics, Random, HimalayaUI
isdefined(@__MODULE__, :fresh_db) || include("test_ingestion_core.jl")

module DemoGenAS
    include(joinpath(@__DIR__, "..", "scripts", "make_demo_data.jl"))
end

module DemoGenWithout
    include(joinpath(@__DIR__, "..", "scripts", "make_demo_data.jl"))
    write_autosampler!(out, seed, images) = (acquisitions = 0, runs = 0, calibration = 0)
end

const _DAS = DemoGenAS

function _das_generate(; seed = 7)
    root = mktempdir()
    redirect_stdout(devnull) do
        _DAS.main([root, "--no-images", "--seed", string(seed)])
    end
    root
end

function _das_truth(root)
    lines = readlines(joinpath(root, _DAS.AUTOSAMPLER_TRUTH_FILE))
    @test lines[1] == _DAS.AUTOSAMPLER_TRUTH_HEADER
    [NamedTuple{Tuple(Symbol.(split(lines[1], ",")))}(Tuple(split(l, ","; limit = 15))) for l in lines[2:end]]
end

_das_files(dir) = sort([relpath(joinpath(d, f), dir) for (d, _, fs) in walkdir(dir) for f in fs])

function _das_ingest(root)
    db = fresh_db()
    eid = HimalayaUI.create_experiment!(db; name = "demo", path = root)
    for (i, setup_root) in enumerate(HimalayaUI.setup_roots(root))
        r = HimalayaUI.resolve_experiment_layout(setup_root)
        HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = r.data_dir, analysis_dir = r.analysis_dir,
            root_dir = setup_root, position = i - 1, mode = r.mode,
            image_pattern = r.image_pattern, metadata_pattern = r.metadata_pattern,
            integration_pattern = r.integration_pattern)
    end
    HimalayaUI.scan_and_group!(db, eid)
    (db, eid)
end

const _DAS_ROWS = """
    SELECT s.id AS id, s.name AS name, e.filename AS filename, e.frames_total AS total,
           e.frames_rejected AS rejected, l.load_index AS load_index, st.mode AS mode, st.id AS setup_id,
           (SELECT COUNT(*) FROM exposures x WHERE x.sample_id = s.id) AS n_exposures
      FROM samples s JOIN setups st ON st.id = s.setup_id
      JOIN exposures e ON e.sample_id = s.id LEFT JOIN loads l ON l.id = s.load_id
     WHERE st.mode = 'autosampler' ORDER BY s.id"""

@testset "demo data: the autosampler setup" begin
    root = _das_generate()
    truth = _das_truth(root)
    by_name = Dict(t.sample => t for t in truth)

    @testset "the demo root holds two setups, the capillary one unchanged in kind" begin
        @test HimalayaUI.setup_roots(root) ==
              [joinpath(root, _DAS.AUTOSAMPLER_DIR), joinpath(root, _DAS.CAPILLARY_DIR)]
        @test !isdir(joinpath(root, "data")) && !isdir(joinpath(root, "analysis"))
        cap = HimalayaUI.resolve_experiment_layout(joinpath(root, _DAS.CAPILLARY_DIR))
        @test cap.mode == "capillary"
        r = HimalayaUI.resolve_experiment_layout(joinpath(root, _DAS.AUTOSAMPLER_DIR))
        @test r.mode == "autosampler"
        @test r.setup_file !== nothing && !r.setup_ambiguous
        @test r.integration_pattern == "{name}_tot.dat"
        @test r.image_pattern == "{name}_0_001.tif"
    end

    @testset "the files have the shape of a real autosampler folder" begin
        data = joinpath(root, _DAS.AUTOSAMPLER_DIR, "data")
        names = readdir(data)
        @test count(f -> occursin(HimalayaUI.AUTOSAMPLER_RUN_FILE, f), names) == 2
        acquisitions = filter(t -> t.kind != "calibration", truth)
        @test length(acquisitions) == 16
        for t in acquisitions
            @test all(f -> isfile(joinpath(data, "$(t.stem)_0_$(lpad(f, 3, '0')).tif")) &&
                           isfile(joinpath(data, "$(t.stem)_0_$(lpad(f, 3, '0')).prp")), 1:_DAS.AS_FRAMES)
            @test occursin(r"^[A-Za-z]+_\d{2}[A-H]_[SB]\d{4}$", t.stem)
        end
        calibration = only(filter(t -> t.kind == "calibration", truth))
        @test isfile(joinpath(data, calibration.stem * "_0_001.tif"))
        tot = joinpath(root, _DAS.AUTOSAMPLER_DIR, "analysis", "automatic_analysis", "tot_files")
        @test sort(readdir(tot)) == sort([t.stem * "_tot.dat" for t in acquisitions])
        logs = joinpath(root, _DAS.AUTOSAMPLER_DIR, "analysis", "automatic_analysis", "sastool")
        @test sort(readdir(logs)) == sort([t.stem * "_0_001.log" for t in acquisitions])
        @test occursin("-i yes i2 ", read(joinpath(root, _DAS.AUTOSAMPLER_DIR, "analysis", "integ.mpp"), String))
        horizontal = Set(match(r"Horizontal position=(\S+)", read(joinpath(data, f), String)).captures[1]
                         for f in names if endswith(f, ".prp"))
        @test length(horizontal) == 1
    end

    @testset "generation is deterministic, and seeds differ" begin
        again = _das_generate()
        auto = joinpath(root, _DAS.AUTOSAMPLER_DIR)
        @test _das_files(auto) == _das_files(joinpath(again, _DAS.AUTOSAMPLER_DIR))
        for f in _das_files(auto)
            @test read(joinpath(auto, f)) == read(joinpath(again, _DAS.AUTOSAMPLER_DIR, f))
        end
        @test read(joinpath(root, _DAS.AUTOSAMPLER_TRUTH_FILE)) == read(joinpath(again, _DAS.AUTOSAMPLER_TRUTH_FILE))
        other = _das_generate(seed = 8)
        stem = first(t.stem for t in truth if t.kind == "sample")
        tot(r) = read(joinpath(r, _DAS.AUTOSAMPLER_DIR, "analysis", "automatic_analysis", "tot_files", stem * "_tot.dat"))
        @test tot(root) != tot(other)
    end

    @testset "the capillary setup does not depend on the autosampler" begin
        bare = mktempdir()
        redirect_stdout(devnull) do
            DemoGenWithout.main([bare, "--no-images", "--seed", "7"])
        end
        files(dir) = sort([relpath(joinpath(d, f), dir) for (d, _, fs) in walkdir(dir) for f in fs])
        capillary = joinpath(root, _DAS.CAPILLARY_DIR)
        @test files(joinpath(bare, _DAS.CAPILLARY_DIR)) == files(capillary)
        @test length(files(capillary)) > 3 * length(_DAS.SAMPLES)
        for f in files(capillary)
            @test read(joinpath(bare, _DAS.CAPILLARY_DIR, f)) == read(joinpath(capillary, f))
        end
        @test read(joinpath(bare, "TRUTH.csv")) == read(joinpath(root, "TRUTH.csv"))
        @test !ispath(joinpath(bare, _DAS.AUTOSAMPLER_TRUTH_FILE))
    end

    @testset "the special cases TRUTH records are in the files" begin
        flagged(flag) = [t for t in truth if flag in split(t.flags)]
        @test length(flagged("i2_collapsed")) == 1
        @test length(flagged("one_frame_left")) == 1
        @test only(flagged("one_frame_left")).frames_rejected == string(_DAS.AS_FRAMES - 1)
        @test !isempty(flagged("several_rejected"))
        @test all(t -> 2 <= parse(Int, t.frames_rejected) < _DAS.AS_FRAMES - 1, flagged("several_rejected"))
        @test !isempty(flagged("no_closing_buffer"))
        runs = unique(t.load for t in truth if t.kind != "calibration")
        @test length(runs) == 3
        missing_closing = unique(t.load for t in flagged("no_closing_buffer"))
        @test length(missing_closing) == 1
        for load in runs
            members = [t for t in truth if t.load == load]
            closing = last(members).kind == "buffer" && length(members) > 1
            @test closing == (load != only(missing_closing))
        end
        repeated = flagged("repeated_well")
        @test !isempty(repeated)
        for t in repeated
            earlier = by_name[t.repeat_of]
            @test earlier.well == t.well && earlier.run == t.run && earlier.stem != t.stem
            @test parse(Int, earlier.load) < parse(Int, t.load)
        end
        firsts = flagged("first_buffer")
        @test length(firsts) == 3
        for t in truth
            t.kind == "sample" && @test by_name[t.default_background].kind == "buffer"
            t.kind == "sample" && @test by_name[t.default_background].load == t.load
            t.kind != "sample" && @test isempty(t.default_background)
        end
        @test all(t -> t.default_background in (f.sample for f in firsts), filter(t -> t.kind == "sample", truth))
    end

    @testset "the app ingests it as an autosampler setup of single-image samples" begin
        (db, eid) = _das_ingest(root)
        rows = Tables.rowtable(DBInterface.execute(db, _DAS_ROWS))
        acquisitions = filter(t -> t.kind != "calibration", truth)
        @test Set(String(r.name) for r in rows) == Set(vcat([t.sample for t in truth]))
        @test all(r -> r.n_exposures == 1, rows)
        @test length(unique(r.setup_id for r in rows)) == 1
        @test length(rows) == length(truth)
        for r in rows
            t = by_name[String(r.name)]
            @test Int(r.load_index) == parse(Int, t.load)
            isempty(t.frames_total) && continue
            @test (Int(r.total), Int(r.rejected)) == (parse(Int, t.frames_total), parse(Int, t.frames_rejected))
        end
        @test !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM setups WHERE experiment_id = ? AND mode = 'capillary'", [eid])))

        ids = Dict(String(r.name) => Int(r.id) for r in rows)
        defaults = HimalayaUI.default_backgrounds(db, eid)
        for t in truth
            expected = isempty(t.default_background) ? nothing : ids[t.default_background]
            @test get(defaults, ids[t.sample], :absent) == expected
        end

        @test HimalayaUI.read_trace_normalisations!(db; experiment_id = eid) >= 1
        st = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT trace_normalisation AS by, normalisation_base AS base FROM setups WHERE mode = 'autosampler'")))
        @test st.by == HimalayaUI.NORMALISED_BY_I2
        @test st.base == _DAS.GRAND_BASE

        monitors = Tables.rowtable(DBInterface.execute(db, """
            SELECT s.name AS name, e.transmitted_counts AS t, e.incident_monitor AS i, e.counting_time AS c,
                   e.setup_id AS setup_id
              FROM exposures e JOIN samples s ON s.id = e.sample_id
              JOIN setups st ON st.id = e.setup_id WHERE st.mode = 'autosampler'"""))
        median_t = HimalayaUI.setup_transmission_median(db, only(unique(m.setup_id for m in monitors)))
        collapsed = [String(m.name) for m in monitors if !ismissing(m.t) &&
                     HimalayaUI.is_i2_collapsed(m.t, m.i, m.c, median_t)]
        @test collapsed == [only(t.sample for t in truth if occursin("i2_collapsed", t.flags))]
    end

    @testset "traces look like divided-by-I_2 solution scattering" begin
        tot = joinpath(root, _DAS.AUTOSAMPLER_DIR, "analysis", "automatic_analysis", "tot_files")
        trace(t) = HimalayaUI.load_dat(joinpath(tot, t.stem * "_tot.dat"))
        buffers = [t for t in truth if t.kind == "buffer"]
        function highq(t)
            q, I, _ = trace(t)
            median(I[q .> 0.25])
        end
        normal = [highq(t) for t in truth if t.kind != "calibration" && !occursin("i2_collapsed", t.flags)]
        @test all(0.9 * _DAS.SOLVENT_FLOOR .< normal .< 1.1 * _DAS.SOLVENT_FLOOR)
        @test std(normal) / mean(normal) < 0.03
        collapsed = only(t for t in truth if occursin("i2_collapsed", t.flags))
        @test highq(collapsed) > 15 * mean(normal)

        for t in buffers
            q, I, σ = trace(t)
            @test I[findfirst(>(0.006), q)] > 20 * median(I[q .> 0.25])
            hi = findall(>(0.3), q)[2:end-1]
            scatter = median(abs(I[k] - (I[k-1] + I[k+1]) / 2) for k in hi)
            @test 0.6 < scatter / (0.83 * median(σ[hi])) < 1.6
        end

        reference = first(b for b in buffers if b.load == "101")
        _, background, _ = trace(reference)
        for t in truth
            t.kind == "sample" && t.load == "101" && !occursin("i2_collapsed", t.flags) || continue
            q, I, _ = trace(t)
            diff = I .- background
            expected = parse(Float64, t.minimum_q)
            window = findall(x -> 0.5 * expected <= x <= 1.5 * expected, q)
            smooth = [mean(diff[max(1, k - 3):min(end, k + 3)]) for k in window]
            found = q[window[argmin(smooth)]]
            @test abs(found - expected) / expected < 0.15
            @test maximum(diff[q .< 0.5 * expected]) > 10 * abs(diff[window[argmin(smooth)]])
        end
    end
end
