using Test, HTTP, JSON3, SQLite
using HimalayaUI

# Filesystem probe routes (spec §6.1). Read-only — no DB writes, no event log,
# no SSE. The DB argument to with_inproc_routes is required by the harness but
# unused by /api/fs/* (these routes never touch the DB).

@testset "GET /api/fs/suggest lists child dirs of a prefix" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            root = mktempdir()
            mkpath(joinpath(root, "alpha"))
            mkpath(joinpath(root, "beta"))
            touch(joinpath(root, "afile.txt"))
            resp = call("GET", "/api/fs/suggest?prefix=$(HTTP.escapeuri(joinpath(root, "a")))")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test any(endswith.(s, "alpha"))
            @test !any(endswith.(s, "afile.txt"))   # files excluded
            @test !any(endswith.(s, "beta"))         # prefix-filtered
        end
    end
end

@testset "GET /api/fs/suggest empty prefix returns empty list" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            resp = call("GET", "/api/fs/suggest")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test length(s) == 0
        end
    end
end

@testset "GET /api/fs/suggest exact dir lists its subdirs" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            root = mktempdir()
            mkpath(joinpath(root, "one"))
            mkpath(joinpath(root, "two"))
            # prefix is the directory itself with trailing slash → isdir branch
            resp = call("GET", "/api/fs/suggest?prefix=$(HTTP.escapeuri(root * "/"))")
            @test resp.status == 200
            s = JSON3.read(resp.body).suggestions
            @test any(endswith.(s, "one"))
            @test any(endswith.(s, "two"))
        end
    end
end

@testset "GET /api/fs/manifest counts by type without writing" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            dir = mktempdir()
            touch(joinpath(dir, "s1.tif")); touch(joinpath(dir, "s1.prp")); touch(joinpath(dir, "s1.dat"))
            touch(joinpath(dir, "s2.tif"))   # missing prp + dat
            q = "path=$(HTTP.escapeuri(dir))&image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            resp = call("GET", "/api/fs/manifest?$q")
            @test resp.status == 200
            m = JSON3.read(resp.body)
            @test m.matched.image == 2
            @test m.matched.metadata == 1
            @test any(u -> u.miss == "metadata", m.unmatched)   # s2 missing prp
            @test isempty(Tables.rowtable(DBInterface.execute(db, "SELECT 1 FROM experiments")))  # no row created
        end
    end
end

@testset "GET /api/fs/manifest flags near-miss image extensions" begin
    # A file that almost matches the image pattern but for a wrong extension/case
    # (.tiff vs .tif, .TIF vs .tif) is a leading cause of a half-failed scan, so the
    # diagnostic surface must see it (§5.5). Emit an `image` miss carrying the
    # offending on-disk filename in `near`.
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            dir = mktempdir()
            touch(joinpath(dir, "good.tif")); touch(joinpath(dir, "good.prp")); touch(joinpath(dir, "good.dat"))
            touch(joinpath(dir, "extra.tiff"))   # near-miss: .tiff superset of .tif
            touch(joinpath(dir, "caps.TIF"))      # near-miss: case
            q = "path=$(HTTP.escapeuri(dir))&image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            m = JSON3.read(call("GET", "/api/fs/manifest?$q").body)
            @test m.matched.image == 1                                       # only good.tif matched exactly
            img_misses = filter(u -> u.miss == "image", m.unmatched)
            @test Set([u.near for u in img_misses]) == Set(["extra.tiff", "caps.TIF"])
            @test all(u -> haskey(u, :file), img_misses)                     # would-be stem present
        end
    end
end

@testset "GET /api/fs/manifest flags near-miss metadata + integration sidecars (6g)" begin
    # An image stem matches but its sidecar is on disk with a wrong-case/extra-char
    # extension (good.PRP for {name}.prp, good.DAT for {name}.dat). The miss row must
    # carry the offending filename in `near` so the scan-failed surface can pair it.
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            data = mktempdir(); analysis = mktempdir()
            touch(joinpath(data, "good.tif"))     # image matches exactly
            touch(joinpath(data, "good.PRP"))     # metadata near-miss: case
            touch(joinpath(analysis, "good.DAT")) # integration near-miss: case (analysis subtree)
            q = "path=$(HTTP.escapeuri(data))&analysis_dir=$(HTTP.escapeuri(analysis))" *
                "&image_pattern=$(HTTP.escapeuri("{name}.tif"))" *
                "&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))" *
                "&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            m = JSON3.read(call("GET", "/api/fs/manifest?$q").body)
            meta_miss  = only(filter(u -> u.miss == "metadata", m.unmatched))
            integ_miss = only(filter(u -> u.miss == "integration", m.unmatched))
            @test meta_miss.file == "good" && meta_miss.near == "good.PRP"
            @test integ_miss.file == "good" && integ_miss.near == "good.DAT"
        end
    end
end

@testset "GET /api/fs/manifest matches integration against analysis_dir" begin
    # Integration (.dat) lives in the analysis subtree, not data_dir (mirrors
    # grouping.jl scan_directory). With analysis_dir supplied, integration is
    # counted there; the per-exposure shortfall surfaces in `unmatched`.
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            data = mktempdir(); analysis = mktempdir()
            for s in ("A", "B", "C")
                touch(joinpath(data, "exp_$s.tif")); touch(joinpath(data, "exp_$s.prp"))
            end
            touch(joinpath(analysis, "exp_A.dat")); touch(joinpath(analysis, "exp_B.dat"))  # C has no .dat
            base = "image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"

            # Without analysis_dir: integration matched against data_dir → 0.
            m0 = JSON3.read(call("GET", "/api/fs/manifest?path=$(HTTP.escapeuri(data))&$base").body)
            @test m0.matched.image == 3
            @test m0.matched.integration == 0

            # With analysis_dir: integration matched against the analysis subtree → 2 (A, B).
            m1 = JSON3.read(call("GET",
                "/api/fs/manifest?path=$(HTTP.escapeuri(data))&analysis_dir=$(HTTP.escapeuri(analysis))&$base").body)
            @test m1.matched.image == 3
            @test m1.matched.integration == 2
            @test any(u -> u.file == "exp_C" && u.miss == "integration", m1.unmatched)
        end
    end
end

@testset "GET /api/fs/manifest 400 for missing dir" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            resp = call("GET", "/api/fs/manifest?path=$(HTTP.escapeuri("/no/such/path/xyz"))")
            @test resp.status == 400
            @test haskey(JSON3.read(resp.body), :error)
        end
    end
end

# ── Helpers for geometry / matched_files manifest tests ──────────────────────
# Inlined from test_ingestion_core.jl (pipeline bucket) — routes bucket does
# NOT include that file; duplicating the two writers avoids cross-bucket deps.

function _fs_write_prp(path::AbstractString;
        beam_energy_ev = 9000.027604502573,
        pipe_length_mm = 1700,
        detector = "Pilatus 1M",
        exposure_time = 15,
        horizontal_position_mm = 58.9)
    open(path, "w") do io
        println(io, "Image file name: $(basename(path)[1:end-4]).tif")
        println(io, "Time this file was written: 26 Apr 2026 18:20:27")
        println(io, "Exposure time=$exposure_time")
        println(io, "Beam energy=$beam_energy_ev eV")
        println(io, "Pipe length=$pipe_length_mm mm")
        println(io, "Horizontal position=$horizontal_position_mm mm")
        println(io, "Detector=$detector")
    end
end

function _fs_write_setup_info(path::AbstractString;
        beam_center_x = 421.409,
        beam_center_y = 836.946,
        mean_distance_mm = 1809.5)
    open(path, "w") do io
        println(io, "File used for calibration:")
        println(io, "../data/agbe_1p7m_S1958_0_001.tif")
        println(io, "")
        println(io, "Beam center is at ($beam_center_x, $beam_center_y)")
        println(io, "")
        println(io, "Mean distance:\t\t $mean_distance_mm mm")
    end
end

@testset "GET /api/fs/manifest includes geometry from PRP + setup_info" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            dir = mktempdir()
            # Write image + PRP pairs (the real parse_prp format).
            for stem in ("a_001", "b_002")
                touch(joinpath(dir, "$stem.tif"))
                _fs_write_prp(joinpath(dir, "$stem.prp");
                    beam_energy_ev = 9000.0, pipe_length_mm = 1700,
                    detector = "Pilatus 1M")
            end
            # Write an analysis/setup_info file (the real parse_setup_info format).
            analysis_dir = joinpath(dir, "analysis")
            mkpath(analysis_dir)
            _fs_write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt");
                beam_center_x = 421.409, beam_center_y = 836.946, mean_distance_mm = 1809.5)

            q = "path=$(HTTP.escapeuri(dir))&image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            resp = call("GET", "/api/fs/manifest?$q")
            @test resp.status == 200
            m = JSON3.read(resp.body)

            # Legacy fields still intact.
            @test m.matched.image == 2
            @test m.matched.metadata == 2

            # Geometry block present and populated from setup_info (beam center) and PRP (energy).
            @test haskey(m, :geometry)
            geo = m.geometry
            @test geo.beam_center_x ≈ 421.409
            @test geo.beam_center_x_source == "setup"
            @test geo.beam_center_y ≈ 836.946
            @test geo.beam_center_y_source == "setup"
            @test geo.flight_path_m ≈ 1809.5 / 1000.0
            @test geo.flight_path_m_source == "setup"
            @test geo.energy_kev ≈ 9.0
            @test geo.energy_kev_source == "prp"
            @test geo.pixel_size_um ≈ 172.0   # Pilatus 1M
            @test geo.pixel_size_um_source == "prp"

            # discrepancies key present (may be empty for a consistent run).
            @test haskey(m, :discrepancies)

            # matched_files: sorted image filenames, capped at 12.
            @test haskey(m, :matched_files)
            mf = m.matched_files
            @test length(mf) == 2
            @test "a_001.tif" in mf
            @test "b_002.tif" in mf
            @test issorted(mf)
        end
    end
end

@testset "GET /api/fs/manifest geometry gracefully degrades when no PRP/setup" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            dir = mktempdir()
            # Only image files — no PRP, no setup_info.
            touch(joinpath(dir, "x_001.tif"))
            touch(joinpath(dir, "x_002.tif"))

            q = "path=$(HTTP.escapeuri(dir))&image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            resp = call("GET", "/api/fs/manifest?$q")
            @test resp.status == 200
            m = JSON3.read(resp.body)

            # Returns 200 (not 500) with all-null geometry.
            @test haskey(m, :geometry)
            geo = m.geometry
            @test geo.beam_center_x === nothing
            @test geo.beam_center_y === nothing
            @test geo.flight_path_m === nothing
            @test geo.energy_kev === nothing
            @test geo.pixel_size_um === nothing

            # matched_files present even with no metadata.
            @test haskey(m, :matched_files)
            @test length(m.matched_files) == 2
        end
    end
end

@testset "GET /api/fs/manifest matched_files capped at 12" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            dir = mktempdir()
            for i in 1:20
                touch(joinpath(dir, "s$(lpad(i, 3, '0')).tif"))
            end
            q = "path=$(HTTP.escapeuri(dir))&image_pattern=$(HTTP.escapeuri("{name}.tif"))&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))&integration_pattern=$(HTTP.escapeuri("{name}.dat"))"
            resp = call("GET", "/api/fs/manifest?$q")
            @test resp.status == 200
            m = JSON3.read(resp.body)
            @test m.matched.image == 20
            @test length(m.matched_files) == 12
            @test issorted(m.matched_files)
        end
    end
end

# ── Structural layout resolver (funnel resolution) ───────────────────────────

"""Build a real-shaped experiment ROOT: root/data/{tif,prp}, root/analysis/setup_info,
   root/analysis/automatic_analysis/tot_files/{dat}. Returns the root path."""
function _fs_make_experiment_root(root::AbstractString)
    data = joinpath(root, "data"); mkpath(data)
    for stem in ("a_001", "b_002")
        touch(joinpath(data, "$stem.tif"))
        _fs_write_prp(joinpath(data, "$stem.prp"))
    end
    analysis = joinpath(root, "analysis"); mkpath(analysis)
    _fs_write_setup_info(joinpath(analysis, "setup_info_20260425_181705.txt"))
    totdir = joinpath(analysis, "automatic_analysis", "tot_files"); mkpath(totdir)
    touch(joinpath(totdir, "a_001_tot.dat"))
    return root
end

@testset "resolve_experiment_layout discovers data/analysis/setup structurally" begin
    mktempdir() do tmp
        root = joinpath(tmp, "1p7m"); mkpath(root)
        _fs_make_experiment_root(root)
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.name == "1p7m"
        @test r.data_dir == joinpath(root, "data")
        # analysis_dir points at the dir holding the integration .dat (the SSRL
        # tot_files convention), discovered from the {base}_tot.dat present there.
        @test r.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.integration_pattern == "{name}_tot.dat"
        @test r.image_pattern == "{name}.tif"        # fixture stems carry no _0_NNN frame
        @test r.metadata_pattern == "{name}.prp"
        @test r.setup_file !== nothing && endswith(r.setup_file, "setup_info_20260425_181705.txt")
        @test r.setup_ambiguous == false
    end
end

@testset "resolve_experiment_layout keys tot patterns off the frame suffix" begin
    # SSRL data: stems carry a _0_001 frame; tot files are per-base ({base}_tot.dat).
    # The resolver must suggest `{name}_0_001.tif` so the scan's stem resolves to base.
    mktempdir() do tmp
        root = joinpath(tmp, "ssrl"); mkpath(root)
        data = joinpath(root, "data"); mkpath(data)
        for stem in ("HA_1_001_S1_0_001", "HA_1_002_S2_0_001")
            touch(joinpath(data, "$stem.tif")); _fs_write_prp(joinpath(data, "$stem.prp"))
        end
        totdir = joinpath(root, "analysis", "automatic_analysis", "tot_files"); mkpath(totdir)
        touch(joinpath(totdir, "HA_1_001_S1_tot.dat")); touch(joinpath(totdir, "HA_1_002_S2_tot.dat"))
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.analysis_dir == totdir
        @test r.image_pattern == "{name}_0_001.tif"
        @test r.metadata_pattern == "{name}_0_001.prp"
        @test r.integration_pattern == "{name}_tot.dat"
    end
end

@testset "resolve_experiment_layout flags no-setup as ambiguous" begin
    mktempdir() do tmp
        root = joinpath(tmp, "bare"); mkpath(root)
        data = joinpath(root, "data"); mkpath(data)
        touch(joinpath(data, "a_001.tif"))  # data, but no analysis/ and no setup_info
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.data_dir == joinpath(root, "data")
        @test r.analysis_dir === nothing
        @test r.setup_file === nothing
        @test r.setup_ambiguous == true
    end
end

@testset "GET /api/fs/manifest uses an explicit setup_file for geometry" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            # data_dir has the PRPs; the setup file lives OUTSIDE it (the real
            # split). Passing setup_file explicitly is what makes geometry derive.
            data = mktempdir()
            for stem in ("a_001", "b_002")
                touch(joinpath(data, "$stem.tif")); _fs_write_prp(joinpath(data, "$stem.prp"))
            end
            setup = joinpath(mktempdir(), "setup_info_20260425_181705.txt")
            _fs_write_setup_info(setup)
            q = "path=$(HTTP.escapeuri(data))&setup_file=$(HTTP.escapeuri(setup))" *
                "&metadata_pattern=$(HTTP.escapeuri("{name}.prp"))"
            resp = call("GET", "/api/fs/manifest?$q")
            @test resp.status == 200
            geo = JSON3.read(resp.body).geometry
            @test geo.beam_center_x ≈ 421.409
            @test geo.beam_center_x_source == "setup"
        end
    end
end

function _fs_make_ssrl_setup(root::AbstractString; stems, traced = stems, distance_mm = 1809.5)
    data = joinpath(root, "data"); mkpath(data)
    for stem in stems
        touch(joinpath(data, "$(stem)_0_001.tif"))
        _fs_write_prp(joinpath(data, "$(stem)_0_001.prp"))
    end
    analysis = joinpath(root, "analysis"); mkpath(analysis)
    _fs_write_setup_info(joinpath(analysis, "setup_info_20260425_181705.txt"); mean_distance_mm = distance_mm)
    totdir = joinpath(analysis, "automatic_analysis", "tot_files"); mkpath(totdir)
    for stem in traced
        touch(joinpath(totdir, "$(stem)_tot.dat"))
    end
    return root
end

@testset "resolve_experiment_layout outvotes a first image that has no trace" begin
    mktempdir() do tmp
        root = joinpath(tmp, "setup1")
        stems = ["AA_cal", "HA_1_001_S1", "HA_1_002_S2", "HA_1_003_S3"]
        _fs_make_ssrl_setup(root; stems, traced = stems[2:end])
        touch(joinpath(root, "analysis", "automatic_analysis", "tot_files", "AA_cal_0_001.dat"))
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.image_pattern == "{name}_0_001.tif"
        @test r.integration_pattern == "{name}_tot.dat"
    end
end

@testset "GET /api/fs/preview lists one setup for a setup folder" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            root = joinpath(mktempdir(), "20261103", "setup1")
            _fs_make_ssrl_setup(root; stems = ["HA_1_001_S1", "HA_1_002_S2"], traced = ["HA_1_001_S1"])
            resp = call("GET", "/api/fs/preview?path=$(HTTP.escapeuri(root))")
            @test resp.status == 200
            setups = JSON3.read(resp.body).setups
            @test length(setups) == 1
            s = only(setups)
            @test s.root == root
            @test s.name == "20261103/setup1"
            @test s.data_dir == joinpath(root, "data")
            @test s.patterns.integration == "{name}_tot.dat"
            @test s.exposures == 2
            @test s.without_trace == 1
            @test s.without_metadata == 0
            @test s.geometry.flight_path_m_source == "setup"
            @test s.duplicate_of === nothing
            @test call("GET", "/api/fs/preview?path=$(HTTP.escapeuri("/no/such/xyz"))").status == 400
        end
    end
end

@testset "GET /api/fs/preview lists every setup under a beamtime folder" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            day = joinpath(mktempdir(), "20261103")
            _fs_make_ssrl_setup(joinpath(day, "setup2"); stems = ["HB_1_001_S9"], distance_mm = 500.0)
            _fs_make_ssrl_setup(joinpath(day, "setup1"); stems = ["HA_1_001_S1", "HA_1_002_S2"])
            mkpath(joinpath(day, "notes")); touch(joinpath(day, "notes", "plan.txt"))
            resp = call("GET", "/api/fs/preview?path=$(HTTP.escapeuri(day))")
            @test resp.status == 200
            setups = JSON3.read(resp.body).setups
            @test [s.name for s in setups] == ["20261103/setup1", "20261103/setup2"]
            @test [s.exposures for s in setups] == [2, 1]
            @test setups[1].setup_file == joinpath(day, "setup1", "analysis", "setup_info_20260425_181705.txt")
            @test setups[2].setup_file == joinpath(day, "setup2", "analysis", "setup_info_20260425_181705.txt")
            @test setups[1].geometry.flight_path_m != setups[2].geometry.flight_path_m
        end
    end
end

@testset "GET /api/fs/preview names the experiment already using a setup" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        day = joinpath(mktempdir(), "20261103")
        root = _fs_make_ssrl_setup(joinpath(day, "setup1"); stems = ["HA_1_001_S1"])
        HimalayaUI.create_experiment!(db; name = "first", path = joinpath(root, "data"),
            data_dir = joinpath(root, "data"), analysis_dir = root)
        with_inproc_routes(db) do call
            s = only(JSON3.read(call("GET", "/api/fs/preview?path=$(HTTP.escapeuri(day))").body).setups)
            @test s.duplicate_of.name == "first"
            empty = mktempdir()
            @test isempty(JSON3.read(call("GET", "/api/fs/preview?path=$(HTTP.escapeuri(empty))").body).setups)
        end
    end
end

@testset "the layout vote samples stems across a large folder" begin
    mktempdir() do tmp
        root = joinpath(tmp, "setup1")
        stems = ["S$(lpad(i, 3, '0'))" for i in 1:30]
        _fs_make_ssrl_setup(root; stems, traced = String[])
        flat = joinpath(root, "analysis", "flat"); mkpath(flat)
        for s in stems[1:5]
            touch(joinpath(flat, "$(s)_0_001.dat"))
        end
        for s in stems[6:end]
            touch(joinpath(root, "analysis", "automatic_analysis", "tot_files", "$(s)_tot.dat"))
        end
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.integration_pattern == "{name}_tot.dat"
    end
end

@testset "a setup copied before any trace was written still gets the SSRL layout" begin
    mktempdir() do tmp
        root = joinpath(tmp, "setup1")
        _fs_make_ssrl_setup(root; stems = ["HA_1_001_S1", "HA_1_002_S2"], traced = String[])
        touch(joinpath(root, "data", "._HA_1_001_S1_0_001.tif"))
        r = HimalayaUI.resolve_experiment_layout(root)
        @test r.analysis_dir == joinpath(root, "analysis", "automatic_analysis", "tot_files")
        @test r.image_pattern == "{name}_0_001.tif"
        @test r.metadata_pattern == "{name}_0_001.prp"
        @test r.integration_pattern == "{name}_tot.dat"
    end
end

@testset "a stray image at the top of a beamtime folder does not hide its setups" begin
    mktempdir() do tmp
        day = joinpath(tmp, "20261103")
        _fs_make_ssrl_setup(joinpath(day, "setup1"); stems = ["HA_1_001_S1"])
        _fs_make_ssrl_setup(joinpath(day, "setup2"); stems = ["HA_1_001_S1"])
        touch(joinpath(day, "stray.tif"))
        @test HimalayaUI.setup_roots(day) == [joinpath(day, "setup1"), joinpath(day, "setup2")]
    end
end

@testset "a duplicate is found through a trailing slash or a symlink" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        day = joinpath(mktempdir(), "20261103")
        root = _fs_make_ssrl_setup(joinpath(day, "setup1"); stems = ["HA_1_001_S1"])
        HimalayaUI.create_experiment!(db; name = "first", path = joinpath(root, "data"),
            data_dir = joinpath(root, "data"), analysis_dir = root)
        link = joinpath(mktempdir(), "linked"); symlink(root, link)
        @test HimalayaUI.experiment_using_dir(db, joinpath(root, "data") * "/").name == "first"
        @test HimalayaUI.experiment_using_dir(db, joinpath(link, "data")).name == "first"
        with_inproc_routes(db) do call
            resp = call("POST", "/api/experiments";
                headers = ["Content-Type" => "application/json"],
                body = Vector{UInt8}(JSON3.write(Dict("data_dir" => joinpath(link, "data") * "/", "name" => "again"))))
            @test resp.status == 409
        end
    end
end

@testset "the duplicate check reads stored folders as text and never visits them" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        gone = joinpath(tmp, "unmounted-share", "setup1", "data")
        HimalayaUI.create_experiment!(db; name = "old", path = gone, data_dir = gone * "/", analysis_dir = gone)
        @test HimalayaUI.experiment_using_dir(db, gone).name == "old"
        real = mkpath(joinpath(tmp, "real", "data"))
        link = joinpath(tmp, "link"); symlink(real, link)
        HimalayaUI.create_experiment!(db; name = "via link", path = link, data_dir = link, analysis_dir = link)
        @test HimalayaUI.experiment_using_dir(db, link).name == "via link"
        @test HimalayaUI.experiment_using_dir(db, real) === nothing
    end
end

@testset "a new experiment stores its folder canonically, so a symlinked spelling is caught" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        real = mkpath(joinpath(tmp, "real", "data"))
        link = joinpath(tmp, "link"); symlink(joinpath(tmp, "real"), link)
        with_inproc_routes(db) do call
            post(path) = call("POST", "/api/experiments";
                headers = ["Content-Type" => "application/json"],
                body = Vector{UInt8}(JSON3.write(Dict("data_dir" => path, "name" => basename(path)))))
            r = post(joinpath(link, "data") * "/")
            @test r.status == 202
            stored = first(Tables.rowtable(DBInterface.execute(db, "SELECT data_dir FROM experiments WHERE id = ?",
                [JSON3.read(r.body).id]))).data_dir
            @test stored == realpath(real)
            @test post(real).status == 409
        end
    end
end
