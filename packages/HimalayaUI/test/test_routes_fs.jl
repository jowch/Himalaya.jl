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

@testset "GET /api/fs/validate gates the picker" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            good = mktempdir()
            touch(joinpath(good, "x.tif"))

            # existing dir → ok
            r1 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))")
            @test r1.status == 200
            body1 = JSON3.read(r1.body)
            @test body1.ok == true
            @test body1.message === nothing

            # nonexistent path → not ok + message
            r2 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good * "_nope"))")
            @test r2.status == 200
            body2 = JSON3.read(r2.body)
            @test body2.ok == false
            @test body2.message !== nothing

            # already-an-experiment → not ok
            HimalayaUI.create_experiment!(db; name="e", path=good,
                data_dir=good, analysis_dir=good)
            r3 = call("GET", "/api/fs/validate?path=$(HTTP.escapeuri(good))")
            @test r3.status == 200
            @test JSON3.read(r3.body).ok == false
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
