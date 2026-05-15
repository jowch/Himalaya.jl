using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# Compare UX Phase A — route-level view-choice round-trip (Tasks A-5/A-5b).
#
# The comparison view choices (`view_grouping_mode`, `view_show_peak_ticks`,
# `view_show_peak_labels`, spec §6.4) must survive the create + submit
# handlers and round-trip back through GET. These tests drive the HTTP
# surface directly (the canonical `with_test_server` pattern) — no helper
# layer between the test and the wire.
# ─────────────────────────────────────────────────────────────────────────────

# Minimal member payload tied to an exposure. The explicit snapshot means
# the route's `_comparison_member_payload` uses it verbatim and never needs
# analysis output — A-5 only exercises the comparison row, not peak data.
_minimal_member_payload(exposure_id) = Dict{Symbol, Any}(
    :exposure_id   => exposure_id,
    :display_order => 0,
    :band_height   => 1.0,
    :y_offset      => 0.0,
    :normalization => "max",
    :snapshot      => Dict(:effective_peaks => Any[],
                           :confirmed_index => nothing,
                           :analysis_inputs_hash => nothing),
)

# Build a DB with one experiment / sample / exposure (id 1000) via raw
# INSERTs — A-5 only needs a valid `exposures` FK target for the member.
function _routes_cmp_db(tmp::String)
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    DBInterface.execute(db, """INSERT INTO experiments
        (id, name, path, data_dir, analysis_dir)
        VALUES (10, 'exp', '/x', '/x/d', '/x/a')""")
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name) VALUES (100, 10, 'sA')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (1000, 100, 'JC001', 1)")
    db
end

@testset "Comparisons routes — view-choice round-trip (Compare UX A-5)" begin

    @testset "submit echoes view_grouping_mode (Compare UX A-5)" begin
        mktempdir() do tmp
            db = _routes_cmp_db(tmp)
            with_test_server(db) do port, base
                # Create a bare comparison + minimal member.
                resp = HTTP.post("$base/api/comparisons",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "t",
                                     :members => [_minimal_member_payload(1000)])))
                created = JSON3.read(resp.body, Dict{Symbol, Any})
                cid  = created[:id]
                hash = created[:content_hash]

                # Submit with view fields.
                resp2 = HTTP.post("$base/api/comparisons/$(cid)/submit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :title => "t",
                        :members => [_minimal_member_payload(1000)],
                        :expected_content_hash => hash,
                        :view_grouping_mode    => "byPhase",
                        :view_show_peak_ticks  => true,
                        :view_show_peak_labels => false,
                    )))
                @test resp2.status == 200
                submitted = JSON3.read(resp2.body, Dict{Symbol, Any})
                @test submitted[:view_grouping_mode]    == "byPhase"
                @test submitted[:view_show_peak_ticks]  == true
                @test submitted[:view_show_peak_labels] == false

                # GET round-trip.
                resp3 = HTTP.get("$base/api/comparisons/$(cid)", ["X-Username" => "alice"])
                fetched = JSON3.read(resp3.body, Dict{Symbol, Any})
                @test fetched[:view_grouping_mode]    == "byPhase"
                @test fetched[:view_show_peak_ticks]  == true
                @test fetched[:view_show_peak_labels] == false
            end
            close(db)
        end
    end

    @testset "create path echoes view fields (Compare UX A-5b)" begin
        mktempdir() do tmp
            db = _routes_cmp_db(tmp)
            with_test_server(db) do port, base
                resp = HTTP.post("$base/api/comparisons",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :title => "fresh",
                        :members => [_minimal_member_payload(1000)],
                        :view_grouping_mode    => "byPhase",
                        :view_show_peak_ticks  => false,
                        :view_show_peak_labels => true,
                    )))
                @test resp.status == 201
                created = JSON3.read(resp.body, Dict{Symbol, Any})
                @test created[:view_grouping_mode]    == "byPhase"
                @test created[:view_show_peak_ticks]  == false
                @test created[:view_show_peak_labels] == true

                resp2 = HTTP.get("$base/api/comparisons/$(created[:id])",
                                 ["X-Username" => "alice"])
                fetched = JSON3.read(resp2.body, Dict{Symbol, Any})
                @test fetched[:view_grouping_mode]    == "byPhase"
                @test fetched[:view_show_peak_ticks]  == false
                @test fetched[:view_show_peak_labels] == true
            end
            close(db)
        end
    end

    @testset "omitted view fields stay null (Compare UX A-5)" begin
        mktempdir() do tmp
            db = _routes_cmp_db(tmp)
            with_test_server(db) do port, base
                resp = HTTP.post("$base/api/comparisons",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "t",
                                     :members => [_minimal_member_payload(1000)])))
                created = JSON3.read(resp.body, Dict{Symbol, Any})
                @test created[:view_grouping_mode]    === nothing
                @test created[:view_show_peak_ticks]  === nothing
                @test created[:view_show_peak_labels] === nothing
            end
            close(db)
        end
    end

end
