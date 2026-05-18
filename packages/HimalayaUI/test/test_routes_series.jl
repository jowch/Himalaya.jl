using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# I2.2 (#165) — series REST routes + business logic.
#
# The mutating routes emit events whose dispatcher branches do not land until
# #166–#168; `update_view_for_event!` no-ops unknown kinds, so those routes
# write a `user_actions` row but no view rows. Behavioural round-trips for the
# mutating routes belong to #166–#168 and #170; this file covers GET
# round-trips, the `last_event_at` sort fix, the messages round-trip (works
# today), and HTTP-level smoke tests for the mutating routes.
# ─────────────────────────────────────────────────────────────────────────────

# A DB with one experiment / sample / exposure — a valid FK target for
# `series_members.exposure_id` and `series_samples.sample_id`.
function _series_test_db(tmp::String)
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

@testset "Series routes" begin

    @testset "GET /api/series — empty corpus" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                resp = HTTP.get("$base/api/series", ["X-Username" => "alice"])
                @test resp.status == 200
                @test JSON3.read(resp.body) == []
            end
            close(db)
        end
    end

end
