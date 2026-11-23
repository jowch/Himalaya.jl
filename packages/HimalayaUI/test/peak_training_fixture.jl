using Test, Random, SQLite, DBInterface
using HimalayaUI

function _pt_trace(rng, q0; n = 700, amp = 400.0)
    q = collect(range(0.02, 0.5; length = n))
    I = @. 3000 * (q / 0.02)^-1.8 + 40
    for (k, a) in zip(1:5, (3.0, 2.5, 3.0, 2.0, 1.5))
        @. I += a * amp * (0.0015^2 / ((q - k * q0)^2 + 0.0015^2))
    end
    σ = sqrt.(I)
    q, I .+ σ .* randn(rng, n), σ
end

function _pt_write_dat(path, q, I, σ)
    open(path, "w") do io
        for i in eachindex(q); println(io, q[i], " ", I[i], " ", σ[i]); end
    end
end

# `experiments` experiments of `samples` samples of `per` exposures each; every exposure has
# 4 adds at its true peaks (plus an exclude on every third one) and a confirmed Lamellar call.
function _pt_db(dir; experiments = 3, samples = 4, per = 3, analyze = true)
    rng = MersenneTwister(11)
    db = HimalayaUI.open_db(joinpath(dir, "h.db"))
    eids = Int[]; sids = Int[]; xids = Int[]
    for x in 1:experiments
        adir = joinpath(dir, "analysis$x"); mkpath(adir)
        xid = HimalayaUI.init_experiment!(db; path = joinpath(dir, "exp$x"), data_dir = joinpath(dir, "data$x"), analysis_dir = adir)
        push!(xids, xid)
        for s in 1:samples
            sid = HimalayaUI.create_sample!(db; experiment_id = xid, name = "S$s")
            push!(sids, sid)
            for k in 1:per
                name = "x$(x)s$(s)k$(k)_tot"
                q0 = 0.05 + 0.02 * rand(rng)
                _pt_write_dat(joinpath(adir, name * ".dat"), _pt_trace(rng, q0)...)
                eid = HimalayaUI.create_exposure!(db; experiment_id = xid, sample_id = sid, filename = name)
                push!(eids, eid)
                analyze && HimalayaUI.analyze_exposure!(db, eid, adir)
                for m in 1:4
                    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)", [eid, m * q0])
                end
                length(eids) % 3 == 0 && DBInterface.execute(db,
                    "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', 0.45)", [eid])
                ix = Int(DBInterface.lastrowid(DBInterface.execute(db,
                    "INSERT INTO indices (exposure_id, phase, basis, kind) VALUES (?, 'Lamellar', ?, 'auto')", [eid, q0])))
                DBInterface.execute(db, "INSERT OR REPLACE INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [eid])
                DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [eid, ix])
                DBInterface.execute(db,
                    "INSERT INTO user_actions (action, entity_type, entity_id) VALUES ('assignment_add', 'exposure', ?)", [eid])
                HimalayaUI.record_review!(db, eid; direct = true)
            end
        end
    end
    (; db, path = db.file, eids, sids, xids)
end
