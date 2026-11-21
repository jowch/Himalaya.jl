using HimalayaUI

# Package loading is the dominant startup cost (~15-20s without a sysimage).
# Exercise the SQLite schema path, the empty-args usage branch, and a real
# ArgParse dispatch against an empty DB (analyze --all errors with no experiments).
mktempdir() do tmp
    db_path = joinpath(tmp, "precompile.db")
    HimalayaUI.open_db(db_path)
    withenv("HIMALAYA_DB_PATH" => db_path) do
        try
            HimalayaUI.main(["analyze", "--all"])
        catch e
            e isa ErrorException && e.msg == "No experiments registered." || rethrow()
        end
    end
end
HimalayaUI.main(String[])
mktempdir(HimalayaUI.run_route_workload)
