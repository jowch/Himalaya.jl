using Test, SQLite, DBInterface
using HimalayaUI: sqlrows

@testset "sqlrows gives every row the declared column types, whatever is NULL" begin
    db = SQLite.DB()
    DBInterface.execute(db, "CREATE TABLE t (i INTEGER, s TEXT, r REAL)")
    for row in ((1, "a", 1.5), (missing, "b", 2.5), (3, missing, missing), (missing, missing, missing))
        DBInterface.execute(db, "INSERT INTO t VALUES (?, ?, ?)", collect(Any, row))
    end
    rows = sqlrows(DBInterface.execute(db, "SELECT i, s, r FROM t ORDER BY rowid"))
    @test eltype(rows) == NamedTuple{(:i, :s, :r), Tuple{Union{Missing, Int64}, Union{Missing, String}, Union{Missing, Float64}}}
    @test isconcretetype(eltype(rows))
    @test isequal([Tuple(r) for r in rows],
                  [(1, "a", 1.5), (missing, "b", 2.5), (3, missing, missing), (missing, missing, missing)])
end
