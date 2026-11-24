using Test, SQLite, DBInterface
using HimalayaUI: sqlrows

function _rows_table(values)
    db = SQLite.DB()
    DBInterface.execute(db, "CREATE TABLE t (i INTEGER, s TEXT, r REAL)")
    for row in values
        DBInterface.execute(db, "INSERT INTO t VALUES (?, ?, ?)", collect(Any, row))
    end
    db
end

_read_rows(db) = sqlrows(DBInterface.execute(db, "SELECT i, s, r FROM t ORDER BY rowid"))

@testset "sqlrows gives every row the declared column types, whatever is NULL" begin
    rows = _read_rows(_rows_table(((1, "a", 1.5), (missing, "b", 2.5), (3, missing, missing))))
    @test eltype(rows) == NamedTuple{(:i, :s, :r), Tuple{Union{Missing, Int64}, Union{Missing, String}, Union{Missing, Float64}}}
    @test isequal([Tuple(r) for r in rows], [(1, "a", 1.5), (missing, "b", 2.5), (3, missing, missing)])
end

@testset "sqlrows falls back to Any for an undeclared column or values that don't match the declared type" begin
    db = _rows_table(((1, "a", 1.5),))
    DBInterface.execute(db, "INSERT INTO t VALUES ('not a number', 'b', 2.5)")
    rows = sqlrows(DBInterface.execute(db, "SELECT i, i + 1 AS j FROM t ORDER BY rowid"))
    @test fieldtypes(eltype(rows)) == (Any, Any)
    @test rows[2].i == "not a number"
end

@testset "sqlrows compiles nothing when it meets a new NULL pattern" begin
    _read_rows(_rows_table(((1, "a", 1.5), (missing, "b", 2.5))))
    unseen = _rows_table(((missing, missing, 9.5), (7, missing, missing)))
    Base.cumulative_compile_timing(true)
    before = Base.cumulative_compile_time_ns()[1]
    _read_rows(unseen)
    compiled_ns = Base.cumulative_compile_time_ns()[1] - before
    Base.cumulative_compile_timing(false)
    @test compiled_ns == 0
end
