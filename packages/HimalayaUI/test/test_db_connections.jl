using Test, HTTP, SQLite, DBInterface, Tables, JSON3
import HimalayaUI
using HimalayaUI: open_db, bind_db!, current_db, db_connection_middleware

function _through_middleware(f, method::AbstractString, target::AbstractString = "/api/x")
    handler = db_connection_middleware(req -> f())
    handler(HTTP.Request(method, target))
end

@testset "GET requests read through a pooled reader connection, not the singleton" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            DBInterface.execute(db, "INSERT INTO users (username) VALUES ('alice')")
            seen = _through_middleware("GET") do
                reader = current_db()
                rows = Tables.rowtable(DBInterface.execute(reader, "SELECT username FROM users"))
                (reader, [String(r.username) for r in rows], current_db() === reader)
            end
            reader, names, stable = seen
            @test reader !== db
            @test names == ["alice"]
            @test stable
            @test current_db() === db

            again = _through_middleware(current_db, "GET")
            @test again === reader
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "a reader connection refuses writes" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            @test_throws SQLite.SQLiteException _through_middleware("GET") do
                DBInterface.execute(current_db(), "INSERT INTO users (username) VALUES ('mallory')")
            end
            @test isempty(Tables.rowtable(DBInterface.execute(db, "SELECT 1 FROM users WHERE username = 'mallory'")))
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "non-GET requests use the singleton and hold the write lock for the whole request" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            for method in ("POST", "PATCH", "DELETE", "PUT")
                conn, held = _through_middleware(method) do
                    (current_db(), islocked(HimalayaUI._DB_WRITE_LOCK) &&
                                   HimalayaUI._DB_WRITE_LOCK.locked_by === current_task())
                end
                @test conn === db
                @test held
            end
            @test !islocked(HimalayaUI._DB_WRITE_LOCK)
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "a GET does not wait for a writer holding the lock" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            inside = Channel{Nothing}(1)
            release = Channel{Nothing}(1)
            writer = @async _through_middleware("POST") do
                put!(inside, nothing)
                take!(release)
                nothing
            end
            take!(inside)
            reader = @async _through_middleware("GET") do
                Tables.rowtable(DBInterface.execute(current_db(), "SELECT COUNT(*) AS n FROM users"))[1].n
            end
            @test timedwait(() -> istaskdone(reader), 5.0) === :ok
            @test fetch(reader) == 0
            @test !istaskdone(writer)
            put!(release, nothing)
            wait(writer)
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "the event stream keeps the singleton without the write lock, and in-memory databases keep the singleton" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            conn, held = _through_middleware("GET", "/api/events") do
                (current_db(), islocked(HimalayaUI._DB_WRITE_LOCK))
            end
            @test conn === db
            @test !held
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
    mem = open_db(":memory:")
    bind_db!(mem)
    try
        @test _through_middleware(current_db, "GET") === mem
    finally
        HimalayaUI.unbind_db!()
        DBInterface.close!(mem)
    end
end

@testset "rebinding closes the old readers" begin
    mktempdir() do tmp
        first_db = open_db(joinpath(tmp, "a.db"))
        bind_db!(first_db)
        old_reader = _through_middleware(current_db, "GET")
        second_db = open_db(joinpath(tmp, "b.db"))
        bind_db!(second_db)
        try
            @test !isopen(old_reader)
            new_reader = _through_middleware(current_db, "GET")
            @test new_reader !== old_reader
            @test new_reader.file == second_db.file
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(first_db)
            DBInterface.close!(second_db)
        end
    end
end

@testset "background work gets a connection of its own, with foreign keys on" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        try
            conn, fk = HimalayaUI.with_background_db(db) do bg
                (bg, Tables.rowtable(DBInterface.execute(bg, "PRAGMA foreign_keys"))[1].foreign_keys)
            end
            @test conn !== db
            @test conn.file == db.file
            @test fk == 1
            @test !isopen(conn)
        finally
            DBInterface.close!(db)
        end
    end
    mem = open_db(":memory:")
    try
        @test HimalayaUI.with_background_db(identity, mem) === mem
    finally
        DBInterface.close!(mem)
    end
end

@testset "a GET never sees another request's uncommitted writes" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            inside = Channel{Nothing}(1)
            release = Channel{Nothing}(1)
            writer = @async _through_middleware("POST") do
                try
                    SQLite.transaction(db) do
                        DBInterface.execute(db, "INSERT INTO users (username) VALUES ('uncommitted')")
                        put!(inside, nothing)
                        take!(release)
                        error("roll back")
                    end
                catch
                end
            end
            take!(inside)
            names = _through_middleware("GET") do
                [String(r.username) for r in Tables.rowtable(DBInterface.execute(current_db(), "SELECT username FROM users"))]
            end
            @test names == String[]
            put!(release, nothing)
            wait(writer)
            @test isempty(Tables.rowtable(DBInterface.execute(db, "SELECT 1 FROM users")))
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "a write answers while a browser holds the event stream open" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir=tmp, analysis_dir=tmp)
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")
        try
            with_test_server(db) do port, base
                opened = Channel{Any}(1)
                stream = @async try
                    HTTP.open("GET", "$base/api/events"; readtimeout = 30) do io
                        startread(io)
                        put!(opened, io)
                        while !eof(io)
                            readavailable(io)
                        end
                    end
                catch
                end
                @test timedwait(() -> isready(opened), 10.0) === :ok
                write = @async HTTP.patch("$base/api/exposures/$e_id/status";
                    body = JSON3.write(Dict(:status => "rejected")),
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    readtimeout = 20, retry = false)
                @test timedwait(() -> istaskdone(write), 15.0) === :ok
                istaskdone(write) && @test fetch(write).status == 200
                isready(opened) && close(take!(opened).stream)
                timedwait(() -> istaskdone(stream), 10.0)
            end
        finally
            DBInterface.close!(db)
        end
    end
end

@testset "a reader left mid-query by one GET still sees later commits in the next" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            for name in ("a", "b", "c")
                DBInterface.execute(db, "INSERT INTO users (username) VALUES (?)", [name])
            end
            first_reader = _through_middleware("GET") do
                first(DBInterface.execute(current_db(), "SELECT username FROM users ORDER BY username"))
                current_db()
            end
            DBInterface.execute(db, "INSERT INTO users (username) VALUES ('d')")
            reader, n = _through_middleware("GET") do
                (current_db(), Tables.rowtable(DBInterface.execute(current_db(), "SELECT COUNT(*) AS n FROM users"))[1].n)
            end
            @test reader === first_reader
            @test n == 4
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "one GET reads one snapshot, even across a commit" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            counts = _through_middleware("GET") do
                count() = Tables.rowtable(DBInterface.execute(current_db(), "SELECT COUNT(*) AS n FROM users"))[1].n
                before = count()
                DBInterface.execute(db, "INSERT INTO users (username) VALUES ('late')")
                (before, count())
            end
            @test counts == (0, 0)
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end

@testset "HEAD requests read through the pool like GET" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        bind_db!(db)
        try
            conn, held = _through_middleware("HEAD") do
                (current_db(), islocked(HimalayaUI._DB_WRITE_LOCK))
            end
            @test conn !== db
            @test !held
        finally
            HimalayaUI.unbind_db!()
            DBInterface.close!(db)
        end
    end
end


@testset "the singleton waits for another writer instead of failing at once" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        try
            @test Tables.rowtable(DBInterface.execute(db, "PRAGMA busy_timeout"))[1].timeout == 5000
        finally
            DBInterface.close!(db)
        end
    end
end
