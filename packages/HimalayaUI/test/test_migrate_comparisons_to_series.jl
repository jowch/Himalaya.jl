using Test, SQLite, DBInterface, Tables, JSON3
using HimalayaUI: create_schema!, migrate_schema!, migrate_comparisons_to_series!,
                  MIGRATION_COMPARISONS_TO_SERIES

@testset "migrate_comparisons_to_series!" begin
    @testset "empty corpus: sentinel written, no series rows, idempotent" begin
        db = SQLite.DB()
        create_schema!(db)
        migrate_schema!(db)   # runs the migration once as part of the sequence

        # Sentinel row present after the full migrate_schema! run.
        sent = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES]))
        @test length(sent) == 1

        # No comparisons → no series.
        n_series = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM series"))[1].n
        @test n_series == 0

        # Calling the migration again is a gated no-op (does not throw, no dup sentinel).
        migrate_comparisons_to_series!(db)
        sent2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM schema_migrations WHERE name = ?",
            [MIGRATION_COMPARISONS_TO_SERIES]))[1].n
        @test sent2 == 1
        close(db)
    end
end
