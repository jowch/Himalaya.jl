@testset "precompile workload drives every hot route and leaves no global state" begin
    mktempdir() do tmp
        @test HimalayaUI.run_route_workload(tmp) === nothing
    end
    @test HimalayaUI._DB_REF[] === nothing
    @test isempty(HimalayaUI.SSE_SUBSCRIBERS[])
    @test isempty(HimalayaUI.OP_LOCKS)
    probe = Oxygen.internalrequest(HTTP.Request("GET", "/api/health"); serialize = true, catch_errors = true)
    @test probe.status == 404
end
