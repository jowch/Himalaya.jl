using Test

@testset "HimalayaUI" begin
    include("test_db.jl")
    include("test_datfile.jl")
    include("test_manifest.jl")
    include("test_pipeline.jl")
    include("test_json.jl")
    include("test_http.jl")
    include("test_routes_health.jl")
end
