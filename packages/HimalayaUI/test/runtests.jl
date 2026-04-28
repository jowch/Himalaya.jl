using Test

@testset "HimalayaUI" begin
    include("test_db.jl")
    include("test_datfile.jl")
    include("test_manifest.jl")
    include("test_pipeline.jl")
    include("test_json.jl")
    include("test_http.jl")
    include("test_routes_health.jl")
    include("test_routes_users.jl")
    include("test_routes_experiments.jl")
    include("test_routes_samples.jl")
    include("test_routes_exposures.jl")
    include("test_routes_image.jl")
    include("test_routes_status.jl")
    include("test_routes_exposures_filter.jl")
    include("test_routes_peaks.jl")
    include("test_routes_messages.jl")
    include("test_routes_trace.jl")
    include("test_routes_analysis.jl")
    include("test_routes_export.jl")
end
