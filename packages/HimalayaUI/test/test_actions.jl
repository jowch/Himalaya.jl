using Test, HTTP
using HimalayaUI

@testset "get_client_id" begin
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", ["X-Client-Id" => "abc-123"], UInt8[])
    ) == "abc-123"
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", Pair{String,String}[], UInt8[])
    ) === nothing
    @test HimalayaUI.get_client_id(
        HTTP.Request("POST", "/x", ["X-Client-Id" => ""], UInt8[])
    ) === nothing
end
