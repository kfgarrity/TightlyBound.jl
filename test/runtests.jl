using TightlyBound
using Test

@testset "TightlyBound.jl" begin
    TightlyBound.f()
    TightlyBound.g()
    TightlyBound.h()
    TightlyBound.i()
    TightlyBound.j()
    @test 3 == 2+1
end
