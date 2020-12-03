using TightlyBound
using Test
using Suppressor


SRCDIR=TightlyBound.SRCDIR

@testset "testing prototypes" begin
    @suppress begin 
        include("$SRCDIR/Prototypes.jl")
        t =  oxidation_guess("Na", "Cl")
        tref = [["Na", "Cl", :core_binary],["Na", "Cl", :A1B1]]
        @test t == tref

        t =  oxidation_guess("Na", "Mg")
        tref = [["Na", "Mg", :core_binary],["Na", "Mg", :metals]]
        @test t == tref

    end
end
