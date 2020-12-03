using TightlyBound
using Test
using Suppressor

for f in readdir("../examples")
    if occursin(".jl", f) && !occursin("~", f)
        @testset "example $f" begin
            @suppress begin 
                include("../examples/$f")
                @test 1 == 1
            end
        end        
    end
end




@testset "rm example pdfs generated" begin

    rm("mgs_compare_dft.pdf")
    rm("mgs_compare_self.pdf")
    rm("mgs_proj.pdf")

    @test 1 == 1
end
