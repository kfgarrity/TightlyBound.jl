using Test
using TightlyBound
using Suppressor

#basic loading
@testset "basic crystal dftout tests" begin
    types=["Li"]
    pos=zeros((1,3))
    A=[ [10.0 0 0]; [0 10.0 0]; [ 0 0 10.0]]

    c=makecrys(A, pos, types)

    @test c.A == A
    @test c.coords == pos
    @test c.types == types

    c2=makecrys("test/POSCAR_testing")
    @test c2.types == ["Cl", "N","N","Sr", "Sr"]

    c3=makecrys("test/qe_testing.in")
    @test c3.types == ["Si", "C"]

    
    forces=ones((1,3))
    energy=0.0
    energy_smear = 0.0
    d=makedftout(c, energy, energy_smear, forces, zeros(3,3))
    @test energy == d.energy
    @test forces == d.forces

    d2=makedftout(A,pos,types, energy, energy_smear, forces, zeros(3,3))
    @test d.forces == d2.forces
    
    #supposed to throw an error
    @suppress begin
        forces2=zeros((2,3))
        @test_throws ErrorException("Error forces") makedftout(c, energy,energy_smear, forces2, zeros(3,3))
    end

    dft_out = TightlyBound.QE.loadXML("test/dimer.save/")

    convert_ryd_ha = 2.0

    @test convert_ryd_ha*-1.428445472333963e1 ≈ dft_out.energy
    @test convert_ryd_ha*-6.425031283667987e-2 ≈ dft_out.bandstruct.efermi
    @test convert_ryd_ha*-1.832485909355843e0 ≈ dft_out.bandstruct.eigs[1,1]
    
end


if false

    #run QE
    @testset "run QE" begin
        types=["Li"]
        pos=zeros((1,3))
        A=[ [6.0 0 0]; [0 6.0 0]; [ 0 0 6.0]]
        
        c=makecrys(A, pos, types)
        
        d = TightlyBound.DFT.runSCF(c, nprocs=4)
        
        @test -14.35143737 ≈ d.energy

    end
end
