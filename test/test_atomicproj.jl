using TightlyBound
using Test
using Suppressor

TESTDIR=TightlyBound.TESTDIR


@testset "test atomicproj" begin

    
    
    @suppress begin

        units_old = TightlyBound.set_units()
        TightlyBound.set_units(both="atomic")


        filname = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/projham.xml.gz"
        tbc_ref = read_tb_crys(filname)

        fil = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/"
        dft = TightlyBound.QE.loadXML(fil*"/qe.save/")

        tbc, tbck, proj_warn = TightlyBound.AtomicProj.projwfx_workf(dft, directory=fil, nprocs=1, writefile=missing, writefilek=missing, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, only_kspace = false, screening = 1.0 )


        
        energy_tot, tbc_, conv_flag = scf_energy(tbc)
        
        @test (energy_tot - dft.atomize_energy) < 1e-4

        energy_tot_ref, tbc_, conv_flag = scf_energy(tbc_ref)
        @test (energy_tot_ref - dft.atomize_energy) < 1e-4

        TightlyBound.set_units(energy = units_old[1], length=units_old[2])
        
    end

    
end
