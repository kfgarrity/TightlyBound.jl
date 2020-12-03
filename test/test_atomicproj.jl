using TightlyBound
using Test
using Suppressor

@testset "test atomicproj" begin

    @suppress begin
        filname = "../test/data_forces/znse.in_vnscf_vol_2/projham.xml.gz"
        tbc_ref = read_tb_crys(filname)

        fil = "../test/data_forces/znse.in_vnscf_vol_2/"
        dft = TightlyBound.QE.loadXML(fil*"/qe.save/")

        tbc, tbck, proj_warn = TightlyBound.AtomicProj.projwfx_workf(dft, directory=fil, nprocs=1, writefile=missing, writefilek=missing, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, only_kspace = false, screening = 0.9 )
        
        energy_tot, tbc_, conv_flag = scf_energy(tbc)
        
        @test (energy_tot - dft.atomize_energy) < 1e-5

        energy_tot_ref, tbc_, conv_flag = scf_energy(tbc_ref)
        @test (energy_tot_ref - dft.atomize_energy) < 1e-5

    end
end
