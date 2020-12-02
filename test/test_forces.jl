using Test
using TightlyBound
using Suppressor

#include("../includes_laguerre.jl")
#include("../Ewald.jl")

function test_force()

    @testset "testing force dimer" begin

        @suppress begin
            ft = open("./data_forces/fil_MgS_dimer", "r"); 
            dirst = readlines(ft); 
            close(ft); 

            #        println(dirst)

            function loaddata(dirs; scf=true)
                tbc_list  = []
                dft_list = []

                for t in dirs
                    #                println(t*"/qe.save")
                    dft = TightlyBound.QE.loadXML(t*"/qe.save")
                    tbc = []
                    tbc_scf = []
                    try
                        if scf
                            tbc_scf = TightlyBound.TB.read_tb_crys("projham_scf.xml.gz", directory=t)
                        else
                            tbc_scf = TightlyBound.TB.read_tb_crys("projham.xml.gz", directory=t)
                        end
                    catch
                        tbc = TightlyBound.AtomicProj.projwfx_workf(dft, directory=t, writefile="projham.xml", skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15)
                        if scf
                            tbc_scf = TightlyBound.SCF.remove_scf_from_tbc(tbc)
                            TightlyBound.TB.write_tb_crys(t*"/projham_scf.xml.gz", tbc_scf)
                        else
                            tbc_scf = tbc
                        end
                        
                    end
                    
                    push!(dft_list, dft)
                    push!(tbc_list, tbc_scf)
                end
                return     tbc_list, dft_list
            end

            for scf = [false true]
                
                tbc_list, dft_list = loaddata(dirst, scf=false);
                database_rec = TightlyBound.FitTB.do_fitting_recursive(tbc_list,dft_list = dft_list,  fit_threebody=false, fit_threebody_onsite=false);


                x = 4;
                smearing = 0.01;  
                en, f_cart,stress = TightlyBound.Force_Stress.get_energy_force_stress(tbc_list[x].crys, database_rec,   smearing = smearing);
                enFD, f_cartFD = TightlyBound.Force_Stress.finite_diff(tbc_list[x].crys, database_rec,1, 3,   smearing = smearing);

                @test abs(f_cartFD - f_cart[1,3]) < 1e-3
                #            println("SCF $scf TEST1 finite diff: ", f_cartFD , " autodiff:   ", f_cart[1,3], " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
                #            println("TEST1 dft ref ", dft_list[x].forces[1,3])

            end
            #        x = 3;
            #        smearing = 0.01;  
            #        en, f_cart,stress = Force_Stress.get_energy_force_stress(tbc_list[x].crys, database_rec,   smearing = smearing);
            #        enFD, f_cartFD = Force_Stress.finite_diff(tbc_list[x].crys, database_rec,1, 3,   smearing = smearing);

            #        @test abs(f_cartFD - f_cart[1,3]) < 1e-3
            #        println("TEST3 finite diff: ", f_cartFD , " autodiff:   ", f_cart[1,3], " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            #        println("TEST3 dft ref ", dft_list[x].forces[1,3])

        end
    end
end

function test_stress()

    @testset "testing force znse" begin
        @suppress begin
        
            ft = open("./data_forces/fil_MgS_znse", "r"); 
            dirst = readlines(ft); 
            close(ft); 

            #        println(dirst)

            function loaddata(dirs; scf=true)
                tbc_list  = []
                dft_list = []

                for t in dirs
                    dft = TightlyBound.QE.loadXML(t*"/qe.save")
                    tbc = []
                    tbc_scf = []
                    try
                        if scf
                            tbc_scf = TightlyBound.TB.read_tb_crys("projham_scf.xml", directory=t)
                        else
                            tbc_scf = TightlyBound.TB.read_tb_crys("projham.xml", directory=t)
                        end
                    catch
                        tbc = TightlyBound.AtomicProj.projwfx_workf(dft, directory=t, writefile="projham.xml", skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15)
                        if scf
                            tbc_scf = TightlyBound.SCF.remove_scf_from_tbc(tbc)
                            TightlyBound.TB.write_tb_crys(t*"/projham_scf.xml", tbc_scf)
                        else
                            tbc_scf = tbc
                        end
                        
                    end
                    
                    push!(dft_list, dft)
                    push!(tbc_list, tbc_scf)
                end
                return     tbc_list, dft_list
            end



            tbc_list, dft_list = loaddata(dirst, scf=false);

            database = TightlyBound.FitTB.do_fitting(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false);
            #        database = FitTB.do_fitting_recursive(tbc_list,dft_list = dft_list,  fit_threebody=true, fit_threebody_onsite=false);

            x = 1;
            smearing = 0.01;  
            en, f_cart, stress = TightlyBound.Force_Stress.get_energy_force_stress(tbc_list[x].crys, database,   smearing = smearing);

            enFD, f_cartFD = TightlyBound.Force_Stress.finite_diff(tbc_list[x].crys, database,1, 3,   smearing = smearing);

            #        println("TEST force finite diff: ", f_cartFD , " autodiff:   ", f_cart[1,3], " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            #        println("TEST dft ref ", dft_list[x].forces[1,3])
            @test abs(f_cartFD - f_cart[1,3]) < 1e-4


            x=1
            enFD, stressFD = TightlyBound.Force_Stress.finite_diff(tbc_list[x].crys, database,1, 1, stress_mode=true,  smearing = smearing);

            #        println("TEST stress11 finite diff: ", stressFD , " autodiff:   ", stress[1,1], " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            #        println("TEST dft ref ", dft_list[x].stress[1,1])
            @test abs(stressFD - stress[1,1]) < 1e-5


            x=1
            enFD, stressFD = TightlyBound.Force_Stress.finite_diff(tbc_list[x].crys, database,1, 2, stress_mode=true,  smearing = smearing);

            #        println("TEST stress12 finite diff: ", stressFD , " autodiff:   ", stress[1,2], " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            #        println("TEST dft ref ", dft_list[x].stress[1,2])
            @test abs(stressFD - stress[1,2]) < 1e-5

            #        println("done")
        end
    end
end


test_force()

#println("sleep ")
#sleep(3)

test_stress()
