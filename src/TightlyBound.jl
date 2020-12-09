"""
holds three body tight binding important stuff
"""
module TightlyBound

include("SetDir.jl")

include("Utility.jl")
include("BandTools.jl")
include("Atomic.jl")
include("Atomdata.jl")
include("Crystal.jl")

using .CrystalMod:crystal
using .CrystalMod:makecrys

export makecrys

include("Ewald.jl")
include("DFToutMod.jl")
using .DFToutMod:dftout
using .DFToutMod:makedftout
export dftout
export makedftout

include("TB.jl")

using .TB:tb_crys
using .TB:tb_crys_kspace

using .TB:Hk
using .TB:calc_bands
using .TB:plot_compare_tb
using .TB:plot_bandstr
using .TB:plot_compare_dft
using .TB:read_tb_crys

export Hk
export calc_bands
export plot_compare_tb
export plot_bandstr
export plot_compare_dft
export read_tb_crys

include("RunDFT.jl")


#include("RunWannier90.jl")
include("AtomicProj.jl")
include("CalcTB_laguerre.jl")
using .CalcTB:calc_tb_fast
export calc_tb_fast

include("SCF.jl")


include("FitTB_laguerre.jl")
include("Force_Stress.jl")


include("ManageDatabase.jl")

export scf_energy
export scf_energy_force_stress
export relax_structure


"""
    relax_structure(c::crystal; mode="vc-relax")

Find the lowest energy atomic configuration of crystal c.

...
# Arguments
- `c::crystal`: the structure to relax, only required argument
- `mode="vc-relax"`: Default (variable-cell relax) will relax structure and cell, anything else will relax structure only.
- `database=missing`: coefficent database, default is to use the pre-fit pbesol database
- `smearing=0.01`: smearing temperature (ryd), default = 0.01
- `grid=missing`: k-point grid, e.g. [10,10,10], default chosen automatically
- `nsteps=100`: maximum iterations
- `update_grid=true`: update automatic k-point grid during relaxation
...
"""
function relax_structure(c::crystal; database=missing, smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true)

    if ismissing(database)
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
    end
    if !ismissing(grid)
        update_grid = false
    end

    cfinal, tbc, energy, force, stress = Force_Stress.relax_structure(c, database, smearing=smearing, grid=grid, mode=mode, nsteps=nsteps, update_grid=update_grid)

    println("Final crystal")
    println(cfinal)

    println()
    println("Relax done")
#    println("Calculate final energy")
#
#    energy_tot, tbc, conv_flag = scf_energy(cfinal; database=database, smearing=0.01, grid = missing)
#    energy_tot, f_cart, stress = scf_energy_force_stress(tbc, database=database, smearing=smearing, grid=grid)
#
#    println("done with all relax")
#    println()

    return cfinal, tbc, energy, force, stress

end

function scf_energy_force_stress(c::crystal; database = missing, smearing = 0.01, grid = missing)
    
    energy_tot, tbc, conv_flag = scf_energy(c; database=database, smearing=smearing, grid = grid)

    if ismissing(database)
        database = ManageDatabase.database_cached
    end

    println()
    println("Calculate Force, Stress")
    
    energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress(tbc, database, do_scf=false, smearing=smearing, grid=grid)

    println("done")
    println("----")

    return energy_tot, f_cart, stress, tbc

end

function scf_energy_force_stress(tbc::tb_crys; database = missing, smearing = 0.01, grid = missing, do_scf=false)
    
    if ismissing(database)
        ManageDatabase.prepare_database(tbc.crys)
        database = ManageDatabase.database_cached
    end

    println()
    println("Calculate Force, Stress (no scf)")
    
    energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress(tbc, database, do_scf=false, smearing=smearing, grid=grid)

    println("done")
    println("----")

    return energy_tot, f_cart, stress, tbc

end



function scf_energy(d::dftout; database = Dict(), smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 50, mix = -1.0, mixing_mode=:pulay)

    return scf_energy(d.crys, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode)

end

function scf_energy(c::crystal; database = missing, smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 50, mix = -1.0, mixing_mode=:pulay)
    println()
    println("Begin scf_energy-------------")
    println()
    if ismissing(database)
        println("Load TB parameters from file")
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
        println()
    end

    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(c, database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix,  mixing_mode=mixing_mode)

    conv_flag = !error_flag
    if tbc.within_fit == false
        conv_flag = false
    end
    if conv_flag == true
        println("scf_energy success, done")
    else
        println("scf_energy error detected somewhere!!!!!!!!!!, done")
    end
    println()

    return energy_tot, tbc, conv_flag

end

function scf_energy(tbc::tb_crys; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 1e-5, iters = 75, mix = -1.0, mixing_mode=:pulay)

    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(tbc; smearing=smearing, grid = grid, e_den0 = e_den0, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode)

    conv_flag = !error_flag
    if tbc.within_fit == false
        conv_flag = false
    end
    if conv_flag == true
        println("success, done")
    else
        println("error detected!!!!!!!!!!, done")
    end
    println()

    return energy_tot, tbc, conv_flag

end




end #end module

