module Relax

using ..SCF:scf_energy
using ..Force_Stress:get_energy_force_stress_fft
using ..Force_Stress:safe_mode_energy
using ..CrystalMod:get_grid
using ..CrystalMod:crystal
using ..CrystalMod:makecrys
using ..CalcTB:calc_tb_fast

using LinearAlgebra
using ..Force_Stress:inv_reshape_vec
using ..Force_Stress:reshape_vec

using Optim
using LineSearches

export relax_structure

"""
    function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true, conv_thr=2e-4)

Relax structure. Primary user function is relax_structure in TightlyBound.jl, which calls this one.
"""
function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true, conv_thr = 2e-4)

    if update_grid==false
        grid = get_grid(crys)
    end

    eden = missing #starts off missing

    #do this ahead of first iteration, to get memory in correct place
    tbc = calc_tb_fast(deepcopy(crys), database)
    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, conv_thr=1e-7)

    if error_flag
        println("warning error computing scf in relax_structure, zeroth iteration")
    end

    eden = deepcopy(tbc.eden)

    A0 = deepcopy(crys.A)
    strain = zeros(3,3)
    
    x0 = inv_reshape_vec(crys.coords, strain, crys.nat, strain_mode=true)
    crys_working = deepcopy(crys)

    fcall = 0
    firstiter = true

    nat = crys.nat

    energy_global = -99.0
    
    function fn(x)
#        println("CALL FN", x)
        coords, strain = reshape_vec(x, nat, strain_mode=true)

        A = A0 * (I(3) + strain)
        
        crys_working.coords = coords
        if  mode == "vc-relax"
            crys_working.A = A
        end

        if crys_working == tbc.crys
            #            println("SKIP")
            #we already have energy from calling grad, we don't need to call again.
            return energy_global
        end

        #        println("FN crys")
#        println(crys_working)
        
        tooshort, energy_short = safe_mode_energy(crys_working, database)

        if tooshort
            return energy_short
        end


        if firstiter
            firstiter = false
        else
            if crys_working != tbc.crys
#                println("yes calc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#                println(crys_working)
                tbc = calc_tb_fast(deepcopy(crys_working), database, verbose=false)
            else
#                println("nocalc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#                println(crys_working)
#                println(tbc.crys)
                
            end
#            println(crys_working)
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=1e-6)
            eden = deepcopy(tbcx.eden)
        end
#        println("fn $energy_tot fnffnfnfnffnfnfnfnfnfnfnfnfnfnfffffff")

        println("ENERGY $energy_tot $energy_global")
        
        energy_global=energy_tot

        return energy_tot

    end

    f_cart_global = []
    stress_global = []

    
    function grad(storage, x)
#        println("CALL GRAD", x)

#        println("typeof x ", typeof(x))
        fcall += 1

        coords, strain = reshape_vec(x, nat, strain_mode=true)

        A = A0 * (I(3) + strain)

        crys_working.coords = coords

        if  mode == "vc-relax"
            crys_working.A = A
        end

#        println("GRAD crys")
#        println(crys_working)

        
#        println("crys_working $fcall")
#        println(crys_working)

        tooshort, energy_short = safe_mode_energy(crys_working, database)
        
#        println("too short ", tooshort)

        if crys_working != tbc.crys && !tooshort
#            println("yes calc forces -----------------------------------------------------------------------------------------")
#            println(crys_working)
            tbc = calc_tb_fast(deepcopy(crys_working), database, verbose=false)
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=1e-6)

            eden = deepcopy(tbcx.eden)

            #        else
#            println("no calc forces ------------------------------------------------------------------------------------------")
#            println(crys_working)
#            println(tbc.crys)
#            println("__")
        end

        energy_global=energy_tot

        energy_tmp,  f_cart, stress =  get_energy_force_stress_fft(tbc, database; do_scf = false, smearing = smearing, grid = grid, vv=[VECTS, VALS, efermi] )
        #energy_tmp,  f_cart, stress =  get_energy_force_stress_fft(tbc, database; do_scf = true, smearing = smearing, grid = grid )

        
        f_cart_global = f_cart
        stress_global = stress
        
        fsum = sum((f_cart).^2)^0.5
        ssum= sum((stress).^2)^0.5

#        println("f_cart")
#        println(f_cart)

        println()
        println("FCALL $fcall en:  $energy_tot (Ryd)  fsum:  $fsum  ssum:  $ssum    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        println()
        
#        println("grad stress: ", stress)
        
#        stress_units = cell_force( crys_working.A, stress)
        #        stress_units = (stress * abs(det(crys_working.A)) *  inv(crys_working.A))'
        #stress_units = (stress * abs(det(crys_working.A)) )'        

#        println("grad stress units: ", -stress_units)

        #        f_crys = f_cart * inv(crys_working.A)
        f_crys = f_cart * crys_working.A

        g = inv_reshape_vec(-f_crys, -stress* abs(det(crys_working.A))  , crys_working.nat, strain_mode=true)

#        println("g ", g)
        
        if mode != "vc-relax" #don't need stress
            g[3*nat+1:end] = zeros(6)
        end

        storage[:] = g

        neaten(storage)

        #        println("ret storage $fcall")
#        println(storage)
    end


    println("starting vec ")
    println(x0)
    println()


#    opts = Optim.Options(g_tol = 7e-4,f_tol = 7e-4, x_tol = 7e-4,
#                         iterations = nsteps,
#                         store_trace = true,
#                         show_trace = false)

    
    opts = Optim.Options(g_tol = conv_thr,f_tol = conv_thr, x_tol = conv_thr,
                             iterations = nsteps,
                             store_trace = true,
                             show_trace = false)


    #res = optimize(fn,grad, x0, ConjugateGradient(), opts)


    #preconditioner, guess for inv Hess is based on the  metric. see qe bfgs_module.f90
    function init(x)
        num = 3*nat + 9

#        factor=0.1
        factor = 1.0
        
        P = zeros(eltype(x), num, num)

        A = crys.A
        g = A'*A
        ginv = inv(g)
        
        for a = 1:nat
            aa = (a-1)*3
            for i = 1:3
                for j = 1:3
                    P[aa+i ,aa+j] = g[i,j] * factor
                end
            end
        end

        vol = abs(det(A))
        
#        if mode == "vc-relax"
        for b = 1:3
            bb = 3 * nat + (b - 1)*3
            for i = 1:3
                for j = 1:3
                    P[bb+i,bb+j] = 0.04 * vol * ginv[i,j] * factor
                end
            end
        end
#    end
        
        return inv(P)
    end


    res = missing
    
    try    
              #        res = optimize(fn,grad, x0, BFGS(  initial_invH = init, linesearch=LineSearches.MoreThuente() ), opts)
              #res = optimize(fn,grad, x0, BFGS(initial_invH = init ), opts)

        res = optimize(fn,grad, x0, BFGS(linesearch=LineSearches.MoreThuente() ), opts)
        
    catch e
        if e isa InterruptException
            println("user interrupt")

            return tbc.crys
        else
            println("unknown error")
            return tbc.crys
            
        end
    end
    
    #res = optimize(fn,grad, x0, BFGS(  initial_invH = init, linesearch=LineSearches.HagerZhang() ), opts)    

    energy = res.minimum
    
    # bad : res = optimize(fn,grad, x0, ConjugateGradient(linesearch=LineSearches.BackTracking(order=3) ), opts)    
    #linesearch=LineSearches.BackTracking(order=2),
    #LineSearches.MoreThuente()
    #LineSearches.HagerZhang()
    #linesearch=LineSearches.BackTracking(order=2)
    #linesearch=LineSearches.StrongWolfe()
    

    println()
    println("res")
    println(res)

    minvec = Optim.minimizer(res)    
#    println("minvec")
#    println(minvec)

    coords, strain = reshape_vec(minvec, nat, strain_mode=true)

    A = A0 *( I(3) + strain)

    cfinal = makecrys(A, coords, crys.types, units="Bohr")
    return cfinal, tbc, energy, f_cart_global, stress_global

#    return res

end

"""
    function neaten(storage)

Detect/Enforce symmetries in forces/stress, tight tolerance. Deal with minor numerical issues causing symmetry breaking
"""
function neaten(storage)
#    println("n before ", storage)
    n = length(storage)
    for i in 1:n-6
        if abs(storage[i]) < 1e-8
            storage[i] = 0.0
        end
    end
        
    for i in n-2:n
        if abs(storage[i]) < 1e-7
            storage[i] = 0.0
        end
    end

    for i in 1:n
        for j in i+1:n
            if abs(storage[j] - storage[i]) < 1e-8
                storage[j] = storage[i]
            elseif abs(storage[j] + storage[i]) < 1e-8
                storage[j] = -storage[i]
            elseif abs(storage[j] + 2.0*storage[i]) < 1e-8
                storage[j] = -2.0*storage[i]
            elseif abs(storage[j] + 0.5*storage[i]) < 1e-8
                storage[j] = -0.5*storage[i]
            elseif abs(storage[j] - 2.0*storage[i]) < 1e-8
                storage[j] = 2.0*storage[i]
            elseif abs(storage[j] - 0.5*storage[i]) < 1e-8
                storage[j] = 0.5*storage[i]
            end
        end
    end
#    println("n after ", storage)
end
end #end module
