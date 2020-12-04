###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module Force_Stress
"""
Scripts to calculate force and stress
"""

#calc for testing non-autodiff forces only 
#using Calculus

using LinearAlgebra
using ForwardDiff
using Optim
using ..CrystalMod:crystal
using ..CrystalMod:makecrys


using ..CalcTB:calc_tb_fast
using ..CalcTB:distances_etc_3bdy
using ..TB:calc_energy_charge_fft
using ..TB:tb_crys
using ..TB:types_energy
using ..TB:make_kgrid
using ..TB:get_dq
using ..TB:get_h1
using ..TB:ewald_energy
using ..Ewald:electrostatics_getgamma
using ..Ewald:estimate_best_kappa
using ..SCF:scf_energy

using ..BandTools:gaussian
using ..CrystalMod:get_grid

export get_energy_force_stress
export relax_structure

function get_energy_force_stress(crys::crystal, database; smearing = 0.01, grid = missing)

    println("crys")
#    tbc = []
    tbc = calc_tb_fast(crys, database)
    return get_energy_force_stress(tbc, database, do_scf=tbc.scf, grid = grid, smearing=smearing)
end

function get_energy_force_stress(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end
    
    kgrid, kweights = make_kgrid(grid)
    nk = size(kgrid)[1]

    tooshort, energy_tot = safe_mode_energy(tbc.crys, database)

    if !(tooshort)
        if !ismissing(vv)
            VECTS, VALS, efermi = vv
            energy_tot = 0.0
        else
        
            #prepare eigenvectors / values
            error_flag = false
            if do_scf
                energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-9)
            else
                energy_tot, efermi, e_den, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc, grid=grid, smearing=smearing)
            end
            if error_flag
                println("warning, trouble with eigenvectors/vals in initial step get_energy_force_stress")
            end
            
        end

        h1, dq = get_h1(tbc)
        
#        println("energy_tot $energy_tot")

        OCCS = gaussian.(VALS.-efermi, smearing)
    
    end


    ct = deepcopy(tbc.crys)
#=
    #reshape function
    function reshape_vec(x, nat)
        T=typeof(x[1])

        print("size x ", size(x))
        x_r = zeros(T, ct.nat, 3)
        for n = 1:ct.nat
            for j = 1:3
                x_r[n,j] = x[3*(n-1) + j]
            end
        end

        x_r_strain = zeros(T, 3,3)


        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = x[3*nat+4]
        x_r_strain[3,2] = x[3*nat+4]

        x_r_strain[1,3] = x[3*nat+5]
        x_r_strain[3,1] = x[3*nat+5]

        x_r_strain[1,2] = x[3*nat+6]
        x_r_strain[2,1] = x[3*nat+6]

        return x_r, x_r_strain
    end

    #reshape function
    function inv_reshape_vec(x, strain, nat)
        T=typeof(x[1])
        x_r = zeros(T, nat*3 + 6)
        for n = 1:nat
            for j = 1:3
                x_r[3*(n-1) + j] = x[n,j] 
            end
        end

        x_r[3*ct.nat + 1] = strain[1,1]
        x_r[3*ct.nat + 2] = strain[2,2]
        x_r[3*ct.nat + 3] = strain[3,3]
        x_r[3*ct.nat + 4] = strain[2,3]+strain[3,2]
        x_r[3*ct.nat + 5] = strain[1,3]+strain[3,1]
        x_r[3*ct.nat + 6] = strain[1,2]+strain[2,1]

        return x_r
    end
=#

    function f(x::Vector)
        T=typeof(x[1])

        x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
        
        A = ct.A * (I(3) + x_r_strain)
        #A = deepcopy(ct.A)

        crys_dual = makecrys( A , ct.coords + x_r, ct.types)


        #this deals with cases where the distances between atoms become very short, which can happen during relaxations
        #currently we just have an artificial repulsive force in this case
        if tooshort
            tooshort, energy_short = safe_mode_energy(crys_dual, database, var_type=T)
            return energy_short
        end


        
        if database["scf"] == true
            scf = true
            kappa = estimate_best_kappa(ct.A)
            gamma_dual = electrostatics_getgamma(crys_dual, kappa=kappa)
        else
            scf = false
            gamma_dual=zeros(T, ct.nat,ct.nat)
        end
        

        tbc_dual = calc_tb_fast(crys_dual, database; verbose=false, var_type=T, use_threebody=true, use_threebody_onsite=true, gamma=gamma_dual)

        nwan = tbc.tb.nwan

        hk = zeros(Complex{T}, nwan, nwan)
        hk0 = zeros(Complex{T}, nwan, nwan)
        hka = zeros(Complex{T}, nwan, nwan)
        sk = zeros(Complex{T}, nwan, nwan)

        twopi_i = -1.0im*2.0*pi

        VALS0 = zeros(T, nk, nwan)


        #analytic fourier transform
        for k = 1:nk
            #vect, vals, hk, sk, vals0 = Hk(hktemp, sktemp, tbc.tb, grid[k,:])


            hk0[:,:] .= 0.0
            sk[:,:] .= 0.0
            
            kmat = repeat(kgrid[k,:]', tbc.tb.nr,1)
            exp_ikr = exp.(twopi_i * sum(kmat .* tbc.tb.ind_arr,dims=2))

            for m in 1:tbc.tb.nwan
                 for n in 1:tbc.tb.nwan
                     hk0[m,n] = tbc_dual.tb.H[m,n,:]'*exp_ikr[:]
                     sk[m,n]  = tbc_dual.tb.S[m,n,:]'*exp_ikr[:]
                 end
            end
            hk0 = 0.5*(hk0 + hk0')
            sk = 0.5*(sk + sk')


            
            if scf
                hk0 = hk0 + h1 .* sk
            end

            for a = 1:nwan
                
                hka[:,:] = hk0 - ( VALS[k,a] )  * sk

                VALS0[k,a] += real.(VECTS[k,:,a]'*hka*VECTS[k,:,a])
            end
        end            


        energy0 = sum(OCCS .* VALS0) / nk * 2.0

        if scf
            eewald, pot = ewald_energy(crys_dual, gamma_dual, dq)
        else
            eewald = 0.0
        end


        etypes = types_energy(tbc.crys)
        

        return energy0 + etypes + eewald 

    end

#    x0 = inv_reshape_vec(ct.coords, ct.nat)

    g = ForwardDiff.gradient(f, zeros(3*ct.nat + 6)  )

    x, stress = reshape_vec(g, ct.nat)

    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)'

    stress = -stress / abs(det(ct.A))

    for i = 1:3
        for j = 1:3
            if abs(stress[i,j]) < 1e-12
                stress[i,j] = 0.0
            end
        end
    end

    return energy_tot,  f_cart, stress


end


#primarily for testing
function finite_diff(crys::crystal, database, ind1, ind2; stress_mode=false, step = 0.0002, smearing = 0.01, grid = missing)
    if ismissing(grid)
        grid = get_grid(crys)
    end


    tbc0 = calc_tb_fast(crys, database, verbose=false)

    energy_tot0, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc0, smearing=smearing, grid=grid)


    if stress_mode == false
        
        println("force mode")

        crys1 = deepcopy(crys)
        cart1 = crys1.coords * crys1.A
        cart1[ind1, ind2] += step

        crys1.coords = cart1 * inv(crys1.A)

        tbc1 = calc_tb_fast(crys1, database, verbose=false)

        energy_tot1, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc1, smearing=smearing, grid=grid)


        crys2 = deepcopy(crys)
        cart2 = crys2.coords * crys2.A
        cart2[ind1, ind2] -= step

        crys2.coords = cart2 * inv(crys2.A)

        tbc2 = calc_tb_fast(crys2, database, verbose=false)

        energy_tot2, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc2, smearing=smearing, grid=grid)
        
        force = - (energy_tot1 - energy_tot2) / (2 * step)

        return energy_tot0, force


    else
        println("stress mode")

        crys1 = deepcopy(crys)

        strain = zeros(3,3)
        strain[ind1,ind2] = step
        strain[ind2,ind1] = step
       
        crys1.A = crys1.A *(I(3) + strain)

        tbc1 = calc_tb_fast(crys1, database, verbose=false)

        energy_tot1, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc1, smearing=smearing, grid=grid)


        crys2 = deepcopy(crys)
        strain = zeros(3,3)
        strain[ind1,ind2] = -step
        strain[ind2,ind1] = -step
       
        crys2.A = crys2.A *(I(3) + strain)

        tbc2 = calc_tb_fast(crys2, database, verbose=false)

        energy_tot2, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc2, smearing=smearing, grid=grid)

        stress = -1.0* (energy_tot1 - energy_tot2) / (2 * step) / abs(det(crys.A))

        if ind1 != ind2
            stress = stress / 2.0
        end

        return energy_tot0, stress
    end
end






######################################################


function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true)

    if update_grid==false
        grid = get_grid(crys)
    end

    eden = missing #starts off missing

    #do this ahead of first iteration, to get memory in correct place
    tbc = calc_tb_fast(deepcopy(crys), database)
    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden)

    if error_flag
        println("warning error computing scf in relax_structure, zeroth iteration")
    end

    eden = deepcopy(tbc.eden)

    x0 = inv_reshape_vec(crys.coords, crys.A, crys.nat, strain_mode=false)
    crys_working = deepcopy(crys)

    fcall = 0
    firstiter = true

    nat = crys.nat

    function fn(x)

        coords, A = reshape_vec(x, nat)
        crys_working.coords = coords
        if  mode == "vc-relax"
            crys_working.A = A
        end

        tooshort, energy_short = safe_mode_energy(crys_working, database)
#        println("fn too short ", tooshort)
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
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false)
            eden = deepcopy(tbc.eden)
        end
        return energy_tot

    end

    function grad(storage, x)
#        println("grad x ", x)
#        println("typeof x ", typeof(x))
        fcall += 1

        coords, A = reshape_vec(x, nat)
        crys_working.coords = coords

        if  mode == "vc-relax"
            crys_working.A = A
        end

#        println("crys_working $fcall")
#        println(crys_working)

        tooshort, energy_short = safe_mode_energy(crys_working, database)
        
#        println("too short ", tooshort)

        if crys_working != tbc.crys && !tooshort
#            println("yes calc forces -----------------------------------------------------------------------------------------")
#            println(crys_working)
            tbc = calc_tb_fast(deepcopy(crys_working), database, verbose=false)
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false)
#        else
#            println("no calc forces ------------------------------------------------------------------------------------------")
#            println(crys_working)
#            println(tbc.crys)
#            println("__")
        end


        energy_tmp,  f_cart, stress =  get_energy_force_stress(tbc, database; do_scf = false, smearing = smearing, grid = grid, vv=[VECTS, VALS, efermi] )

        fsum = sum(abs.(f_cart))
        ssum= sum(abs.(stress))
        
        println("FCALL $fcall en:  $energy_tot  fsum:  $fsum  ssum:  $ssum    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

        stress_units = (stress * abs(det(crys_working.A)) *  inv(crys_working.A))'


        f_crys = f_cart * inv(crys_working.A)

        g = inv_reshape_vec(-f_crys, -stress_units, crys_working.nat, strain_mode=false)
        
        if mode != "vc-relax" #don't need stress
            g[3*nat+1:end] = zeros(9)
        end
        storage[:] = g
#        println("ret storage $fcall")
#        println(storage)
    end

#=
    function fn_fake(x)
        println("fn_fake ", x[1]*2)
        return sum(x.^2)
    end

    function grad_fake(storage, x)
        storage[:] =  2.0 * x

        println("grad_fake  ", x[1]*2)
        println("grad_fakex ", x)
        println("grad_fakes ", storage)
        return storage

    end
=#


#    opts = Optim.Options(g_tol = 1e-3,f_tol = 1e-3, x_tol = 1e-3,
#                         iterations = 10,
#                         store_trace = true,
#                         show_trace = false)


#    println("X000000000000000000, ", x0)
#    println()
#    println(fn(x0))
#    println("X000000000000000001, ", x0)
#    storage = zeros(size(x0))
#    println(grad(x0, storage))

#    println("storage")
#    println(storage)


    

#    println("starting vec ")
#    println(x0)
#    println()


    opts = Optim.Options(g_tol = 7e-4,f_tol = 7e-4, x_tol = 7e-4,
                         iterations = nsteps,
                         store_trace = true,
                         show_trace = false)


    res = optimize(fn,grad, x0, ConjugateGradient(), opts)

#for testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    storage = zeros(size(x0))
#    g  = grad(storage, x0)

#    println()

#    res = Calculus.gradient(fn, x0)

#    println("grad x0")
#    println(storage)
#    println("calc")
#    println(res)

#    println([storage res])
#
#    return [storage res]
#end for testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#    res = optimize(fn, x0, LBFGS())
#
#    res = optimize(fn,grad, x0, LBFGS())

#    res = optimize(fn, x0, ConjugateGradient())

#    res = optimize(fn,grad, x0, LBFGS(), opts)

#    res = optimize(fn, x0, LBFGS(), opts)

#    f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
#    function g!(storage, x)
 #       println("g $x");
 #       storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
 #       storage[2] = 200.0 * (x[2] - x[1]^2)
 #       println("g storage", storage)
 #   end

#    res = optimize(f, g!,  [0.0, 0.0], LBFGS())

    println()
    println("res")
    println(res)

    minvec = Optim.minimizer(res)    
#    println("minvec")
#    println(minvec)

    coords, A = reshape_vec(minvec, nat)
    cfinal = makecrys(A, coords, crys.types)
    return cfinal

#    return res

end

    
function reshape_vec(x, nat; strain_mode=false)
#    println("RESHAPEVEC RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRr ", strain_mode)

    T=typeof(x[1])
    
#    println("size x ", size(x))
    x_r = zeros(T, nat, 3)
    for n = 1:nat
        for j = 1:3
            x_r[n,j] = x[3*(n-1) + j]
        end
    end
    
    x_r_strain = zeros(T, 3,3)

    if length(x) == 3*nat+6 && strain_mode
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = 0.5*x[3*nat+4]
        x_r_strain[3,2] = 0.5*x[3*nat+4]

        x_r_strain[1,3] = 0.5*x[3*nat+5]
        x_r_strain[3,1] = 0.5*x[3*nat+5]

        x_r_strain[1,2] = 0.5*x[3*nat+6]
        x_r_strain[2,1] = 0.5*x[3*nat+6]
    elseif length(x) == 3*nat+6 
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = x[3*nat+4]
        x_r_strain[3,2] = x[3*nat+4]

        x_r_strain[1,3] = x[3*nat+5]
        x_r_strain[3,1] = x[3*nat+5]

        x_r_strain[1,2] = x[3*nat+6]
        x_r_strain[2,1] = x[3*nat+6]

    elseif length(x) == 3*nat+9
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[1,2] = x[3*nat+2]
        x_r_strain[1,3] = x[3*nat+3]
        x_r_strain[2,1] = x[3*nat+4]
        x_r_strain[2,2] = x[3*nat+5]
        x_r_strain[2,3] = x[3*nat+6]
        x_r_strain[3,1] = x[3*nat+7]
        x_r_strain[3,2] = x[3*nat+8]
        x_r_strain[3,3] = x[3*nat+9]
    else
        println("I'm confusing about the length reshape_vec $nat ", length(x) )
    end

    return x_r, x_r_strain
end

function inv_reshape_vec(x, strain, nat; strain_mode=true)
    T=typeof(x[1])
    if strain_mode
        x_r = zeros(T, nat*3 + 6)
    else
        x_r = zeros(T, nat*3 + 9)
    end

    for n = 1:nat
        for j = 1:3
            x_r[3*(n-1) + j] = x[n,j] 
        end
    end
    if strain_mode
        x_r[3*nat + 1] = strain[1,1]
        x_r[3*nat + 2] = strain[2,2]
        x_r[3*nat + 3] = strain[3,3]
        x_r[3*nat + 4] = (strain[2,3]+strain[3,2])
        x_r[3*nat + 5] = (strain[1,3]+strain[3,1])
        x_r[3*nat + 6] = (strain[1,2]+strain[2,1])
    else
        x_r[3*nat+1] = strain[1,1]  
        x_r[3*nat+2] = strain[1,2]  
        x_r[3*nat+3] = strain[1,3]  
        x_r[3*nat+4] = strain[2,1]  
        x_r[3*nat+5] = strain[2,2]  
        x_r[3*nat+6] = strain[2,3]  
        x_r[3*nat+7] = strain[3,1]  
        x_r[3*nat+8] = strain[3,2]  
        x_r[3*nat+9] = strain[3,3]  
    end

    return x_r
end


function safe_mode_energy(crys::crystal, database; var_type=Float64)

    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types = distances_etc_3bdy(crys,10.0, 0.0, var_type=var_type)
    nkeep = size(R_keep_ab)[1]
    
    energy = 1.0
    tooshort = false

    for a1 = 1:crys.nat
        t1 = crys.types[a1]
        for a2 = 1:crys.nat
            t2 = crys.types[a2]
            dmin = database[(t1,t2)].min_dist * 0.979
            for c = 1:nkeep
                cind = R_keep_ab[c,1]
                if dist_arr[a1,a2,cind,1] < dmin && dist_arr[a1,a2,cind,1] > 1e-7
                    tooshort = true
                    energy += (dist_arr[a1,a2,cind,1] - dmin)^2 + abs(dist_arr[a1,a2,cind,1] - dmin)
                    if var_type == Float64
                        println("WARNING, SAFE MODE $a1 $t1 $a2 $t2 $c ", dist_arr[a1,a2,cind,1])
                    end
                end
            end
        end
    end
    if tooshort
        println("WARNING, safe mode activated, minimum distances < fitting data * 0.98")
    end
    return tooshort, energy

end


end #end module
