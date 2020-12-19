###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module Ewald
"""
Electrostatics
"""

#include("Atomdata.jl")
using ..Atomdata:atoms

using LinearAlgebra
using ..CrystalMod:crystal
###using ..CalcTB:distances_etc_3bdy
using SpecialFunctions




function getU(types)

    U = zeros(length(types))
    for (c,t) in enumerate(types)
        U[c] = atoms[t].U
    end

    return U

end

function get_onsite(crys::crystal, U::Array{Float64,1})

    return diagm(U)
#    gamma_onsite = zeros(Float64, crys.nat, crys.nat)
#
#    for i = 1:crys.nat
#        gamma_onsite[i,i] = U[i]
#    end##
#
#    return gamma_onsite

end

function electrostatics_getgamma(crys::crystal;  kappa=missing, noU=false, onlyU=false, screening = 1.0)
#noU and onlyU are for testing purposes

#R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero = distances_etc_3bdy(crys,cutoff2X, 0.0)

    if ismissing(kappa)
#        kappa_default = 0.25
        kappa = estimate_best_kappa(crys.A)
    end
    kappa = Float64(kappa)

    if noU
        println("noU - FOR TESTING")
        U = zeros(Float64, crys.nat)
    else
        U = getU(crys.types)
        U = U * screening
    end

    starting_size_rspace = 1
    starting_size_kspace = 1

    gamma_onsiteU = get_onsite(crys, U)

    gamma_rs, gamma_U = real_space(crys, kappa, U, starting_size_rspace)
#    println("gamma_rs")
#    println(gamma_rs)
    gamma_k = k_space(crys, kappa, starting_size_kspace)
#    println("gamma_k")
#    println(gamma_k)

    #self
    gamma_self = zeros(Float64, crys.nat, crys.nat)
    for i = 1:crys.nat
        gamma_self[i,i] -= kappa / sqrt(pi) * 2.0
    end

    if false #for debugging
        println("gamma_rs")
        println(gamma_rs)
        println("gamma_k")
        println(gamma_k)
        println("gamma_self")
        println(gamma_self)
        println("gamma_onsiteU")
        println(gamma_onsiteU)
        println("gamma_U")
        println(gamma_U)
        println()
        println("only 1/r")
        println(gamma_rs + gamma_k + gamma_self)
    end

    #rydberg units
    e2 = 2.0

    gamma_tot = e2*(gamma_rs + gamma_k + gamma_self + gamma_U) + gamma_onsiteU

    if onlyU #for debugging
        gamma_tot = gamma_onsiteU
    end

    
    return gamma_tot

end

function real_space(crys::crystal, kappa::Float64, U::Array{Float64}, starting_size_rspace=2)
    
    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)

    gamma_U_tot = zeros(T, crys.nat, crys.nat)
    gamma_U_new = zeros(T, crys.nat, crys.nat)

    Uconst = zeros(T, crys.nat, crys.nat)

    if sum(abs.(U)) > 1e-5
        useU = true
        for i = 1:crys.nat
            Fi = sqrt( 8* log(2)/pi ) / (U[i]/2.0) #rydberg units, U in Ryd , we divide by 2, and multiply by e^2 layer. For Hartree formula, see koskinen comp mater sci 47 (2009) 237
            for j = 1:crys.nat
#                Uconst[i,j] = sqrt(pi/2 * (U[i]^2 * U[j]^2 / (U[i]^2 + U[j]^2)))

                Fj = sqrt( 8* log(2)/pi ) / (U[j]/2.0)
                Uconst[i,j] = sqrt(4 * log(2) / (Fi^2 + Fj^2))
            end
        end
    else
        useU = false
    end

    R = zeros(T, 1,3)
    
    coords_cart = crys.coords * crys.A
    
    converged = false
    R0 = false
    
    newcontr = 0.0
    newcontrU = 0.0
    
    for N = starting_size_rspace:20
        gamma_ij_new[:,:] .= 0.0
        gamma_U_new[:,:] .= 0.0
        for x = -N:N
            for y = -N:N
                for z = -N:N
                    if x != 0 || y != 0 || z != 0
                        R0 = false
                    else
                        R0 = true
                    end
                    R[1,:] = [x,y,z] 
                    R[1,:] = R * crys.A
                    if first_iter == true || (abs(x) > (N-1) || abs(y) > (N-1) || abs(z) > (N-1) )
                        for i = 1:crys.nat
                            for j = 1:crys.nat
                                r = sum((coords_cart[i,:] - coords_cart[j,:] - R[1,:]).^2)^0.5
                                if i != j ||  R0 == false
                                    gamma_ij_new[i,j] += erfc( kappa * r) / r

                                    if useU
                                        gamma_U_new[i,j] += -erfc( Uconst[i,j] * r) / r  #see eq 3 in prb 66 075212, or koskinen comp mater sci 47 (2009) 237
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        newcontr = sum(abs.(gamma_ij_new))
        newcontrU = sum(abs.(gamma_U_new))

#        println("new $N real_space $newcontr $newcontrU")

        if newcontr < 1e-7 && newcontrU < 1e-7
#            println("real_space YES converged $N : $newcontr")
            converged = true
            break
#        else
#            println("real_space NOT converged $N : $newcontr")
        end
        first_iter = false
        
        gamma_ij_tot += gamma_ij_new
        gamma_U_tot += gamma_U_new

    end

    if converged == false
        println("WARNING, real_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("real_space NOT converged  : $newcontr    $newcontrU")
        println(crys)
        
    end

    return gamma_ij_tot, gamma_U_tot

end


function estimate_best_kappa(A)

    a1 = sqrt(sum(A[1,:].^2))
    a2 = sqrt(sum(A[2,:].^2))
    a3 = sqrt(sum(A[3,:].^2))

    a = minimum([a1,a2,a3])

    B = inv(A)

    b1 = sqrt(sum(B[1,:].^2))
    b2 = sqrt(sum(B[2,:].^2))
    b3 = sqrt(sum(B[3,:].^2))
    
    b = minimum([b1,b2,b3])

    tot = Float64[]
    kappa_test = [0.00002 0.0001 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5 0.65 0.75 0.80 0.90 1.0 1.1 1.3 2.0 5.0 10.0 20.0 50.0 100.0]
    for kappa = kappa_test
        rs = erfc( kappa * a) / a * 20
        ks = exp(-b^2 / 4.0 / kappa^2  * (2 * pi)^2 ) / (2*pi*b)^2
        push!(tot, rs + ks)
#        println("$kappa rs $rs $ks $ks  tot ",rs + ks)
    end
    i = argmin(tot)
    kappa = kappa_test[i]
#    println("bestkappa :", kappa)
    
    return kappa
end

function k_space(crys::crystal, kappa, starting_size_kspace=2)

    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)
    K = zeros(T, 1,3)
    
    coords_cart = crys.coords * crys.A
    
    converged = false
    
    B = inv(crys.A)'
    vol = abs(det(crys.A))

    newcontr = 0.0
    
    for N = starting_size_kspace:25
        gamma_ij_new[:,:] .= 0.0
        for kx = -N:N
            for ky = -N:N
                for kz = -N:N
                    if kx == 0 && ky == 0 && kz == 0
                        continue
                    end
                    K[1,:] = [kx,ky,kz] 
                    K = K * B
                    k2 = sum(K.^2) * (2*pi)^2
                    factor_k = (2*pi) * exp(-k2 / 4.0 / kappa^2  ) / k2

                    if first_iter == true || abs(kx) > (N-1) || abs(ky) > (N-1) || abs(kz) > (N-1) 
                        for i = 1:crys.nat
                            for j = 1:crys.nat
                                kr = (K * (coords_cart[i,:] - coords_cart[j,:]) )[1]

                                exp_c = exp(2*pi*im*kr)
                                temp = real(exp_c )
                                gamma_ij_new[i,j] += factor_k * temp

                            end
                        end
                    end
                    
                end
            end
        end
        first_iter = false
        newcontr = sum(abs.(gamma_ij_new))
#        println("kspace $N newcontr $newcontr")
        if newcontr < 1e-7
#            println("k_space YES converged $N : $newcontr")
            converged = true
            break
#        else
#            println("k_space NOT converged $N : $newcontr")
        end

        gamma_ij_tot += 2.0*gamma_ij_new /vol
    end

    if converged == false
        println("WARNING, k_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("newcontr $newcontr")
        println(crys)
    end

    return gamma_ij_tot

end



end #end module
