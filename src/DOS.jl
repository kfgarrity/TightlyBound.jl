###include("Crystal1563.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module DOS
"""
DOS
"""

using LinearAlgebra
using Plots
using Base.Threads
using ..CrystalMod:get_grid
using ..TB:calc_energy_fft
using ..TB:tb_crys
using ..CrystalMod:crystal
using ..CrystalMod:orbital_index
using ..TB:summarize_orb



"""
    function simple_dos(tbc::tb_crys; grid=missing, smearing=0.01, npts=300)

Simple Gaussian DOS, mostly for testing.

return energies, dos
"""
function projection(tbc::tb_crys, vects, sk3, grid; ptype=missing)

    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)

    nk = size(vects)[1]
    
    if ismissing(ptype)
        if length(Set(tbc.crys.types)) != 1
            ptype=:atomic
        else
            ptype=:orbs
        end
        println("ptype ", ptype)
    end

    names = []
    PROJ = []       
    
    if ptype == :atomic || ptype == "atomic" || ptype == :atoms || ptype == "atoms"
        for ti in tbc.crys.stypes
            
            proj_inds = Int64[]
            for n = 1:tbc.tb.nwan
                at,t, orb = ind2orb[n]
                sorb = summarize_orb(orb)
                if t == ti
                    push!(proj_inds, n)
                end
            end

            push!(names, String(ti))
            push!(PROJ, proj_inds)
        end

    else
        
        for o in [:s, :p, :d, :f]
            proj_inds = Int64[]
            for n = 1:tbc.tb.nwan
                at,t, orb = ind2orb[n]
                sorb = summarize_orb(orb)
                if o == sorb
                    push!(proj_inds, n)
                end
            end
            if length(proj_inds) > 0
                push!(names, String(o))
                push!(PROJ, proj_inds)
            end

        end
    end

    temp = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    proj = zeros(nk, tbc.tb.nwan, length(PROJ))

    for (pind, proj_inds) in enumerate(PROJ)
        c=0
        for k1 = 1:grid[1]
            for k2 = 1:grid[2]
                for k3 = 1:grid[3]
                    c += 1
                    for p in proj_inds
                        for a = 1:tbc.tb.nwan
                            for b = 1:tbc.tb.nwan
                                temp[a,b] = vects[c, a,p]'*sk3[a, b, k1, k2, k3]* vects[c, b,p]
                            end
                        end
                        temp = temp + conj(temp)
                        proj[c,:, pind] += 0.5*sum(real(temp), dims=1)[:]
                    end
                end
            end
        end
                        
    end

    return proj, names
    
    
end
            
        
    
    

function simple_dos(tbc::tb_crys; grid=missing, smearing=0.01, npts=300, proj_type=missing)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.6))
        println("grid $grid")
    end

    etot, efermi, vals, vects, sk3 = calc_energy_fft(tbc, grid=grid, smearing=smearing, return_more_info=true)

    #prelim
    vals = vals .- efermi

    nk = size(vals)[1]
    
    vmin = minimum(vals)
    vmax = maximum(vals)
    r = vmax - vmin

    energies = collect(vmin - r*0.05 : r*1.1 / npts    : vmax + r*0.05 + 1e-7)
    
    dos = zeros(length(energies))

    if ismissing(proj_type) ||  proj_type != "none" 
        do_proj=true
        proj, names =  projection(tbc, vects, sk3, grid, ptype=proj_type)
        nproj = size(proj)[3]

        pdos = zeros(length(energies),nproj)

    else
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
    end
    


    
    for (c,e) in enumerate(energies)
        dos[c] = sum(exp.( -0.5 * (vals[:,:] .- e).^2 / smearing^2 ) )
    end
    dos = dos / smearing / (2.0*pi)^0.5 / nk
    
    if do_proj
        for i = 1:nproj
            for (c,e) in enumerate(energies)
                pdos[c, i] = sum(proj[:,:,i] .*  exp.( -0.5 * (vals[:,:] .- e).^2 / smearing^2 ) )
            end
        end
        pdos = pdos / smearing / (2.0*pi)^0.5 / nk
    end

    

    println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )
    
    return energies, dos, pdos

    
end

"""
    function setup_tetra(grid)

Setup simple tetrahedron method
based on tetra.f90 in QE. Simple tetra method.

"""
function setup_tetra(grid)

    ntetra = prod(grid) * 6

    tetra_map = zeros(Int64, 4, ntetra)

    kpts = zeros(Int64, prod(grid), 3)
    kpts_dict = Dict()
    
    c=0    
    for k1 = 0:grid[1]-1
        for k2 = 0:grid[2]-1
            for k3 = 0:grid[3]-1
                c+=1
                kpts[c,:] = [k1+1, k2+1, k3+1]
                kpts_dict[ [k1+1, k2+1, k3+1]] = c

            end
        end
    end

    c=1
    
    for k1 = 0:grid[1]-1
        for k2 = 0:grid[2]-1
            for k3 = 0:grid[3]-1
                
                k1a = (k1) % grid[1] + 1
                k2a = (k2) % grid[2] + 1
                k3a = (k3) % grid[3] + 1 

                k1b = (k1+1) % grid[1] + 1
                k2b = (k2+1) % grid[2] + 1
                k3b = (k3+1) % grid[3] + 1 

                #corners of cube (eight)
                n1 = kpts_dict[[k1a, k2a, k3a]]
                n2 = kpts_dict[[k1a, k2a, k3b]]
                n3 = kpts_dict[[k1a, k2b, k3a]]
                n4 = kpts_dict[[k1a, k2b, k3b]]
                n5 = kpts_dict[[k1b, k2a, k3a]]
                n6 = kpts_dict[[k1b, k2a, k3b]]
                n7 = kpts_dict[[k1b, k2b, k3a]]
                n8 = kpts_dict[[k1b, k2b, k3b]]

                #six tetra per cube
                tetra_map[1, c] = n1
                tetra_map[2, c] = n2
                tetra_map[3, c] = n3
                tetra_map[4, c] = n6

                tetra_map[1, c+1] = n2
                tetra_map[2, c+1] = n3
                tetra_map[3, c+1] = n4
                tetra_map[4, c+1] = n6
                
                tetra_map[1, c+2] = n1
                tetra_map[2, c+2] = n3
                tetra_map[3, c+2] = n5
                tetra_map[4, c+2] = n6
                
                tetra_map[1, c+3] = n3
                tetra_map[2, c+3] = n4
                tetra_map[3, c+3] = n6
                tetra_map[4, c+3] = n8
                
                tetra_map[1, c+4] = n3
                tetra_map[2, c+4] = n6
                tetra_map[3, c+4] = n7
                tetra_map[4, c+4] = n8
                
                tetra_map[1, c+5] = n3
                tetra_map[2, c+5] = n5
                tetra_map[3, c+5] = n6
                tetra_map[4, c+5] = n7
                
                c+=6

            end
        end
    end

    return kpts, kpts_dict, tetra_map

end


function dos(tbc::tb_crys; grid=missing, npts=300, proj_type=missing)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.6))
        println("grid $grid")
    end

    @time etot, efermi, vals, vects,sk3 = calc_energy_fft(tbc, grid=grid, return_more_info=true)

    vals = vals .- efermi
    
    nk = size(vals)[1]
    
    vmin = minimum(vals)
    vmax = maximum(vals)
    r = vmax - vmin

    energies = collect(vmin - r*0.05 : r*1.1 / npts    : vmax + r*0.05 + 1e-7)
    
    dos = zeros(size(energies))


    if ismissing(proj_type) ||  proj_type != "none" 
        do_proj=true
        proj, names =  projection(tbc, vects, sk3, grid, ptype=proj_type)
        nproj = size(proj)[3]

        pdos = zeros(length(energies),nproj)

    else
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
    end

    
    @time kpts, kpts_dict, tetra_map = setup_tetra(grid)

#    e=zeros(4)

    ntetra = nk*6
    norm = Float64(1/ntetra)
    
    e=zeros(4, tbc.tb.nwan)
    p=zeros(4, tbc.tb.nwan)
    ex=zeros(4)
    ktet = zeros(Int64, 4)

    eps = 1e-7
    
    #    @threads for (c,en) in enumerate(energies)
    @time for nt = 1: ntetra
        ktet[:] .= tetra_map[:, nt]
        
        e[:,:] .= vals[ktet, :]        

        if do_proj
            p[:,:] .= proj[ktet,:]
        end
        
        for nband = 1:tbc.tb.nwan
            ex[:] = e[:,nband]

            if do_proj
                
            else
                sort!(ex)
            end

#            println(ex)
            for c in 1:length(energies)
                en = energies[c]
                
                if en < ex[4] && en >= ex[3]
                    
                    dos[c] += 3.0 * (ex[4] - en)^2 / (ex[4] - ex[1] + eps) / (ex[4] - ex[2] + eps) / (ex[4] - ex[3] + eps)
                    
                elseif en < ex[3] && en >= ex[2]
                    
                    dos[c] += 1.0 / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps) * (3.0 * (ex[2] - ex[1]) + 6.0 * (en - ex[2]) - 3.0 * (ex[3] - ex[1] + ex[4] - ex[2]) / (ex[3]-ex[2] + eps) / (ex[4] - ex[2] + eps) * (en -ex[2])^2)
                        
                elseif en < ex[2] && en >= ex[1]
                    
                    dos[c] += 3.0 * (en-ex[1])^2 / (ex[2] - ex[1] + eps) / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps)
                    
                end
                    
            end

                
        end
        
    end
    
    dos = dos * norm

    println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )

    
    return energies, dos

end

end #end module
