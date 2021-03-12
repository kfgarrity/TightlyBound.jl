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

using ..TightlyBound:convert_energy
using ..TightlyBound:convert_dos
using ..TightlyBound:global_energy_units

"""
    function projection(tbc::tb_crys, vects, sk3, grid; ptype=missing)

Figures out the projections.
`ptype` can be `:atomic` or `:orbs` for atom projection or orbital projection (:s,:p,:d)
Default is to choose `:atomic` except for elemental systems.
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
        println("Projection type: ", ptype)
    end

    names = []
    PROJ = []       
    pwan=[]
    
    if ptype == :atomic || ptype == "atomic" || ptype == :atoms || ptype == "atoms" || ptype == :atom || ptype == "atom" || ptype == :Atomic || ptype == "Atomic"
        for ti in Set(tbc.crys.stypes)
            
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
            push!(pwan, length(proj_inds))
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
                push!(pwan, length(proj_inds))                
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

    return proj, names, pwan
    
    
end
            
        
    
    

"""
    function gaussian_dos(tbc::tb_crys; grid=missing, smearing=0.02, npts=500, proj_type=missing, do_display=true)

Simple Gaussian DOS, mostly for testing.

- `npts` is number of energies
- `proj_type` can be `"none"`, `"atomic"`, or `"orbital"`
- `do_display=false` will suppress the actual plot

The combination of `smearing` and `grid` are important to get converged results.

See also `dos`

`return energies, dos, projected_dos, pdos_names`
"""
function gaussian_dos(tbc::tb_crys; grid=missing, smearing=0.02, npts=300, proj_type=missing, do_display=true)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.4))
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
        proj, names, pwan =  projection(tbc, vects, sk3, grid, ptype=proj_type)
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

    energies = convert_energy(energies)
    dos = convert_dos(dos)
    pdos = convert_dos(pdos)
    
    plot_dos(energies, dos, pdos, names, do_display=do_display)


    
    return energies, dos, pdos, names

    
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


"""
    function dos(tbc::tb_crys; grid=missing, npts=300, proj_type=missing, do_display=true)

DOS, using tetrahedral integration

- `grid` is the k-point grid. Defaults to 1.6 times the default grid for energy integration.
- `npts` is number of energies
- `proj_type` can be `"none"`, `"atomic"`, or `"orbital"`. Defaults to `atomic` if more than one atom type.
- `do_display=false` will suppress plotting

`return energies, dos, projected_dos, pdos_names`
"""
function dos(tbc::tb_crys; grid=missing, npts=300, proj_type=missing, do_display=true)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.4))
        println("grid $grid")
    end


    etot, efermi, vals, vects,sk3 = calc_energy_fft(tbc, grid=grid, return_more_info=true)

    vals = vals .- efermi
    
    nk = size(vals)[1]
    
    vmin = minimum(vals)
    vmax = maximum(vals)
    r = vmax - vmin

    energies = collect(vmin - r*0.05 : r*1.1 / npts    : vmax + r*0.05 + 1e-7)
    
    dos = zeros(size(energies))
    dos_id = zeros( length(energies), nthreads())

    if ismissing(proj_type) ||  (proj_type != "none"   && proj_type != :none)
        do_proj=true
        #println("Projection")
        proj, names, pwan =  projection(tbc, vects, sk3, grid, ptype=proj_type)
        nproj = size(proj)[3]

        pdos = zeros(length(energies),nproj)
        pdos_id = zeros(length(energies),nproj, nthreads())

    else
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
    end

    #println("setup tetra")
    kpts, kpts_dict, tetra_map = setup_tetra(grid)

    #    e=zeros(4)

    ntetra = nk*6
    norm = Float64(1/ntetra)
    

    #    ex=zeros(4)
#    px=zeros(4)    
#    ktet = zeros(Int64, 4)

    
    eps = 1e-7

    #println("calc dos tetra")

    range = 1:length(energies)
    
    @threads for nt = 1: ntetra
        id = threadid()
        ex=zeros(4)
        px=zeros(4)    
        e=zeros(4, tbc.tb.nwan)
        p=zeros(4, tbc.tb.nwan, nproj)
        ktet = zeros(Int64, 4)
        wt = zeros(4)


        ktet[:] .= tetra_map[:, nt]
        
        e[:,:] .= vals[ktet, :]        

        if do_proj
            p[:,:, :] .= proj[ktet,:, :]
            f14 = 0.0
            f24 = 0.0
            f35 = 0.0
            G = 0.0
            f13 = 0.0
            f31 = 0.0 
            f14 = 0.0 
            f41 =  0.0
            f23 = 0.0
            f32 = 0.0
            f24 = 0.0
            f42 = 0.0
            
            f12 =  0.0
            f21 = 0.0
            

        end
        
        for nband = 1:tbc.tb.nwan

            ex[:] = e[:,nband]


            if do_proj==false  ####
                
                sort!(ex)
                


                inds = (ex[1] .< energies .< ex[4] )
                
                for c in range[inds]

                    en = energies[c]
                    
                    if en < ex[4] && en >= ex[3]
                        
                        dos_id[c, id] += 3.0 * (ex[4] - en)^2 / (ex[4] - ex[1] + eps) / (ex[4] - ex[2] + eps) / (ex[4] - ex[3] + eps)
                        
                    elseif en < ex[3] && en >= ex[2]
                        
                        dos_id[c, id] += 1.0 / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps) * (3.0 * (ex[2] - ex[1]) + 6.0 * (en - ex[2]) - 3.0 * (ex[3] - ex[1] + ex[4] - ex[2]) / (ex[3]-ex[2] + eps) / (ex[4] - ex[2] + eps) * (en -ex[2])^2)
                        
                    elseif en < ex[2] && en >= ex[1]
                        
                        dos_id[c, id] += 3.0 * (en-ex[1])^2 / (ex[2] - ex[1] + eps) / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps)
                        
                    end
                    
                end

            elseif do_proj == true ####

                
                perm = sortperm(ex)
                ex[:] = ex[perm]

                inds = (ex[1] .< energies .< ex[4] )
                
                for c in range[inds]
#                for c in 1:length(energies)
                    en = energies[c]

                    if en < ex[4] && en >= ex[1]
                    
                        if en < ex[4] && en >= ex[3]
                            
                            f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                            f24 = (en-ex[4])/(ex[2]-ex[4] - eps)
                            f34 = (en-ex[4])/(ex[3]-ex[4] - eps)
                            
                            G = 3.0 * f14 * f24 * f34 / (ex[4] - en + eps)
                            
                            @inbounds wt[1] = f14 / 3.0
                            @inbounds wt[2] = f24 / 3.0
                            @inbounds wt[3] = f34 / 3.0
                            @inbounds wt[4] = (3.0 - f14 - f24 - f34) / 3.0
                            
                        elseif en < ex[3] && en >= ex[2]
                            
                            f13 = (en-ex[3])/(ex[1]-ex[3] - eps)
                            f31 = 1.0 - f13
                            f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                            f41 = 1.0 - f14
                            f23 = (en-ex[3])/(ex[2]-ex[3] - eps)
                            f32 = 1.0 - f23
                            f24 = (en-ex[4])/(ex[2]-ex[4] - eps)
                            f42 = 1.0 - f24
                            
                            G   =  3.0 * (f23*f31 + f32*f24)
                            
                            @inbounds wt[1]  =  f14 / 3.0 + f13*f31*f23 / (G + eps)
                            @inbounds wt[2]  =  f23 / 3.0 + f24*f24*f32 / (G + eps)
                            @inbounds wt[3]  =  f32 / 3.0 + f31*f31*f23 / (G + eps)
                            @inbounds wt[4]  =  f41 / 3.0 + f42*f24*f32 / (G + eps)
                            
                            G   =  G / (ex[4]-ex[1] + eps)
                            
                            
                        elseif en < ex[2] && en >= ex[1]
                            
                            
                            f12 = (en-ex[2])/(ex[1]-ex[2] - eps)
                            f21 = 1.0 - f12
                            f13 = (en-ex[3])/(ex[1]-ex[3] - eps)
                            f31 = 1.0 - f13
                            f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                            f41 = 1.0 - f14
                            
                            G  =  3.0 * f21 * f31 * f41 / (en-ex[1] + eps)
                            
                            @inbounds wt[1] =  (f12 + f13 + f14) / 3.0
                            @inbounds wt[2] =  f21  / 3.0
                            @inbounds wt[3] =  f31  / 3.0
                            @inbounds wt[4] =  f41  / 3.0
                            
                            
                        end
                        
                        
                        
                        @inbounds dos_id[c, id] += sum(wt) * G
                        
                        for projind = 1:nproj
                            px[:] = @view p[perm,nband, projind]
                            @inbounds pdos_id[c,projind, id] += ( px' * wt )  * G
                        end
                        
                    end
                    
                end

            end #### if do_proj

        end

    end

    
    dos[:] = sum(dos_id, dims=2) * norm
    dos = max.(dos, 0.0)
    
    #println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )
    correct =  tbc.tb.nwan / ( sum(dos) * (energies[2]-energies[1]) ) 

    dos = dos * correct
        
    if nproj > 0
        pdos = max.(pdos, 0.0)
        pdos[:,:] = sum(pdos_id, dims=3) * norm
        pdos = pdos * correct

        #println("Int PDOS " , sum(dos, dims=1) * (energies[2]-energies[1]) )

        t = sum(pdos, dims=1) * (energies[2]-energies[1])
        fix = pwan' ./ t
        pdos = pdos .* fix
        pdos2 = zeros(size(pdos))
        
    else
        pdos2=missing
    end

    
    #simple smoothing
    dos2 = zeros(size(dos))


    for i = 3:(length(dos)-2)
        dos2[i] =  0.125*dos[i-2] + 0.25*dos[i-1] + 0.25*dos[i] + 0.25*dos[i+1] + 0.125*dos[i+2]
        if nproj>0
            #            pdos2[i,:] = 0.25*pdos[i-1,:] + 0.5*pdos[i,:] + 0.25*pdos[i+1,:]
            pdos2[i,:] = 0.125*pdos[i-2,:] + 0.25*pdos[i-1,:] + 0.25*pdos[i,:] + 0.25*pdos[i+1,:] + 0.125*pdos[i+2,:]            
        end
    end

    if nproj > 0
        pdos2 = max.(pdos2, 0.0)
        pdos2 = convert_dos(pdos2)
    end
    
    energies = convert_energy(energies)
    dos2 = convert_dos(dos2)

    plot_dos(energies, dos2, pdos2, names, do_display=do_display)

    
    return energies, dos2, pdos2, names
    
end

"""
    function plot_dos(energies, dos, pdos, names; filename=missing, do_display=true)

Does the actual DOS plotting, called by `dos` or `gaussian_dos`
"""
function plot_dos(energies, dos, pdos, names; filename=missing, do_display=true)

    if ismissing(names)
        plot(legend=false, grid=false, framestyle=:box)
    else 
        plot(legend=true, grid=false, framestyle=:box)
    end

    plot!(energies, dos, color="black", lw=4, label="Total", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=12)

    if !ismissing(names)
        colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
        for i in 1:size(pdos)[2]
            color = colors[i%7+1]
#            println("i $i $color", names[i])
            plot!(energies, pdos[:,i], color=color, lw=3, label=names[i])
        end
    end    

    if global_energy_units == "eV"
        ylabel!("DOS  ( 1 / eV )", guidefontsize="16")
        xlabel!("Energy - E_F ( eV )", guidefontsize="16")
    else
        ylabel!("DOS  ( 1 / Ryd. )", guidefontsize="16")
        xlabel!("Energy - E_F ( Ryd. )", guidefontsize="16")
    end

    xl1 = minimum(energies)
    xl2 = min(10.0, maximum(energies))

    xlims!(xl1, xl2)
    
    if do_display
        display(plot!([0,0], [0, maximum(dos) * 1.1], color="black", linestyle=:dash, label=""))
    else

    plot!([0,0], [0, maximum(dos) * 1.1], color="black", linestyle=:dash, label="")

    end        

    if !ismissing(filename)
        savefig(filename)
    end
    
end


end #end module
