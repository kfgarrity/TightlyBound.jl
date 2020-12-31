
module CrystalMod


using LinearAlgebra
using Printf
using ..Atomdata:atoms
using GZip

export crystal
export makecrys
export generate_supercell
export generate_random
export write_poscar
export write_efs

#holds data
mutable struct crystal{T}

    A::Array{T,2}
    coords::Array{T,2}
    types::Array{String,1}
    nat::Int64

end




#printing
Base.show(io::IO, c::crystal) = begin

    for i in 1:3
        @printf(io, "A%.1i=     %.5f  %.5f  %.5f\n", i, c.A[i,1], c.A[i,2], c.A[i,3])
    end
    @printf(io, "\n")
    for i in 1:c.nat
        @printf(io, "%-3s  %.5f  %.5f  %.5f\n", c.types[i], c.coords[i,1], c.coords[i,2], c.coords[i,3])        
    end
    
end   


Base.:+(c1::crystal, c2::crystal) = begin

    if c1.nat != c2.nat
        println("nat not same ", c1.nat, c2.nat)
        return c1
    end
    A = c1.A .+ c2.A
    coords = c1.coords .+ c2.coords

    types = deepcopy(c1.types)
    
    return makecrys(A, coords, types)
end   

Base.:(==)(c1::crystal, c2::crystal) = begin

    if c1.nat != c2.nat
        return false
    end
    
    if c1.types != c2.types
        return false
    end
    
    if sum(abs.(c1.A - c2.A)) > 1e-11
        return false
    end
    
    if sum(abs.(c1.coords - c2.coords)) > 1e-11
        return false
    end
    
    return true

end   


Base.:*(a::Real, c::crystal) = begin

    A = a * c.A 
    coords = deepcopy( c.coords)

    types = deepcopy(c.types)
    
    return makecrys(A, coords, types)
end   

Base.:*(c::crystal,a::Real) = begin
    return a*c
end

Base.:*(a::Array{Int}, c::crystal) = begin

    A_new = zeros(3,3)
    A_new[1,:] = c.A[1,:] *a[1]
    A_new[2,:] = c.A[2,:] *a[2]
    A_new[3,:] = c.A[3,:] *a[3]

    n = prod(a)
    coords_new = zeros(n*c.nat, 3)
    types_new = String[]

    counter = 0
    for x = 1:a[1]
        for y = 1:a[2]
            for z = 1:a[3]
                for i = 1:c.nat
                    counter +=1
                    coords_new[counter,:] = (c.coords[i,:] + [x-1, y-1, z-1]) ./ [a[1], a[2], a[3]]
                    push!(types_new, c.types[i])
                end
            end
        end
    end

    return makecrys(A_new, coords_new, types_new)
end   

Base.:*(c::crystal,  a::Array{Int} ) = begin
    return a*c
end


#setup crystal
"""
    makecrys(A,coords,types)

Return a crystal object from 3×3 lattice A (Bohr), nat × 3 coords in crystal units, and nat element strings

Note: also export-ed directly from TightlyBound for convenience


```julia-repl
julia> makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0.0 0.0 0.0], ["H"])
A1=     10.00000  0.00000  0.00000
A2=     0.00000  10.00000  0.00000
A3=     0.00000  0.00000  10.00000

H    0.00000  0.00000  0.00000
```
"""
function makecrys(A,coords,types)
    T = typeof(coords[1,1])

    #entering integer coords or A messes everything up later
    if T == Int64
        coords = Float64.(coords)
        T = typeof(coords[1,1])
    end
    if typeof(A[1,1]) == Int64
        A = Float64.(A)
    end

    if length(size(types)) == 2  # if are accidently given a 1 x nat types array 
        types = types[:]
    end

    nat = length(types)

    if nat != size(coords,1)
        println(types)
        println(coords)
        error("Error initilzing crystal, nat doesn't match", nat, size(coords,1)  )
        
        return 
    end

    if (3,3) != size(A)
        println(A)
        error("Error initilzing crystal, A wrong", size(A))
        
        return 
    end

    if 3 != size(coords,2)
        println(coords)
        error("Error initilzing crystal, coords wrong", size(coords))
        
        return 
    end
    
    A = convert(Array{T,2}, A)
    coords = convert(Array{T,2}, coords)

    if size(coords)[1] != size(types)[1]
        error("Error making crys, types and coords sizes don't match")
        return 
    end

    nat = size(coords)[1]
    return crystal{T}(A, coords, types, nat)
end


"""
    makecrys(filename::String)

Read filename, return crystal object.
File can be POSCAR, or simple quantum espresso inputfile


"""
function makecrys(filename::String)

    if !isfile(filename)
        println("error - tried to open $filename to load crystal, file does not exist")
    end

    f = gzopen(filename, "r")
    lines = readlines(f)
    close(f)
    
    c = makecrys(lines)
    return c
    
end


"""
    makecrys(lines::Array{String,1})

Read string array, return crystal object.
string array can be POSCAR, or simple quantum espresso inputfile
"""
function makecrys(lines::Array{String,1})

    intype="POSCAR"
    for line in lines
        sp = split(line)
        if size(sp)[1] >= 1 && sp[1] == "&control"
            intype="QE"
        end
    end

    if intype == "POSCAR"

        A, coords, types = parsePOSCAR(lines)

    elseif intype == "QE"

        A, coords, types = parseQEinput(lines)
        
    end




    
    c = makecrys(A, coords,types)
    return c
    
    end

function parseARRfloat(sp)

    return map(x->parse(Float64,x),sp)

end

function parseARRint(sp)

    return map(x->parse(Int64,x),sp)

end

function parsePOSCAR(lines)

    
    title = lines[1]

    
    a = parse(Float64, split(lines[2])[1])
    
    A = zeros(3,3)
    A[1,:] = parseARRfloat(split(lines[3]))
    A[2,:] = parseARRfloat(split(lines[4]))
    A[3,:] = parseARRfloat(split(lines[5]))

    A = A * a / 0.529177 

    thetypes = split(lines[6])



    thenumbers = parseARRint(split(lines[7]))



    
    if size(thetypes)[1] != size(thenumbers)[1]

        return -1,-1,-1
    end
    
    nat = sum(thenumbers)


    dc = split(lines[8])[1][1]
    cart = false
    if dc == 'c' || dc == 'C' 
        cart = true
    end

    
    types = String[]
    for (t,n) in zip(thetypes,thenumbers)

        for i = 1:n
            push!(types, t)
        end
    end
    
    coords = zeros(nat,3)
    for i = 1:nat

        coords[i,:] = parseARRfloat(split(lines[8+i]))
    end

    if cart
        coords = (coords / 0.529177 )* inv(A)
    end

    return A, coords, types
end
    
function parseQEinput(lines)

    posind = -1
    Aind = -1
    types = String[]
    A=zeros(3,3)
    coords = Float64[]
    nat = -1
    units = 1.0
    for line in lines
        liner = replace(replace(line, "," => " "), "=" => " ")

        sp = split(liner)
        
        if size(sp)[1] >= 1

            if posind >= 0
                posind += 1
                coords[posind,:] = parseARRfloat(sp[2:4])
                push!(types, sp[1])
                if posind >= nat
                    posind = -1
                end
            end

            if Aind >= 0
                Aind += 1
                A[Aind,:] = parseARRfloat(sp[1:3])
                if Aind >= 3
                    Aind = -1
                end
            end
            
            if sp[1] == "nat"
                nat = parse(Int64, sp[2])
                coords = zeros(nat,3)
                
            elseif  sp[1] == "ATOMIC_POSITIONS"
                posind = 0
                if len(sp) > 1
                    if sp[2] != "crystal" && sp[2] != "(crystal)"
                        println("warning, only crystal coords supported !!!!!!!!!!!!!!!!!!")
                    end
                end
            elseif sp[1] == "CELL_PARAMETERS"
                Aind = 0
                if len(sp) > 1
                    if sp[2] == "angstrom" || sp[2] == "Angstrom" || sp[2] == "(angstrom)"
                        units = 0.529177
                    else 
                        println("warning alat or other CELL_PARAMETERS not supported !!!!!!!!!!!!!!!!!!!!!!!")
                    end
                end
            end
                
                   
        end 
    end

    A = A * units

    

    return A, coords, types

    
end

function generate_supercell(crys, cell)

    if typeof(cell) == Int
        cell = [cell,cell,cell]
    end

    A = copy(crys.A)
    A[1,:] *= cell[1]
    A[2,:] *= cell[2]
    A[3,:] *= cell[3]

    cells = prod(cell)

    println("cells ", cell, " " , cells)
    
    coords = zeros(crys.nat*cells, 3)

    types=String[]
    
    c = 0
    for i in 0.0:cell[1]-1
        for j in 0.0:cell[2]-1
            for k in 0.0:cell[3]-1
                for at in 1:crys.nat
                    c+= 1

                    coords[c,:] = (crys.coords[at,:] + [i, j, k]) ./ cell'
                    push!(types, crys.types[at])
                end
            end
        end
    end

    csuper = makecrys(A, coords, types)
    return csuper
end

function generate_random(crys, amag, strain_mag)


    st = (rand(3,3) .- 0.5) * strain_mag
    st = (st + st')/2.0
    
    A = crys.A * (I + st)
    coords_real = (crys.coords * A) + (rand(crys.nat, 3) .- 0.5)*amag
    coords = coords_real * inv(A)

    return makecrys(A, coords, copy(crys.types))
    
end

function write_poscar(crys, filename)

    fil = open(filename, "w")
    write(fil, "title kfg\n")
    write(fil, "1.0000000\n")

    Aang = crys.A * 0.529177
    
    write(fil, "$(Aang[1,1])  $(Aang[1,2])  $(Aang[1,3]) \n")
    write(fil, "$(Aang[2,1])  $(Aang[2,2])  $(Aang[2,3]) \n")
    write(fil, "$(Aang[3,1])  $(Aang[3,2])  $(Aang[3,3]) \n")

    tstr = ""
    numstr = ""
    for t in crys.types
        tstr = tstr * t * " " 
        numstr = numstr * "1 "
    end
    write(fil,tstr*'\n')
    write(fil,numstr*'\n')
 
    write(fil,"Direct\n")
    for at in 1:crys.nat
        write(fil, "$(crys.coords[at,1])  $(crys.coords[at,2]) $(crys.coords[at,3])\n")
    end
    close(fil)
    
    return 0
    
end


function write_efs(crys, energy, forces, stress, filename)

    fil = open(filename, "w")

    write(fil,"Program PWSCF fake\n")
    write(fil,"number of atoms/cell = $(crys.nat) \n")
#    ntypes = pos.shape[0]
    ntypes = 1
    write(fil,"number of types = 1\n")
    write(fil,"celldm(1)= 1.00\n")

    write(fil,"a(1) = ( $(crys.A[1,1]) $(crys.A[1,2]) $(crys.A[1,3]) ) \n")
    write(fil,"a(2) = ( $(crys.A[2,1]) $(crys.A[2,2]) $(crys.A[2,3]) ) \n")
    write(fil,"a(3) = ( $(crys.A[3,1]) $(crys.A[3,2]) $(crys.A[3,3]) ) \n")

    write(fil,"     site n.     atom                  positions (cryst. coord.)\n")
    for na in 1:crys.nat
        write(fil,"         1       $(crys.types[na])   tau(   $na) = (  $(crys.coords[na,1]) $(crys.coords[na,2]) $(crys.coords[na,3])  )\n")
    end
    write(fil,"\n")

    write(fil,"!    total energy              =     $energy Ry\n")

                   

    write(fil,"     Forces acting on atoms (Ry/au):\n")
    for na in 1:crys.nat
        write(fil,"     atom    $na type  1   force =   $(forces[na,1]) $(forces[na,2]) $(forces[na,3]) \n")
    end
        
    write(fil,"The non-local\n")

    write(fil,"          total   stress  (Ry/bohr**3)                   (kbar)     P=  ???\n")
    write(fil,"$(stress[1,1]) \t $(stress[1,2]) \t  $(stress[1,3])   0 0 0 \n")
    write(fil,"$(stress[2,1]) \t $(stress[2,2]) \t  $(stress[2,3])   0 0 0 \n")
    write(fil,"$(stress[3,1]) \t $(stress[3,2]) \t  $(stress[3,3])   0 0 0 \n")


    write(fil,"JOB DONE.\n")

    
    
    close(fil)

end

function get_grid(c, kden=55.0)

    B = transpose(inv(c.A))
#    kden = 55.0 

#    b1 = 1.0/norm(c.A[1,:])
#    b2 = 1.0/norm(c.A[2,:])
#    b3 = 1.0/norm(c.A[3,:])

    b1 = norm(B[1,:])
    b2 = norm(B[2,:])
    b3 = norm(B[3,:])


    k1 = convert(Int, round(kden * b1))
    k2 = convert(Int, round(kden * b2))
    k3 = convert(Int, round(kden * b3))

    if k1%2 == 1
        k1 = k1 + 1
    end
    if k2%2 == 1
        k2 = k2 + 1
    end
    if k3%2 == 1
        k3 = k3 + 1
    end
    k1 = max(k1, 2)
    k2 = max(k2, 2)
    k3 = max(k3, 2)

    k1 = min(k1, 14)
    k2 = min(k2, 14)
    k3 = min(k3, 14)

    
    kpoints = [k1, k2, k3]
#    println("get grid ", kpoints)
    return kpoints
    
#    R = [0,0,0]
#    for i = 1:3
#        R[i] = Int64(round(20.0/sum(c.A[1,:].^2)^0.5))
#        R[i] = max(R[i], 2)
#    end
#    println("get grid ", R)
#    return R

end    

function orbital_index(c::crystal)

    ind2orb = Dict()
    orb2ind = Dict()

    atomtypes = []
    ntypes=0
    for t in c.types
        if !(t in atomtypes)
            ntypes += 1
            push!(atomtypes, t)
        end
    end

    wan_counter = 1

    nsemi = 0
    nval = 0
    nwan = 0
    ind = 0

    etotal = 0.0
    
#    for at in atomtypes
#        for (i, t) in enumerate(c.types)
#            if t == at
    for (i,t) in enumerate(c.types)
        atom = atoms[t]
        nsemi += atom.nsemicore
        nval += atom.nval
        nwan += atom.nwan

        etotal += atom.total_energy
        
        for o in atom.orbitals
            if o == :s
                ind += 1
                ind2orb[ind] = [i,t, :s]
            elseif o == :p
                ind += 1
                ind2orb[ind] = [i,t, :pz]
                ind += 1
                ind2orb[ind] = [i,t, :px]
                ind += 1
                ind2orb[ind] = [i,t, :py]
            elseif o == :d
                ind += 1
                ind2orb[ind] = [i,t, :dz2]
                ind += 1
                ind2orb[ind] = [i,t, :dxz]
                ind += 1
                ind2orb[ind] = [i,t, :dyz]
                ind += 1
                ind2orb[ind] = [i,t, :dx2_y2]
                ind += 1
                ind2orb[ind] = [i,t, :dxy]
            elseif o == :f
                ind += 1
                ind2orb[ind] = [i,t, :fz3]
                ind += 1
                ind2orb[ind] = [i,t, :fxz2]
                ind += 1
                ind2orb[ind] = [i,t, :fyz2]
                ind += 1
                ind2orb[ind] = [i,t, :fz_x2y2]
                ind += 1
                ind2orb[ind] = [i,t, :fxyz]
                ind += 1
                ind2orb[ind] = [i,t, :fx_x2_3y2]
                ind += 1
                ind2orb[ind] = [i,t, :fy_3x2_y2]
            else
                error("orbitals ????")
            end

        end
        if ind != round(nwan/2)
            error("counting error", ind, " " , nwan)
        end
        
    end
            
    nwan = 0
    for (i,t) in enumerate(c.types)
        atom = atoms[t]
        orb2ind[i] = 1+nwan : nwan + Int(round(atom.nwan/2))
        nwan += Int(round(atom.nwan/2))
    end

    #reverse dictionary
#    for ind = 1:convert(Int, round(nwan/2))
#        orb2ind[ind2orb[ind]] = ind
#        println(ind, " ", ind2orb[ind][1], " ", ind2orb[ind][2])
#        
#    end

    return ind2orb, orb2ind, etotal, nval
    
end

end #end module

#using .CrystalMod
#Base.:(==)(x::crystal, y::crystal) = (x.A == y.A && x.coords == y.coords && x.types == y.types && x.nat==y.nat)
