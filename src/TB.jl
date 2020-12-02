###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module TB
"""
Scripts to run tight-binding from wannier90
"""

using LinearAlgebra
using SpecialFunctions
#using EzXML
using XMLDict
using GZip
using Printf
using PyPlot
using FFTW
using JLD
using Base.Threads

#include("Atomdata.jl")
using ..Atomdata:atoms
#include("Commands.jl")

using ..DFToutMod:bandstructure
using ..DFToutMod:dftout
using ..AtomicMod:atom
using ..CrystalMod:crystal
using ..CrystalMod:makecrys
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..BandTools:calc_fermi
using ..BandTools:band_energy
using ..BandTools:gaussian
using ..BandTools:smearing_energy

using ..CrystalMod:get_grid
using ..Ewald:electrostatics_getgamma
using ..CrystalMod:orbital_index

export tb
export tb_crys
export tb_k
export tb_crys_kspace
export read_tb_crys
export read_tb_crys_kspace
export write_tb_crys
export write_tb_crys_kspace
export make_tb_crys
export make_tb_crys_kspace
export make_tb
export make_tb_k
export load_hr_dat
export Hk

export calc_bands
export plot_compare_tb
export plot_bandstr
export plot_compare_dft



symbol_dict = Dict()
#        1   2    3    4    5     6    7     8      9   
list = [:s, :pz, :px, :py, :dz2,:dxz,:dyz,:dx2_y2,:dxy,:fz3,:fxz2,:fyz2,:fz_x2y2,:fxyz,:fx_x2_3y2,:fy_3x2_y2]

for i = 1:16
    symbol_dict[list[i]] = i
end


    
mutable struct tb{T}

    #    H::Array{Complex{Float64},3}
    H::Array{Complex{T},3}    
    ind_arr::Array{Int64,2}
    r_dict::Dict
    nwan::Int
    nr::Int
    nonorth::Bool
    #    S::Array{Complex{Float64},3}
    S::Array{Complex{T},3}    
    scf::Bool
    h1::Array{T,2} #scf term
              
end

Base.show(io::IO, h::tb) = begin
    println(io)
    nwan=h.nwan
    nr=h.nr
    nonorth = h.nonorth
    scf = h.scf
    println(io, "tight binding real space object; nwan = $nwan, nr = $nr, nonorth = $nonorth, scf = $scf" )
    println(io)
    
end   

mutable struct tb_crys{T}

    tb::tb
    crys::crystal
    nelec::Float64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    eden::Array{Float64,1}
    within_fit::Bool
end



Base.show(io::IO, x::tb_crys) = begin
    println(io)
    println(io, x.crys)
    println(io)
    println(io, "nelec: ", x.nelec)
    println(io)
    println(io, x.tb)    
    println(io)
    
end   


mutable struct tb_k{T}


    Hk::Array{Complex{T},3}    
    K::Array{Float64,2}
    kweights::Array{Float64,1}
    k_dict::Dict
    nwan::Int64
    nk::Int64
    nonorth::Bool
    Sk::Array{Complex{T},3}    
    scf::Bool
    h1::Array{T,2} #scf term
    grid::Array{Int64,1}

end


Base.show(io::IO, h::tb_k) = begin
    println(io)
    nwan=h.nwan
    nk=h.nk
    nonorth = h.nonorth
    scf = h.scf
    grid = h.grid
    println(io, "tight binding k-space object; nwan = $nwan, nk = $nk, nonorth = $nonorth, scf = $scf, grid = $grid" )
    println(io)
    
end   



mutable struct tb_crys_kspace{T}

    tb::tb_k
    crys::crystal
    nelec::Float64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    eden::Array{Float64,1}

end



Base.show(io::IO, x::tb_crys_kspace) = begin
    println(io)
    println(io, x.crys)
    println(io)
    println(io, "nelec: ", x.nelec)
    println(io)
    println(io, x.tb)    
    println(io)
    
end   


function get_el(tbc::tb_crys, n1,n2,na,nb,r)

    ind = 1
    try
        ind=tbc.tb.r_dict[r]
        
    catch
        try
            ind=tbc.tb.r_dict[r']
        catch
            println("missing")
            return missing, missing, missing, missing
        end
    end
    
    R = tbc.tb.ind_arr[ind,:]
    dR = tbc.crys.A'*(-tbc.crys.coords[na,:] + tbc.crys.coords[nb,:] + R) 
    dist = sum(dR.^2)^0.5

    
    return tbc.tb.H[n1,n2,ind], tbc.tb.S[n1,n2,ind], dist, dR
    
end

#function get_el(tbc::tb_crys, n1,n2,r)
#
#    return get_el(tbc.tb, n1,n2,r)
#    
#end

function read_tb_crys(filename; directory=missing)
"""
get tbc object from xml file, written by write_tb_crys (see below)
"""
    if ismissing(directory)
        println("read ", filename)
    else
        println("read ", directory*"/"*filename)
    end

    #    try
    f = missing
    if !ismissing(directory)
        filename=directory*"/"*filename
    end
    if !isfile(filename)
        if isfile(filename*".xml")
            filename=filename*".xml"
        elseif isfile(filename*".gz")
            filename=filename*".gz"
        elseif isfile(filename*".xml.gz")
            filename=filename*".xml.gz"
        end
    end
    if !isfile(filename)
        println("warning error read_tb_crys $filename $directory not found")
    else
        println("found $filename")
    end

    try
        f = gzopen(filename, "r")
    catch
        println("error opening $filename")
    end
    
    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]

    ##crys
    nat = parse(Int64, d["crystal"]["nat"])

    A = parse_str_ARR_float(d["crystal"]["A"])
    A = reshape(A, 3,3)'
    
    coords = parse_str_ARR_float(d["crystal"]["coords"])
    
    coords = reshape(coords, 3,nat)'
    
    types = split(d["crystal"]["types"])
    
    crys = makecrys(A,coords,types)
    ##
    ##nelec
    nelec = parse(Float64,d["nelec"])
    dftenergy = parse(Float64,d["dftenergy"])        


    #SCF stuff is optional
    scf = false
    gamma = missing
    eden = missing

    if "scf" in keys(d) scf = parse(Bool,d["scf"])  end
    
    if "gamma" in keys(d)
        gamma = parse_str_ARR_float(d["gamma"])
        s = Int64(round(sqrt(size(gamma)[1])))
        gamma = reshape(gamma, s, s)
        
    end
    if "eden" in keys(d)
        eden  = parse_str_ARR_float(d["eden"])
    end
    

    ##tb
    
    nr = parse(Int64,d["tightbinding"]["nr"])
    nwan = parse(Int64,d["tightbinding"]["nwan"])
    
    nonorth = parse(Bool,d["tightbinding"]["nonorth"])
        
    ind_arr = parse_str_ARR_float(d["tightbinding"]["ind_arr"])
    ind_arr = reshape(ind_arr,3, nr)'
    
    r_dict = Dict()
    for i in 1:nr
        r_dict[copy(ind_arr[i,:])] = i
    end
    #    r_dict = Dict(d["tightbinding"]["r_dict"])

    if "scf" in keys(d["tightbinding"])
        tb_scf = parse(Bool,d["tightbinding"]["scf"])
    else
        tb_scf = false
    end
    if "h1" in keys(d["tightbinding"])
        h1 = parse_str_ARR_float(d["tightbinding"]["h1"])
        h1 = reshape(h1, nwan, nwan)
    else
        h1 = missing
    end
    if tb_scf == false
        h1 =  missing
    end

    
    function readstr(st)

        H = zeros(Complex{Float64}, nwan,nwan,nr)
        S = zeros(Complex{Float64}, nwan,nwan,nr)        
        
        lines = split(st, "\n")
        for line in lines
            sp = split(line)
            if length(sp) == 7
                m = parse(Int64,sp[2])
                n = parse(Int64,sp[3])
                r = parse(Int64,sp[1])
                H[m,n,r] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                S[m,n,r] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
            end
        end

        return H, S
    end
        
    H,S = readstr(d["tightbinding"]["H"])

#    println("read ")
#    println("h1")
#    println(h1)

    if nonorth
        tb = make_tb(H, ind_arr, r_dict, S, h1=h1)
    else
        tb = make_tb(H, ind_arr, r_dict, h1=h1)
    end    

    tbc = make_tb_crys(tb, crys, nelec, dftenergy, scf=scf, eden=eden, gamma=gamma)

    return tbc
        
#    catch
#        println("ERROR: Failed to read ", filename)
#        return -1
#    end
    
end



function read_tb_crys_kspace(filename; directory=missing)
"""
get tbc object from xml file, written by write_tb_crys (see below)
"""
    if ismissing(directory)
        println("read ", filename)
    else
        println("read ", directory*"/"*filename)
    end

    f = missing
    if !ismissing(directory)
        filename=directory*"/"*filename
    end
    if !isfile(filename)
        if isfile(filename*".xml")
            filename=filename*".xml"
        elseif isfile(filename*".gz")
            filename=filename*".gz"
        elseif isfile(filename*".xml.gz")
            filename=filename*".xml.gz"
        end
    end
    if !isfile(filename)
        println("warning error read_tb_crys $filename $directory not found")
    else
        println("read_tb_crys_kspace found $filename")
    end

    try
        f = gzopen(filename, "r")
    catch
        println("error opening $filename")
    end


    #    try
#    f = missing
#    try
#        if ismissing(directory)
#            f = gzopen(filename, "r")
#        else
#            f = gzopen(directory*"/"*filename, "r")
#        end
#    catch
#        filename=filename*".xml"
#            println("trying ", filename)
#        if ismissing(directory)
#            f = gzopen(filename, "r")
#        else
#            f = gzopen(directory*"/"*filename, "r")
#       end
#    end
    
    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]

    ##crys
    nat = parse(Int64, d["crystal"]["nat"])

    A = parse_str_ARR_float(d["crystal"]["A"])
    A = reshape(A, 3,3)'
    
    coords = parse_str_ARR_float(d["crystal"]["coords"])
    
    coords = reshape(coords, 3,nat)'
    
    types = split(d["crystal"]["types"])
    
    crys = makecrys(A,coords,types)
    ##
    ##nelec
    nelec = parse(Float64,d["nelec"])
    dftenergy = parse(Float64,d["dftenergy"])        


    #SCF stuff is optional
    scf = false
    gamma = missing
    eden = missing

    if "scf" in keys(d) scf = parse(Bool,d["scf"])  end
    
    if "gamma" in keys(d)
        gamma = parse_str_ARR_float(d["gamma"])
        s = Int64(round(sqrt(size(gamma)[1])))
        gamma = reshape(gamma, s, s)
        
    end
    if "eden" in keys(d)
        eden  = parse_str_ARR_float(d["eden"])
    end
    

    ##tb
    
    nk = parse(Int64,d["tightbinding"]["nk"])
    nwan = parse(Int64,d["tightbinding"]["nwan"])
    
    nonorth = parse(Bool,d["tightbinding"]["nonorth"])
        
    kind_arr = parse_str_ARR_float(d["tightbinding"]["kind_arr"])
    kind_arr = reshape(kind_arr,3, nk)'

    kweights = parse_str_ARR_float(d["tightbinding"]["kweights"])

    if "grid" in keys(d["tightbinding"])
#        grid = parse(Int64,d["tightbinding"]["grid"])
#        grid = parse_str_ARR_float(d["tightbinding"]["kweights"])
        grid = Int64.(parse_str_ARR_float(d["tightbinding"]["grid"]))

    else
        grid = [0,0,0]
    end
    
    k_dict = Dict()
    for i in 1:nk
        k_dict[copy(kind_arr[i,:])] = i
    end
    #    r_dict = Dict(d["tightbinding"]["r_dict"])

    if "scf" in keys(d["tightbinding"])
        tb_scf = parse(Bool,d["tightbinding"]["scf"])
    else
        tb_scf = false
    end
    if "h1" in keys(d["tightbinding"])
        h1 = parse_str_ARR_float(d["tightbinding"]["h1"])
        h1 = reshape(h1, nwan, nwan)
    else
        h1 = missing
    end
    if tb_scf == false
        h1 =  missing
    end

    
    function readstr(st)

        H = zeros(Complex{Float64}, nwan,nwan,nk)
        S = zeros(Complex{Float64}, nwan,nwan,nk)        
        
        lines = split(st, "\n")
        for line in lines
            sp = split(line)
            if length(sp) == 7
                m = parse(Int64,sp[2])
                n = parse(Int64,sp[3])
                r = parse(Int64,sp[1])
                H[m,n,r] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                S[m,n,r] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
            end
        end

        return H, S
    end
        
    Hk,Sk = readstr(d["tightbinding"]["Hk"])

#    println("read ")
#    println("h1")
#    println(h1)

    tb = make_tb_k(Hk, kind_arr, kweights, Sk, h1=h1, grid=grid, nonorth=nonorth)

    tbck = make_tb_crys_kspace(tb, crys, nelec, dftenergy, scf=scf, eden=eden, gamma=gamma)

    return tbck
        
#    catch
#        println("ERROR: Failed to read ", filename)
#        return -1
#    end
    
end


function write_tb_crys(filename, tbc::tb_crys)
    """
    write xml tb_crys object
    """

    doc = XMLDocument()
    root = ElementNode("root")
    setroot!(doc, root)

    crystal = ElementNode("crystal")
    link!(root, crystal)
    
    addelement!(crystal, "A", arr2str(tbc.crys.A))
    addelement!(crystal, "nat", string(tbc.crys.nat))
    addelement!(crystal, "coords", arr2str(tbc.crys.coords))
    addelement!(crystal, "types", str_w_spaces(tbc.crys.types))

    addelement!(root, "nelec", string(tbc.nelec))
    addelement!(root, "dftenergy", string(tbc.dftenergy))
    addelement!(root, "scf", string(tbc.scf))
    addelement!(root, "gamma", arr2str(tbc.gamma))
    addelement!(root, "eden", arr2str(tbc.eden))

    tightbinding = ElementNode("tightbinding")
    link!(root, tightbinding)
    
    addelement!(tightbinding, "nr", string(tbc.tb.nr))
    addelement!(tightbinding, "nwan", string(tbc.tb.nwan))
    addelement!(tightbinding, "nonorth", string(tbc.tb.nonorth))
    addelement!(tightbinding, "ind_arr",arr2str(tbc.tb.ind_arr))
    addelement!(tightbinding, "r_dict",string(tbc.tb.r_dict))

    addelement!(tightbinding, "scf",string(tbc.tb.scf))
    addelement!(tightbinding, "h1",arr2str(tbc.tb.h1))
    
    function makestr(H,S, nonorth, ind_arr)

        nw = size(H)[1]
        nr = size(H)[3]
#        st = ""
        io_tmp = IOBuffer() #writing to iostream is much faster than making a giant string 
                            #strings in julia are immutable, so the entire thing had to be copied over and over

        for r = 1:nr
            for m = 1:nw
                for n = 1:nw            
                    if nonorth
#                        t= @sprintf("% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),real(S[m,n,r]),imag(S[m,n,r]))
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),real(S[m,n,r]),imag(S[m,n,r]))
                    else
                        if (m == n) && ind_arr[r,1] == 0 && ind_arr[r,2] == 0 && ind_arr[r,3] == 0 
                            @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),1.0,0.0)
                        else
                            @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),0.0,0.0)
                        end
                    end
#                    st = st*t

                end
            end
        end
        st= String(take!(io_tmp))
        close(io_tmp)

        return st

    end

    st = makestr(tbc.tb.H, tbc.tb.S, tbc.tb.nonorth, tbc.tb.ind_arr)
    
    addelement!(tightbinding, "H", st)
    
#    println(doc)
#
#    println()

    #prettyprint
    io=open(filename, "w")
    prettyprint(io, doc);
    close(io)

    #  write(io, doc);
#    write(filename, doc);    

    return Nothing
    
#    return doc
    
end

function write_tb_crys_kspace(filename, tbc::tb_crys_kspace)
    """
    write xml tb_crys kspace object
    """

    doc = XMLDocument()
    root = ElementNode("root")
    setroot!(doc, root)

    crystal = ElementNode("crystal")
    link!(root, crystal)
    
    addelement!(crystal, "A", arr2str(tbc.crys.A))
    addelement!(crystal, "nat", string(tbc.crys.nat))
    addelement!(crystal, "coords", arr2str(tbc.crys.coords))
    addelement!(crystal, "types", str_w_spaces(tbc.crys.types))

    addelement!(root, "nelec", string(tbc.nelec))
    addelement!(root, "dftenergy", string(tbc.dftenergy))
    addelement!(root, "scf", string(tbc.scf))
    addelement!(root, "gamma", arr2str(tbc.gamma))
    addelement!(root, "eden", arr2str(tbc.eden))

    tightbinding = ElementNode("tightbinding")
    link!(root, tightbinding)
    
    addelement!(tightbinding, "nk", string(tbc.tb.nk))
    addelement!(tightbinding, "nwan", string(tbc.tb.nwan))
    addelement!(tightbinding, "nonorth", string(tbc.tb.nonorth))
    addelement!(tightbinding, "kind_arr",arr2str(tbc.tb.K))
    addelement!(tightbinding, "kweights",arr2str(tbc.tb.kweights))
    addelement!(tightbinding, "k_dict",string(tbc.tb.k_dict))

    addelement!(tightbinding, "scf",string(tbc.tb.scf))
    addelement!(tightbinding, "h1",arr2str(tbc.tb.h1))

    addelement!(tightbinding, "grid",arr2str(tbc.tb.grid))
    
    function makestr(H,S, nonorth, ind_arr)

        nw = size(H)[1]
        nr = size(H)[3]
#        st = ""
        io_tmp = IOBuffer() #writing to iostream is much faster than making a giant string 
                            #strings in julia are immutable, so the entire thing had to be copied over and over

        for r = 1:nr
            for m = 1:nw
                for n = 1:nw            
                    if nonorth
#                        t= @sprintf("% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),real(S[m,n,r]),imag(S[m,n,r]))
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),real(S[m,n,r]),imag(S[m,n,r]))
                    else
                        if (m == n) && ind_arr[r,1] == 0 && ind_arr[r,2] == 0 && ind_arr[r,3] == 0 
                            @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),1.0,0.0)
                        else
                            @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r]), imag(H[m,n,r]),0.0,0.0)
                        end
                    end
#                    st = st*t

                end
            end
        end
        st= String(take!(io_tmp))
        close(io_tmp)

        return st

    end

    st = makestr(tbc.tb.Hk, tbc.tb.Sk, tbc.tb.nonorth, tbc.tb.K)
    
    addelement!(tightbinding, "Hk", st)
    
#    println(doc)
#
#    println()

    #prettyprint
    io=open(filename, "w")
    prettyprint(io, doc);
    close(io)

    #  write(io, doc);
#    write(filename, doc);    

    return Nothing
    
#    return doc
    
end
    
function make_tb_crys(ham::tb,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, within_fit=true, screening=1.0)

    T = typeof(crys.coords[1,1])

    if ismissing(eden)
        if scf == false
            eden = zeros(ham.nwan)
        else
            eden = get_neutral_eden(crys, ham.nwan)
        end
    end
    
    if ismissing(gamma) 
#        println("ewald time")
        gamma = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end

#    println("type scf " , typeof(scf))
#    println("type gamma " , typeof(gamma))
#    println("type eden " , typeof(eden))
    
#    return tb_crys{T}(ham,crys,nelec, dftenergy, scf, gamma, eden)
    return tb_crys{T}(ham,crys,nelec, dftenergy, scf, gamma, eden, within_fit)
end

function make_tb_crys_kspace(hamk::tb_k,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, screening=1.0)

    T = typeof(crys.coords[1,1])

    if ismissing(eden)
        if scf == false
            eden = zeros(hamk.nwan)
        else
            eden = get_neutral_eden(crys, hamk.nwan)
        end
    end
    
    if ismissing(gamma) 
        gamma = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end

    return tb_crys_kspace{T}(hamk,crys,nelec, dftenergy, scf, gamma, eden)
end



function make_tb(H, ind_arr, r_dict::Dict; h1=missing)
    nw=size(H,2)
    if nw != size(H)[1]
        s=size(H)
        error("error in make_tb $s")
    end
    nr=size(H,3)

    
    if size(ind_arr) != (nr,3)
        error("make_tb ind_arr ", size(ind_arr))
    end

    T=typeof(real(H[1,1,1]))
    
    S =zeros(Complex{T}, size(H))

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end
    
    return tb{T}(H, ind_arr, r_dict,nw, nr, false, S, scf, h1)
end

function make_tb(H, ind_arr, r_dict::Dict, S; h1=missing)
    nw=size(H,1)
    if nw != size(H,1) 
        exit("error in make_tb ", size(H))
    end
    nr=size(H,3)

    if size(ind_arr) != (nr,3) 
        error("make_tb ind_arr size ", size(ind_arr), size(H))
    end

    if size(S) != (nw,nw,nr)
        error("make_tb size S doesn't match size H", size(S), size(H))
    end
    
    T=typeof(real(H[1,1,1]))

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end
    
    return tb{T}(H, ind_arr, r_dict,nw, nr, true, S, scf, h1)
end





function make_tb_k(Hk, K, kweights, Sk; h1=missing, grid=[0,0,0], nonorth=true)
    nw=size(Hk,1)
    nk=size(Hk,3)

#    println("make_tb_k $nw $nk xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx") 
#    for k in 1:nk
#        println("K $k ", K[k,:])
#    end

    if size(K) != (nk,3) 
        error("make_tb_k kind_arr size ", size(K), size(Hk))
    end

    if size(Sk) != (nw,nw,nk)
        error("make_tb_k size S doesn't match size H", size(Sk), size(Hk))
    end
    
    T=typeof(real(Hk[1,1,1]))

    K2 = zeros(Float64,size(K))
    k_dict = Dict()  
    for i in 1:nk
        K2[i,:] = Int64.(round.(K[i,:] * 100000000))/100000000
        k_dict[K2[i,:]] = i
    end

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end
#=
    println("make_tb_k")
    println(typeof(Hk))
    println(typeof(K2))
    println(typeof(kweights))
    println(typeof(k_dict))
    println(typeof(nw))
    println(typeof(nk))
    println(typeof(true))
    println(typeof(Sk))
    println(typeof(scf))
    println(typeof(h1))
    println(typeof(grid))
=#
    return tb_k{T}(Hk, K2, kweights, k_dict,nw, nk, nonorth, Sk, scf, h1, grid)
end


function make_rdict(ind_arr)
    r_dict = Dict()
    for c = 1:size(ind_arr)[1]
        r_dict[ind_arr[c,:]] = c
    end
    return r_dict
end

        


function make_tb(H, ind_arr, S; h1=missing)

    r_dict = make_rdict(ind_arr)
    
    nw=size(H,1)
    if nw != size(H,1) 
        exit("error in make_tb ", size(H))
    end
    nr=size(H,3)

    if size(ind_arr) != (nr,3) 
        error("make_tb ind_arr size ", size(ind_arr), size(H))
    end

    if size(S) != (nw,nw,nr)
        error("make_tb size S doesn't match size H", size(S), size(H))
    end
    
    vartype=typeof(H[1,1,1])

    if ismissing(h1)
        h1 =zeros(Float64, nw, nw)
        scf = false
    else
        scf = true
    end

    
    return tb{vartype}(H, ind_arr, r_dict,nw, nr, true, S, scf, h1)
end



function load_hr_dat(filename, directory="")
"""
Load hr.dat file from wannier90
"""

    convert_ev_ryd = 1.0/13.605693122
    if directory != ""
        f=gzopen("$directory/$filename", "r")
    else
        f=gzopen("$filename", "r")
    end
    
    hr_file = readlines(f)
    close(f)

    println(hr_file[1]) #header
    nwan = parse(Int,hr_file[2])          
    nr = parse(Int,hr_file[3])          

    println("nwan:", nwan)
    println("nr:", nr)    


    #parse weights
    if nr % 15 == 0
        wlines = convert(Int, floor(nr / 15))
    else
        wlines = convert(Int, floor(nr / 15)) + 1
    end
    println("wlines ", wlines)
    weights = zeros(nr)
    for i = 4:(3+wlines)
        line = hr_file[i]

        ind_start = (i-4)*15 + 1
        ind_end = (i-3)*15

#        println(line)
#        println(i, " ", ind_start," ", ind_end)

        if ind_end > nr
            ind_end = nr
        end
        weights[ind_start:ind_end] = map(x->parse(Float64,x),split(line))

    end

#    println(weights)

    #now parse main Ham

    H = zeros(Complex{Float64},nwan,nwan,nr)

    r_dict = Dict()
    ind_arr = zeros(nr,3)
    
    nind = 0

    nonorth=false
    sp = split(hr_file[4+wlines])
    if length(sp) == 9
        nonorth = true
        S = zeros(Complex{Float64}, nwan,nwan, nr)
        
    end


    
    for i = (4+wlines):length(hr_file)
        sp = split(hr_file[i])

        rind = [parse(Int, sp[1]), parse(Int, sp[2]),parse(Int, sp[3])]
        
        if rind in keys(r_dict)
            ind = r_dict[rind]
        else
            nind += 1
            r_dict[rind] = nind
            ind = nind
            ind_arr[ind,:] = rind
        end

        
        
        nw1 = parse(Int64,sp[4])
        nw2 = parse(Int64,sp[5])
        
        h = (parse(Float64, sp[6]) + parse(Float64, sp[7])*im) / weights[ind] * convert_ev_ryd
        H[nw1,nw2,ind] = h

        if nonorth
            s = (parse(Float64, sp[8]) + parse(Float64, sp[9])*im) / weights[ind]
            S[nw1,nw2,ind] = s

        end
        
    end

#    println(H[1,1,1])
#    println(H[2,2,2])

    if nonorth
        return make_tb(H, ind_arr, r_dict, S)
    else
        return make_tb(H, ind_arr, r_dict)
    end
end

function Hk(h::tb_crys_kspace, kpoint; scf=missing)

    return Hk(h.tb, kpoint, scf=scf)

end

function Hk(h::tb_k, kpoint; scf=missing) 

    if ismissing(scf)
        scf = h.scf
    end
    
    kpoint = vec(kpoint)

    hk = zeros(Complex{Float64}, h.nwan, h.nwan)
    sk = zeros(Complex{Float64}, h.nwan, h.nwan)
    hk0 = zeros(Complex{Float64}, h.nwan, h.nwan)
    
    if kpoint in keys(h.k_dict)
        ind = h.k_dict[kpoint]
    else
        ind = 0
        for i = 1:h.nk
#            println(kpoint')
#            println(size(kpoint'))
#            println(h.K[i,:])
#            println(size(h.K[i,:]))
            if sum(abs.(kpoint - h.K[i,:])) < 5e-4
                ind = i
                break
            end
        end
        if ind == 0
            println("error, k not found!!! , $kpoint tb_k")
            println(tb)
            for i = 1:h.nk
                println(h.K[i,:])
            end

            return
        end
    end

    hk0[:,:] = h.Hk[:,:,ind]
    hk0 = 0.5*(hk0 + hk0')
    
    sk[:,:] = h.Sk[:,:,ind]
    sk = 0.5*(sk + sk')
    
    if scf
#        println("add SCF")
        hk = hk0 + 0.5 * sk .* (h.h1 + h.h1')
    else
        hk[:,:] = hk0[:,:]
    end
    hk = 0.5*(hk[:,:] + hk[:,:]')
    
    nw=size(hk)[1]
    vects = zeros(nw,nw)
    vals = zeros(nw)
    vals0 = zeros(nw)

    try
        if h.nonorth
            sk = 0.5*(sk[:,:] + sk[:,:]')
            F=eigen(hk[:,:], sk[:,:])
        else
            F=eigen(Hermitian(hk)) #orthogonal
        end

        vects = F.vectors
        vals = real(F.values)
        vals0 = real.(diag(vects'*hk0*vects))
        
    
    catch
        println("warning eigen failed, ", kpoint)
        vects = collect(I(nw))
        vals = 1000.0 * ones(nw)
        vals0 = 1000.0 * ones(nw)
        
    end
    return vects, vals, hk, sk, vals0
    
end


function Hk(hk,sk, h::tb, kpoint)

    kpoint = vec(kpoint)

#    println("Hk kpoint : $kpoint")
    
#    println("Hk nonorth: " , h.nonorth)
    
#    println("hk", typeof(hk), size(hk))

#    println("zero")

    hk0 = zeros(Complex{Float64}, size(hk))
    
    fill!(hk, zero(Float64))

    if h.nonorth  
        fill!(sk, zero(Float64))
    end
        
##    Hk = zeros(Complex{Float64}, h.nwan, h.nwan)

    twopi_i = -1.0im*2.0*pi

#    println("repeat")
    kmat = repeat(kpoint', h.nr,1)
    
#    println("exp")
    exp_ikr = exp.(twopi_i * sum(kmat .* h.ind_arr,dims=2))

    for m in 1:h.nwan
        for n in 1:h.nwan        
            hk0[m,n] = h.H[m,n,:]'*exp_ikr[:]
            if h.nonorth
                sk[m,n] = h.S[m,n,:]'*exp_ikr[:]
            end
        end
    end
    hk0 = 0.5*(hk0 + hk0')
    
    if h.scf
#        println("add SCF")
        hk = hk0 + sk .* h.h1
    else
        hk[:,:] = hk0[:,:]
    end
    hk = 0.5*(hk[:,:] + hk[:,:]')
    
#    ex = 0.0 + im*0.0
#    for n = 1:h.nr
 #       ex = exp(-twopi_i*(transpose(h.ind_arr[n,:])*kpoint))
 #       hk .+= ex .* h.H[:,:,n]
 #       if h.nonorth
 #           sk .+= ex .* h.S[:,:,n]
 #       end
 #   end
    #end

    
    

    nw=size(hk)[1]
    vects = zeros(nw,nw)
    vals = zeros(nw)
    vals0 = zeros(nw)

    try
        if h.nonorth
            sk = 0.5*(sk[:,:] + sk[:,:]')
            F=eigen(hk[:,:], sk[:,:])
        else
            #        println("orth")
            #        println(typeof(hk))
            hk = 0.5*(hk[:,:] + hk[:,:]')            
            F=eigen(hk[:,:]) #orthogonal
        end

        vects = F.vectors
        vals = real(F.values)
        vals0 = real.(diag(vects'*hk0*vects))


        
    catch
        println("warning eigen failed, ", kpoint)
        vects = collect(I(nw))
        vals = 1000.0 * ones(nw)
        vals0 = 1000.0 * ones(nw)
        
    end
    
    #    F=eigen(hk, sk)

#    println("hk")
#    println(hk)

#    println("hk vals")
#    println(vals)

    return vects, sort!(vals), hk, sk, vals0
    
end

function compare_bands(h::tb, bs::bandstructure; energy_lim=missing, doplot=true)

    if ismissing(energy_lim)
        println("comparing below Ef ", bs.efermi)
        energy_lim = bs.efermi
    end
    
    tb_vals = calc_bands(h, bs.kpts)

    max_diff = 0.0
    badk = -1
    badval = -999.0

    vals1 = zeros(bs.nks, h.nwan)
    vals2 = zeros(bs.nks, h.nwan)

    m1 = minimum(bs.eigs)
    m2 = minimum(tb_vals)
    
    for k = 1:bs.nks
        for w = 1:h.nwan
            
            val = tb_vals[k,w]

            vals1[k,w] = bs.eigs[k,w] - m1
            vals2[k,w] = val - m2
            
            if val < energy_lim
                m = minimum(abs.( (bs.eigs[k,:] .- m1) .- (val .- m2)))
                if m > max_diff
                    max_diff = m
                    badk = k
                    badval = w
                end
            end
        end
    end
    if doplot
        plot(vals1, "g.", MarkerSize=8)
        plot(vals2, "y.", MarkerSize=4)
    end
    println("max abs diff $max_diff at kpoint $badk and value $badval")
    if max_diff > 1e-3
        println("comparison")
        println("tb  ", tb_vals[badk,:])
        println("dft ", bs.eigs[badk,:])    
    end
    return max_diff
end


function calc_bands(tbc::tb_crys, kpoints::Array{Float64,2})
#    if tbc.scf == true
#        h1, dq = get_h1(tbc, tbc.eden)
#    else
#        h1 == missing
#    end
    return calc_bands(tbc.tb, kpoints)
end


function calc_bands(h, kpoints::Array{Float64,2})
"""
calculate band structure at group of points
"""
#    println("calc_bands")
#    println("h.scf ", h.scf)
#    println("temp arrays")

#    T=Float64
#    
#    hktemp= zeros(Complex{T}, h.nwan, h.nwan)
#    sktemp= zeros(Complex{T}, h.nwan, h.nwan)

    nk = size(kpoints,1)

#    println("vals time")
    
    Vals = zeros(Float64, nk,h.nwan)
    
    for i = 1:nk

        vect, vals, hk, sk, vals0 = Hk(h, kpoints[i,:])
        Vals[i,:] = vals
#        if sum(abs.(kpoints[i,:])) == 0
#            println("k0")
#            println(vals)
#        end

    end

    return Vals

end    

function plot_compare_tb(h1::tb_crys, h2::tb_crys; h3=missing, kpath=[0.5 0 0 ; 0 0 0; 0.5 0.0 0.5], names = missing, npts=30, efermi = missing, yrange=missing, plot_hk=false,  align=false)
    if ismissing(h3)
        plot_compare_tb(h1.tb, h2.tb, h3=missing, kpath=kpath, names = names, npts=npts, efermi = efermi, yrange=yrange, plot_hk=plot_hk, align=align)
    else
        plot_compare_tb(h1.tb, h2.tb, h3=h3.tb, kpath=kpath, names = names, npts=npts, efermi = efermi, yrange=yrange, plot_hk=plot_hk, align=align)
    end
end


function plot_compare_tb(h1::tb, h2::tb; h3=missing, kpath=[0.5 0 0 ; 0 0 0; 0.5 0.0 0.5], names = missing, npts=30, efermi = missing, yrange=missing, plot_hk=false, align=false)
#    println("kpath ", kpath)
    plot_bandstr(h1, kpath=kpath, names = names, npts=npts, efermi = efermi, colors="-g.", MarkerSize=10, yrange=yrange, plot_hk=plot_hk, align=align)
    plot_bandstr(h2, kpath=kpath, names = names, npts=npts, efermi = efermi, colors="--y.", MarkerSize=6, yrange=yrange, plot_hk=plot_hk, align=align)    
    if !ismissing(h3)
        plot_bandstr(h3, kpath=kpath, names = names, npts=npts, efermi = efermi, colors="--m.", MarkerSize=4, yrange=yrange, plot_hk=plot_hk, align=align)
    end

end

function summarize_orb(orb::Symbol)
    if orb == :s
        return :s
    elseif orb == :p || orb == :px || orb == :py || orb == :pz
        return :p
    elseif orb == :d || orb == :dz2 || orb == :dxz || orb == :dyz || orb == :dx2_y2 || orb == :dxy
        return :d
    else
        println("not coded yet $orb summarize_orb")
    end
end

    
function plot_bandstr(h::tb_crys; kpath=[0.5 0 0 ; 0 0 0; 0.5 0.0 0.5], names = missing, npts=30, efermi = missing, colors="-b.", MarkerSize=8, yrange=missing, plot_hk=false, align = false, proj_types = missing, proj_orbs = missing, proj_nums=missing)


    if !ismissing(proj_types) || !ismissing(proj_orbs) || !ismissing(proj_nums)
        
        if ismissing(proj_types)
            proj_types = h.crys.types
        end
        if ismissing(proj_orbs)
            proj_orbs = [:s, :p, :d, :f]
        end
        if ismissing(proj_nums)
            proj_nums = collect(1:h.crys.nat)
        end
        if typeof(proj_types) == String
            proj_types = [proj_types]
        end
        if typeof(proj_orbs) == Symbol
            proj_orbs = [proj_orbs]
        end
        if typeof(proj_nums) == Int64
            proj_nums = [proj_nums]
        end

        ind2orb, orb2ind, etotal, nval = orbital_index(h.crys)
        
        proj_inds = Int64[]
        for n = 1:h.tb.nwan
            at,t, orb = ind2orb[n]
            sorb = summarize_orb(orb)
            if (t in proj_types) && (orb in proj_orbs || sorb in proj_orbs) && (at in proj_nums)
                push!(proj_inds, n)
            end
        end
        println("proj_inds ", proj_inds)
    else
        proj_inds = missing
    end

    plot_bandstr(h.tb; kpath=kpath, names = names, npts=npts, efermi = efermi, colors=colors, MarkerSize=MarkerSize, yrange=yrange, plot_hk=plot_hk, align=align, proj_inds=proj_inds)
    
end


function get_kpath(kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5], names = missing, npts=30)

    NK = size(kpath)[1]
    K = zeros(npts *(NK-1)+1 , 3)

    nk = 0

    locs = []
    
    for kp = 1:(NK-1)

        push!(locs,nk)
        
        for x = 0:1/npts:(1-1e-5)
            nk += 1
            kt = kpath[kp,:] * (1-x) .+ kpath[kp+1,:] * x
#            println(kt)
            K[nk,:] = kt
        end

    end

#final point
    kp = NK-1
    x = 1.0
    nk += 1
    kt = kpath[kp,:] * (1-x) .+ kpath[kp+1,:] * x
#    println(kt)
    K[nk,:] = kt
    
    push!(locs,nk-1)

#    println("locs: ", locs)
    
    return K, locs
    
end



function plot_bandstr(h::tb; kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5; 0 0.5 0.5; 0 0 0 ;0 0 0.5], names = missing, npts=30, efermi = missing, colors="-b.", MarkerSize=8, yrange=missing, plot_hk=false, align=false, proj_inds=missing)
#function plot_bandstr( kpath; names = missing, npts=30, efermi = missing)

    NK = size(kpath)[1]

    K, locs = get_kpath(kpath, names , npts)
    nk = size(K)[1]
    
    proj = zeros(nk, h.nwan)

    if plot_hk != false

        if plot_hk == :Seig ||  plot_hk == :Heig
            vals = zeros(nk, h.nwan)
        else
            vals = zeros(nk, h.nwan*h.nwan)
        end
            
        for i = 1:nk
            vect, vals_t, hk, sk, vals0 = Hk(h, K[i,:])
            if plot_hk == :Hreal
                vals[i,:] = real(hk[:])
            elseif plot_hk == :Himag
                vals[i,:] = imag(hk[:])
            elseif plot_hk == :Sreal
                vals[i,:] = real(sk[:])
            elseif plot_hk == :Simag
                vals[i,:] = imag(sk[:])
            elseif plot_hk == :Seig
                F=eigen(sk)
                vals[i,:] = real(F.values)
            elseif plot_hk == :Heig
                F=eigen(hk)
                vals[i,:] = real(F.values)
            else
                vals = calc_bands(h, K)
            end
        end

    elseif !ismissing(proj_inds )
        vals = zeros(nk, h.nwan)
        temp = zeros(Complex{Float64}, h.nwan, h.nwan)
        println("proj inds $proj_inds")
        for i = 1:nk
            vect, vals_t, hk, sk, vals0 = Hk(h, K[i,:])
            for p in proj_inds
                for a = 1:h.nwan
                    for b = 1:h.nwan
                        temp[a,b] = vect[a,p]'*sk[a, b]* vect[b,p]
                    end
                end
                temp = temp + conj(temp)
                proj[i,:] += 0.5*sum(real(temp), dims=1)[:]
            end
            #            println(vect[proj_inds,:] .* conj(vect[proj_inds,:]))
#            println()
#            proj[i,:] = sum(vect[proj_inds,:]*vect[proj_inds, :]', dims=1)
            #            proj[i,:] = real(diag(vect[proj_inds,:]'*sk[proj_inds, proj_inds]*vect[proj_inds,:]))
            vals[i,:] = vals_t
        end
                
        
    else
        vals = calc_bands(h, K)
    end

    if align
        vmin = minimum(vals)
        vals = vals .- vmin
    end

    if ismissing(proj_inds)
        plot(vals, colors, MarkerSize=MarkerSize)
    else
        X=repeat(0:(nk-1), 1,h.nwan)
        plot(X, vals, Color="grey", zorder=1)
        scatter(X, vals, 20.0*ones(size(vals)), proj, zorder=2, cmap="viridis", edgecolors="face", alpha=0.6)
#        @save "t.jld" vals proj
    end


    if (ismissing(names))
        names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
        names = [names;names;names;names;names;names;names; names; names; names; names;names;names;names ]
        names = names[1:NK]
    end

    if !(ismissing(yrange))
        if size(yrange) == (1,2)
            yrange=yrange'
        end
        ylim(yrange[1], yrange[2])

    end

    if ismissing(yrange)
        r = maximum(vals) - minimum(vals) 
        yrange = [minimum(vals)-0.05*r,  maximum(vals)+0.05*r]
    end
    
    ylim(yrange[1], yrange[2])
    xlim(1, NK)

    if length(locs) < 10
        for l in locs
            plot([l, l], yrange, "--k")
        end
    end

    println(locs)
    println(names)    
    xticks(locs, names)
    
end


#function calc_bands(h::tb, kpoints::Array{Float64,2})
#"""
#calculate band structure at group of points
#"""
#
##    println("temp arrays")
#    hktemp= zeros(Complex{Float64}, h.nwan, h.nwan)
#    sktemp= zeros(Complex{Float64}, h.nwan, h.nwan)
#
#    nk = size(kpoints,1)
#
##    println("vals time")
#    
#    Vals = zeros(Float64, nk,h.nwan)
#    
#    for i = 1:nk
#
#        vect, vals, hk, sk = Hk(hktemp, sktemp, h, kpoints[i,:])
#        Vals[i,:] = vals
#
#    end
#
#    return Vals
#
#end    
#


function Hk(h::tb, kpoint)

    T=typeof(real(h.H[1,1,1]))
    #if we don't have a temp array
    hktemp= zeros(Complex{T}, h.nwan, h.nwan)
    sktemp= zeros(Complex{T}, h.nwan, h.nwan)    
    return Hk(hktemp, sktemp, h, kpoint)
end

function Hk(tbc::tb_crys, kpoint )

    
    return Hk(tbc.tb, kpoint)
    
end


function calc_energy_fft(tbc::tb_crys; grid=missing, smearing=0.01, returnef=false)

    etypes = types_energy(tbc.crys)


    hk3, sk3 = myfft_R_to_K(tbc, grid)

    ret =  calc_energy_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, returnef=returnef)
    if returnef
        etot = etypes + ret[1]
        ef = ret[2]

        return etot, ef
    end
    return ret + etypes

end

function calc_energy_fft_band(hk3, sk3, nelec; smearing=0.01, returnef=false, h1 = missing)
#h1 is the scf contribution

    grid = size(hk3)[3:5]
    nk = prod(grid)
    nwan = size(hk3)[1]
    VALS = zeros(Float64, nk, nwan)
    c=0

    thetype=typeof(real(sk3[1,1,1,1,1]))
    sk = zeros(Complex{thetype}, nwan, nwan)
    hk = zeros(Complex{thetype}, nwan, nwan)

    for k1 = 1:grid[1]
        for k2 = 1:grid[2]
            for k3 = 1:grid[3]
                c += 1
                try
                    
                    sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                    hk[:,:] = 0.5*(hk3[:,:,k1,k2,k3] + hk3[:,:,k1,k2,k3]')

                    if !ismissing(h1)
                        hk += sk .* h1
                    end

                    vals, vects = eigen(hk, sk)
                    VALS[c,:] = vals
                catch
                    println("error calc_energy_fft $k1 $k2 $k3")
                    sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                    valsS, vectsS = eigen(sk)
                    println(valsS)
                end
            end
        end
    end

    band_en, efermi = band_energy(VALS, ones(nk), nelec, smearing, returnef=true)
    energy_smear = smearing_energy(VALS, ones(nk), efermi, smearing)

    if returnef
        return  band_en + energy_smear, efermi
    else
        return  band_en + energy_smear
    end

end 


function calc_energy_charge_fft_band(hk3, sk3, nelec; smearing=0.01, h1 = missing)

#    println("memory")
    if true

        grid = size(hk3)[3:5]
        #    print("calc_energy_charge_fft_band grid $grid")
        nk = prod(grid)
        nwan = size(hk3)[1]
        VALS = zeros(Float64, nk, nwan)
        VALS0 = zeros(Float64, nk, nwan)
        c=0

        thetype=typeof(real(sk3[1,1,1,1,1]))
        sk = zeros(Complex{thetype}, nwan, nwan)
        hk = zeros(Complex{thetype}, nwan, nwan)
        hk0 = zeros(Complex{thetype}, nwan, nwan)

        VECTS = zeros(Complex{thetype}, nk, nwan, nwan)
        cVECTS = zeros(Complex{thetype}, nk, nwan, nwan)
        SK = zeros(Complex{thetype}, nk, nwan, nwan)

        error_flag = false

    end

    #println("grid")
    
    for k1 = 1:grid[1]
        for k2 = 1:grid[2]
            for k3 = 1:grid[3]
                c += 1
                try
                    hk0[:,:] = 0.5*(hk3[:,:,k1,k2,k3] + hk3[:,:,k1,k2,k3]')
                    sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')

                    if !ismissing(h1)
                        hk = hk0 + sk .* h1
                    else
                        hk = hk0
                    end
                    

                    vals, vects = eigen(hk, sk)
                    VALS[c,:] = vals
                    VALS0[c,:] = real.(diag(vects'*hk0*vects))

                    VECTS[c,:,:] = vects
                    SK[c,:,:] = sk
                catch
                    if error_flag == false
                        println("error calc_energy_fft $k1 $k2 $k3")
                        sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                        valsS, vectsS = eigen(sk)
                        println(valsS)
                    end
                    error_flag=true
                end
            end
        end
    end
        

    maxSK = maximum(abs.(SK), dims=1)[1,:,:]


    energy, efermi = band_energy(VALS, ones(nk), nelec, smearing, returnef=true)

    occ = gaussian.(VALS.-efermi, smearing)

    energy_smear = smearing_energy(VALS, ones(nk), efermi, smearing)
        #    println("energy smear $energy_smear")

    energy0 = sum(occ .* VALS0) / nk * 2.0
    energy0 += energy_smear

        #    println("sum occ ", sum(occ), "  ", sum(occ) / (grid[1]*grid[2]*grid[3]))#


    denmat = zeros(Float64, nwan, nwan)
        #    denmatS = zeros(Complex{Float64}, nwan, nwan)


#=        println("cd")
        TEMP = zeros(Complex{Float64}, nwan, nwan) 
        #    println("occ")
        #    for nk = 1:(grid[1]*grid[2]*grid[3])
        #        println(nk)
        #        println(VALS[nk,:])
        #        println(occ[nk,:])
        #    end
        @time for nk = 1:(grid[1]*grid[2]*grid[3])
            for a = 1:nwan
                for i = 1:nwan
                    for j = 1:nwan
                        TEMP[i,j] = VECTS[nk,i,a]' * SK[nk,i,j] * (VECTS[nk,j,a])
                    end
                end
                TEMP =   (TEMP + conj(TEMP))
                #            println("typeof temp ", typeof(TEMP), " " , size(TEMP))
                #            println("typeof denmat ", typeof(denmat), " " , size(denmat))
                denmat += 0.5*occ[nk,a] * real.( TEMP)
                #            denmatS += 0.5*occ[nk,a] * TEMP * SK[nk,:,:]'

            end
        end
        
        denmat = denmat / (grid[1]*grid[2]*grid[3])
        #    denmatS = denmatS / (grid[1]*grid[2]*grid[3])

        println("cd 2")
        TEMP = zeros(Complex{Float64}, nwan, nwan) 
        denmat2 = zeros(Float64, nwan, nwan)

        @time if true
#            occ_p = permutedims(occ, [2 1])
#            VECTS_p = permutedims(VECTS, [3 1 2])
            cVECTS = conj(VECTS)

            tij = 0.0+0.0im
            TEMP[:,:] .= 0.0
            for i = 1:nwan
                for j = 1:nwan
                    tij = 0.0
                    for a = 1:nwan
                        #for nk = 1:(grid[1]*grid[2]*grid[3])
                        tij += sum(occ[:,a].* cVECTS[:,i,a]  .* (VECTS[:,j,a] .* SK[:,i,j]))
#                        end
                        #tij = sum(occ_p[:,nk] .* cVECTS_p[:,nk,i]  .* VECTS_p[:,nk,j])

                    end
                    TEMP[i,j] += tij

                end
            end
            TEMP =   (TEMP + conj(TEMP))
            denmat2 += 0.5* real.( TEMP)
        end
        denmat2 = denmat2 / (grid[1]*grid[2]*grid[3])
        println("diff denmat2 " , sum(abs.(denmat - denmat2)))

=#
#    println("charge")

    if true
        TEMP = zeros(Complex{Float64}, nwan, nwan) 
        pVECTS = permutedims(VECTS, [1,3,2])
        pVECTS_C = conj(pVECTS)

#        cVECTS = conj(VECTS)
        grid3 = prod(grid)
        t_temp = zeros(Complex{Float64},grid3 )
#        println("loop")
#        println(size(occ))
#        println(size(pVECTS_C))
#        println(size(SK))
#        println("nwan $nwan grid $grid grid3 $grid3, nk $nk")
#        count = 0
        @threads for i = 1:nwan
            for j = 1:nwan
                if maxSK[i,j] > 1e-7
#                temp = 0.0
#                for k = 1:nk
#                    for a = 1:nwan
#                        temp += occ[k,a]* pVECTS_C[k,a,i]  * pVECTS[k,a,j] * SK[k,i,j]
#                    end
#                end
#                TEMP[i,j] = temp

  #              if i == 1 && j < 3
#                for k = 1:grid3
 #                   t_temp[k] = sum(occ[k,:].* cVECTS[k,i,:]  .* VECTS[k,j,:])
 #               end
  #              TEMP[i,j] += sum(t_temp .* SK[:,i,j])
 #               else

#                   t_temp[:] = sum(occ[:,:].* cVECTS[:,i,:]  .* VECTS[:,j,:], dims=2)
#                    TEMP[i,j] += sum(t_temp .* SK[:,i,j])

#                t_temp[:] = 
                    TEMP[i,j] += sum( sum(occ[:,:].* pVECTS_C[:,:,i]  .* pVECTS[:,:,j], dims=2) .* SK[:,i,j])


#                    t_temp[:] = sum(occ[:,:].* pVECTS_C[:,:,i]  .* pVECTS[:,:,j], dims=2)
#                    TEMP[i,j] += sum(t_temp .* SK[:,i,j])
 #               else
#                    count += 1
                end

#                end
            end
        end

        TEMP =   (TEMP + conj(TEMP))
        denmat += 0.5* real.( TEMP)
        
        denmat = denmat / (grid[1]*grid[2]*grid[3])
        chargeden = sum(denmat, dims=1)
    end    


    return energy0, efermi, chargeden[:], VECTS, VALS, error_flag
    

end 

function calc_energy_charge_fft(tbc::tb_crys; grid=missing, smearing=0.01)

    etypes = types_energy(tbc.crys)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end
#    println("calc_energy_charge_fft grid $grid")
    hk3, sk3 = myfft_R_to_K(tbc, grid)

    if tbc.scf
        h1 = tbc.tb.h1
        echarge, pot = ewald_energy(tbc)
    else
        h1 = missing
        echarge = 0.0
    end
    eband, efermi, chargeden, VECTS, VALS, error_flag  =  calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1 = h1)


#    println("energy comps $eband $etypes $echarge")
    energy = eband + etypes + echarge

    return energy, efermi, chargeden, VECTS, VALS, error_flag

end


function calc_energy(tbc::tb_crys; smearing=0.01, returnk=false)
    """
    calculate the energy from a kgrid
    """
    kgrid = get_grid(tbc.crys)
    
    return calc_energy(tbc, kgrid, smearing=smearing,returnk=returnk)

end 

function types_energy(tbc::tb_crys)

    return types_energy(tbc.crys)
    
end

function types_energy(tbc::tb_crys_kspace)

    return types_energy(tbc.crys)
    
end

function types_energy(c::crystal)
    
    return types_energy(c.types)
end

function types_energy(types)

    et = 0.0
    for t in types
        et += atoms[t].energy_offset
    end
    return et

end

function calc_energy(h::tb_crys, kgrid; smearing=0.01, returnk=false)
    """
    calculate the energy from a kgrid
    """
    if ismissing(kgrid)
        kgrid = get_grid(tbc.crys)
    end

    etypes = types_energy(h.crys)
    if h.scf
        energy_charge, pot = ewald_energy(h)
    else
        energy_charge = 0.0
    end

    ret = calc_energy_band(h.tb, h.nelec, kgrid, smearing=smearing,returnk=returnk)
#    if returnk
#        eband = ret[1]
#        ek = ret[2]
#        return eband+etypes, ek
#    end

    return ret + etypes + energy_charge

end 

function calc_energy_band(h::tb, nelec, kgrid; smearing=0.01, returnk=false)

    kpts, kweights = make_kgrid(kgrid)

#    println(size(kpts), " ksize " , sum(kweights))
    eigs = calc_bands(h, kpts)
    en, efermi =  band_energy(eigs, kweights, nelec, smearing, returnef=true)
    e_smear  = smearing_energy(eigs, kweights, efermi)

    return en+e_smear
    
end 

function make_kgrid(kgrid)
    kpts = zeros(Float64, prod(kgrid),3)
    kweights = ones(Float64, prod(kgrid))/prod(kgrid) * 2.0
    c=0
    for k1 = 0:kgrid[1]-1
        for k2 = 0:kgrid[2]-1
            for k3 = 0:kgrid[3]-1
                c+=1
                kpts[c,1] = Float64(k1)/Float64(kgrid[1])
                kpts[c,2] = Float64(k2)/Float64(kgrid[2])
                kpts[c,3] = Float64(k3)/Float64(kgrid[3])
            end
        end
    end
    return kpts, kweights

end


#=
function trim_shift(tbc::tb_crys, tol=0.0002)

    energy_orig = tbc.dftenergy
    trim(tbc.tb, tol)
    energy_new = calc_energy(tbc)
    shift = (energy_orig - energy_new)/tbc.nelec
    c_zero = tbc.tb.r_dict[[0,0,0]]
    for i = 1:tbc.tb.nwan
        tbc.tb.H[i,i,c_zero] += shift
    end 
    energy_final = calc_energy(tbc)
    println("trim shift energy_orig $energy_orig energy_trimmed $energy_new energy_final $energy_final")
end
=#
    

function trim(h::tb, tol=0.0002)


    c=0
    keep = []
    for r = 1:h.nr
        if sum(abs.(h.H[:,:,r]) .> tol) > 0 || h.ind_arr[r,:] == [0, 0, 0]
            push!(keep, r)
            c += 1
        end
    end
    println("trim, out of ", h.nr, " keep $c")

    h.H = h.H[:,:,keep]
    if h.nonorth
        h.S = h.S[:,:,keep]
    end
    h.nr = c
    h.ind_arr = h.ind_arr[keep,:]

    h.r_dict = Dict()
    for i in 1:h.nr
        h.r_dict[copy(h.ind_arr[i,:])] = i
    end


    
    
end    

function renormalize_tb(d::dftout, h::tb)
"""
Changes the tight binding matrix elements so that the total band energy 
equals the total atomization energy 


"""

    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)

    eigs = calc_bands(h, d.bandstruct.kpts)
    
    band_en = band_energy(eigs, d.bandstruct.kweights, nval)

    println("original band energy ", band_en)
    
    etot_dft = d.energy

    atomization_energy = etot_dft - etotal_atoms

    println("atomization_enegy: ", atomization_energy)

    shift = (atomization_energy - band_en)/nval

    println("shift $shift nval $nval")
    
    ind = h.r_dict[[0,0,0]]

    println("ind $ind ", h.ind_arr[ind,:])
    
    for i = 1:h.nwan
        h.H[i,i,ind] = h.H[i,i,ind] + shift
    end

    eigs = calc_bands(h, d.bandstruct.kpts)
    band_en_new = band_energy(eigs, d.bandstruct.kweights, nval)

    println("new band_energy ", band_en_new)
    
end

function organizedata(tbc::tb_crys)
    
    return organizedata(tbc.crys, tbc.tb)

end

function organizedata(crys::crystal, h::tb)

    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    data = []

    data_arr = zeros(h.nwan*h.nwan*h.nr,12)
    
    At = transpose(crys.A)

    counter = 0
    ind_dir = Dict()
    for n = 1:h.nwan
        (na, nta, norb) = ind2orb[n]
        for m = 1:h.nwan        
            (ma, mta, morb) = ind2orb[m]            
            if na == ma
                counter += 1
                ind_dir[(n,m)] = counter
            end
        end
    end
    
    data_onsite = zeros(counter,12)
    data_onsite[:,7] .= 10000000000.0
    
    lmn = zeros(3)
    c=0
    c2 = 0

    c_zero = 0
    
    for n = 1:h.nwan
        (na, nta, norb) = ind2orb[n]
        for m = 1:h.nwan        
            (ma, mta, morb) = ind2orb[m]            
            for i in 1:h.nr

                c+=1

                #                println("n m i $n $m $i")
                
                R = h.ind_arr[i,:]

                if R[1] == 0 && R[2] == 0 && R[3] == 0
                    c_zero = i
                end
                
#                dR = At*(-crys.coords[na,:] .+ crys.coords[ma,:] .+ R) 
                dR = At*(crys.coords[na,:] .- crys.coords[ma,:] .+ R) 

                dist = sum(dR.^2)^0.5

                if dist > 1e-7
                    c2+=1
                    if dist > 1e-7
                        lmn[:] = dR/dist
                    else
                        lmn .= 0.0
                    end
                    push!(data, [na, norb, ma, morb, h.H[ n,m,i], dist, lmn])
                    
                    data_arr[c2, 1] = na
                    data_arr[c2, 2] = symbol_dict[norb]              
                    data_arr[c2, 3] = ma
                    data_arr[c2, 4] = symbol_dict[morb]
                    data_arr[c2, 5] = real(h.H[ n,m, i])
                    data_arr[c2, 6] = imag(h.H[ n,m, i])
                    data_arr[c2, 7] = dist
                    data_arr[c2, 8] = lmn[1]
                    data_arr[c2, 9] = lmn[2]
                    data_arr[c2, 10] = lmn[3]
                    data_arr[c2, 11] = real(h.S[ n,m, i])
                    data_arr[c2, 12] = imag(h.S[ n,m, i])
                end
            end
        end
    end
    for (n,m) in keys(ind_dir)
        c2 = ind_dir[(n,m)]
        (na, nta, norb) = ind2orb[n]
        (ma, mta, morb) = ind2orb[m]            

        data_onsite[c2,1] = na
        data_onsite[c2,3] = ma
        data_onsite[c2,2] = symbol_dict[norb]
        data_onsite[c2,4] = symbol_dict[morb]

        data_onsite[c2, 5] = real(h.H[ n,m, c_zero])
        data_onsite[c2, 6] = imag(h.H[ n,m, c_zero])

        data_onsite[c2, 11] = real(h.S[ n,m, c_zero])
        data_onsite[c2, 12] = imag(h.S[ n,m, c_zero])
        
        for x in 1:crys.nat
            for rx = [-1,0,1]
                for ry = [-1,0,1]
                    for rz = [-1,0,1]
                        R = [rx,ry,rz]
                        dR = At*(-crys.coords[na,:] .+ crys.coords[x,:] .+ R) 
                        dist = sum(dR.^2)^0.5
                        if dist > 1e-7 && dist <= data_onsite[c2,7]
                            data_onsite[c2,7] = dist
                            data_onsite[c2,8] = x
                        end
                    end
                end
            end
        end
    end
    return data_onsite, data_arr
                
end

                
            
function myfft_R_to_K(tbc, grid=missing)
"""
Does Fourier Transform R->k. 
"""

    
#    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(tbc.crys)
    
    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end

#    println("grid ", grid)

#    kgrid = make_kgrid(grid)

    nwan = tbc.tb.nwan

    nr = prod(grid)

#    if tbc.tb.nonorth

    ham_R = zeros(Complex{Float64},  nwan, nwan, grid[1], grid[2], grid[3])
    S_R = zeros(Complex{Float64}, nwan, nwan,  grid[1], grid[2], grid[3])
        
    ind = zeros(Int64, 3)
    new_ind = zeros(Int64, 3)

    for c in 1:size(tbc.tb.ind_arr)[1]

        ind[:] = tbc.tb.ind_arr[c,:]
        
        new_ind[1] = mod(ind[1], grid[1])+1
        new_ind[2] = mod(ind[2], grid[2])+1
        new_ind[3] = mod(ind[3], grid[3])+1

#        println(ind, " new ", new_ind)
        
        ham_R[:,:,new_ind[1], new_ind[2], new_ind[3]] += tbc.tb.H[:,:,c]
        S_R[:,:,new_ind[1], new_ind[2], new_ind[3]] += tbc.tb.S[:,:,c]

    end
    
    hamK3 = fft(ham_R, [3,4,5])
    SK3 = fft(S_R, [3,4,5])
    
    return hamK3, SK3

end    
            
function myfft(crys, nonorth, grid, kpts,ham_kS, Sk=missing)

    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(crys)

    
    nwan = size(ham_kS)[1]
    nks  = size(ham_kS)[3]

    if length(size(ham_kS)) == 5 #we don't have to calculate, already in fft form
        println("don't calc")
        ham_k_fftw = ham_kS
        S_k_fftw = Sk

    else
        
        ham_k_fftw = zeros(Complex{Float64}, nwan, nwan, grid[1], grid[2], grid[3])
        S_k_fftw = zeros(Complex{Float64}, nwan, nwan, grid[1], grid[2], grid[3])
        
        for k = 1:nks
            #        println("k $k ", p.bs.kpts[k,:], " " , grid)
            kn = Int.(round.(kpts[k,:] .* grid)).+1
            for i in 1:3
                if kn[i] <= 0
                    kn[i] = kn[i] + grid[i]
                end
            end
            if nonorth
                #            println("kn ", kn)
                ham_k_fftw[:,:,kn[1], kn[2], kn[3]] = ham_kS[:,:,k]
                S_k_fftw[:,:,kn[1], kn[2], kn[3]] = Sk[:,:,k]            
            else
                ham_k_fftw[:,:,kn[1], kn[2], kn[3]] = ham_kS[:,:,k]
            end
        end
        
        println("ham_k_fftw ", size(ham_k_fftw))
    end    


    #does the actual ifft
    hamR3 = ifft(ham_k_fftw, [3,4,5])
    if nonorth
        SR3 = ifft(S_k_fftw, [3,4,5])    
    end


#    println("myfft SR3 0 0 0 ")
#    println(SR3[:,:,1,1,1])
#    println()
    
#    twopi_i = 1.0im*2.0*pi

    
    r_dict = Dict()

    R_grid, R_int_grid, sym_R = get_sym_R(crys, grid)
    ngrid2 = size(R_grid)[1]
    
#    print("R_int_grid")
#    print(R_int_grid[1:20, :])
#    println()
    
    ind_arr = R_grid
    ham_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)

    if nonorth
        S_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)
    end
    
    sym = zero(Float64)
    exp_val = zero(Complex{Float64})
    c=0

    sym_mat = zeros(Float64, ngrid2)

#    println("R_int_grid")
    
    for c = 1:ngrid2
        r_dict[ind_arr[c,:]] = c
        rint = R_int_grid[c,:]
        #println("$c ", R_int_grid[c,:] )
        #put ifft'd arrays into new arrays
        for nw_1 = 1:nwan
            a1 = wan_atom[nw_1]
            for nw_2 = 1:nwan
                a2 = wan_atom[nw_2]
                if nonorth
                    S_r[nw_1,nw_2,c] += SR3[nw_1,nw_2,rint[1], rint[2], rint[3]] * sym_R[a1,a2,c]
                end
                ham_r[nw_1,nw_2,c] += hamR3[nw_1,nw_2,rint[1], rint[2], rint[3]] * sym_R[a1,a2,c]
#                println("adding $a1 $a2 $c ", hamR3[a1,a2,rint[1], rint[2], rint[3]], " " , sym_R[a1,a2,c]
#                        )
            end
        end
    end

    if nonorth
        return ham_r, S_r, r_dict, ind_arr
    else
        return ham_r, r_dict, ind_arr
    end             

             
end

function get_sym_R(crys, grid, sss = 1.0)


    grid2 = [0,0,0]
    for i = 1:3
        if grid[i] == 1
            grid2[i] = 3
        elseif grid[i]%2 == 0
            grid2[i] = Int(grid[i]/2)+4
        else        
            grid2[i] = Int((grid[i]+1)/2)+4
        end
    end

 #   println("grid2: " , grid2)

    A = crys.A
    At = transpose(crys.A)
#    At = crys.A
    c=0


    r = zeros(Float64, 1,3)
    R = zeros(Float64, 1, 3)
    rint = zeros(Int64, 1, 3)

    dR = zeros(Float64, 1, 3)
    
    dist_arr = zeros(5^3)

    ngrid2 = (grid2[1]*2+1)*(grid2[2]*2+1)*(grid2[3]*2+1)
#    println("ngrid2: $ngrid2")
    sym_R = zeros(Float64,  crys.nat, crys.nat, ngrid2)

    R_grid = zeros(Int64, ngrid2, 3)
    R_int_grid = zeros(Int64, ngrid2, 3)

    
    for a = 1:crys.nat
        for b = 1:crys.nat        
            c=0
            for r1 = -grid2[1]:(grid2[1])
                r[1] = r1
#                rint[1]= Int.(-r[1] .+ 1)
                for r2 = -grid2[2]:(grid2[2])
                    r[2] = r2
#                    rint[2]= Int.(-r[2] .+ 1)
                    for r3 = -grid2[3]:(grid2[3])

                        c+=1

                        r[3] = r3
#                        rint[3]= Int.(-r[3] .+ 1)

#####                   #prepare Rgrids
                        if (a == 1 && b == 1)

                            R_grid[c,:] = [r1,r2,r3]
                            
                            rint[:] = Int.(r[:] .+ 1)
                            for i = 1:3
                                if rint[i] <= 0
                                    rint[i] = rint[i] + grid[i]
                                elseif rint[i] <= 0
                                    rint[i] = rint[i] + grid[i]
                                elseif rint[i] > grid[i]
                                    rint[i] = rint[i] - grid[i]
                                end
                            end
                            R_int_grid[c,:] = rint[:]
                        end
###############################
                        
#                        dR[:] = At*(-crys.coords[a,:] + crys.coords[b,:] + r) A
                        dR[:] = (crys.coords[[a],:] - crys.coords[[b],:] + r) * A
                        dist = sum(dR.^2)^0.5
                        cd = 0
#                        if dist < 7.0
#                            println("AP $a $b [$r1 $r2 $r3] $dist $dR $r")
#                        end
                        #                        for R1 = [0]
#                            for R2 = [0]                            
#                        for R1 = [-2*grid[1], -grid[1], 0, grid[1], 2*grid[1]]
#                            for R2 = [-2*grid[1], -grid[2], 0, grid[2], 2*grid[2]]
#                                for R3 = [-2*grid[3], -grid[3], 0, grid[3], 2*grid[3]]
                        for R1 = -2:2
                            for R2 = -2:2
                                for R3 = -2:2
                                    R[1,:] = [R1*grid[1] R2*grid[2] R3*grid[3]]
                                    cd += 1
                                    dR[:] = (crys.coords[[a],:] - crys.coords[[b],:] + r + R) * A
                                    dist_arr[cd] = sum(dR.^2)^0.5
                                end
                            end
                        end
                        if isapprox(dist, minimum(dist_arr))
                            cd = 0
                            for i = 1:125
                                if isapprox(dist, dist_arr[i])
                                    cd += 1
                                end
                            end
                            sym_R[a,b,c] = 1.0/Float64(cd)
                        end

                        
                    end
                end
            end
        end
    end


    for a = 1:crys.nat
        for b = 1:crys.nat        
#            for i in 1:ngrid2
#                println("sym_R[$a,$b,$i] = ", sym_R[a,b,i], "  R= ", R_grid[i,:])                
#            end
            println("SUM $a $b : ", sum(sym_R[a,b,:]))
        end
    end

    keep = []
    nkeep = 0
    for i = 1:ngrid2 
        if maximum(sym_R[:,:,i]) > 1e-5
            push!(keep, i)
            nkeep += 1
        end
    end
    println("nkeep : $nkeep")
    R_grid = R_grid[keep,:]
    R_int_grid = R_int_grid[keep,:]
    sym_R = sym_R[:,:,keep]

    return R_grid, R_int_grid, sym_R

end    

function tb_indexes(d::dftout)
    crys = d.crys
    return tb_indexes(crys)
end

function tb_indexes(crys::crystal)

    semicore = []
    wan = []

    wan_atom = Dict()
    atom_wan = Dict()
    
    n = 1
    n2 = 1
    c=0
    for (at_num, t) in enumerate(crys.types)
#        println(t, " ",atoms[t].nsemicore/2," ",atoms[t].nwan/2)
        
        n2 += Int(atoms[t].nsemicore/2)
        if n2 > n
            append!(semicore, collect(n:n2-1))
        end
        n = n2
        n2 += Int(atoms[t].nwan/2)
        if n2 > n
            append!(wan, collect(n:n2-1))
        end
        n = n2

        atom_wan[at_num] = []
        for nw = 1:Int(atoms[t].nwan/2)
            c += 1
            wan_atom[c] = at_num
            append!(atom_wan[at_num], c)
        end
        
    end

    nwan = length(wan)
    nsemi = length(semicore)

    
    
    return wan, semicore, nwan, nsemi, wan_atom, atom_wan
end

function symm_by_orbitals(crys::crystal, mat)
    c=0
    for (at_num, t) in enumerate(crys.types)
        for o in atoms[t].orbitals
            if o == :s
                c += 1
            elseif o == :p
                t = sum(mat[c+1:c+3])/3.0
                mat[c+1:c+3] .= t
                c += 3
            elseif o == :d
                t = sum(mat[c+1:c+5])/5.0
                mat[c+1:c+5] .= t
                c += 5
            elseif o == :f
                t = sum(mat[c+1:c+7])/7.0
                mat[c+1:c+7] .= t
                c += 7
            else
                println("WARNING symm_by_orbitals expecting orbital got : $o")
                c += 1
            end
        end
    end
    return mat
end


function plot_compare_dft(tbc, bs; tbc2=missing, align_min=true)

    kpts = bs.kpts
    kweights = bs.kweights

    vals = calc_bands(tbc.tb, kpts)
    

    nelec = tbc.nelec

    nsemi = 0
    for t in tbc.crys.types
        atom = atoms[t]
        nsemi += atom.nsemicore
    end
    nsemi = Int64(nsemi / 2)

    efermi_dft = calc_fermi(bs.eigs[:,1:end], kweights, nelec+nsemi*2)
    efermi_tbc = calc_fermi(vals, kweights, nelec)

    if !ismissing(tbc2)
        vals2 = calc_bands(tbc2.tb, kpts)
        efermi_tbc2 = calc_fermi(vals2, kweights, nelec)
        println("efermi_dft $efermi_dft efermi_tbc $efermi_tbc efermi_tbc2 $efermi_tbc2")

        en_dft = band_energy(bs.eigs[:,1:end], kweights, nelec+nsemi)
        en_1 = band_energy(vals, kweights, nelec)
        en_2 = band_energy(vals2, kweights, nelec)
        println("energy_dft $en_dft en_1 $en_1 en_2 $en_2")

        e_smear_dft  = smearing_energy(bs.eigs[:,nsemi+1:end], kweights, efermi_dft)
        e_smear_1  = smearing_energy(vals, kweights, efermi_tbc)
        e_smear_2  = smearing_energy(vals2, kweights, efermi_tbc2)
        
        println("smear_dft $e_smear_dft sm_1 $e_smear_1 sm_2 $e_smear_2")

    end


    nk = size(kpts)[1]

#    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(tbc.crys)
    println("nsemi ", nsemi)
    
    if align_min
        a = minimum(bs.eigs[:,nsemi+1:end])
        b = minimum(vals)
        plot(bs.eigs .- a, "-c", LineWidth=1.5)    
        plot(vals .- b, "-b", LineWidth=1.0)

        plot([0, nk], [efermi_tbc - b, efermi_tbc - b], "--r")
        plot([0, nk], [efermi_dft - a, efermi_dft - a], ":c")

        if !ismissing(tbc2)
            c = minimum(vals2)
            plot(vals2 .- c, "--m.", MarkerSize=3)
        end


    else
    
        vbmD, cbmD = find_vbm_cbm(bs.eigs[:,1:end] , efermi_dft)
        vbm, cbm = find_vbm_cbm(vals , efermi_tbc)

        plot(bs.eigs[:,1:end] .- vbmD , "-c", LineWidth=1.5)    
        plot(vals .- vbm, "-b", LineWidth=1.0)


        if !ismissing(tbc2)
            vbm2, cbm2 = find_vbm_cbm(vals2 , efermi_tbc2)
            plot(vals2 .- vbm2,  "--m.", MarkerSize=3)
        end


        plot([0, nk], [efermi_tbc, efermi_tbc] .- vbm, "--k")
        plot([0, nk], [efermi_dft, efermi_dft] .- vbmD, ":r")

        ylim(minimum(vals .- vbm)*1.05, maximum(vals .- vbm) * 1.05)

    end

end

function find_vbm_cbm(eigs, fermi)

    vbm = -100000.0
    cbm = 100000.0
    for eig in eigs[:]
        if eig < fermi && eig  > vbm
            vbm = eig
        end

        if eig > fermi && eig  < cbm
            cbm = eig
        end
    end
    return vbm, cbm
end


function ewald_energy(tbc::tb_crys, delta_q=missing)
    
    gamma = tbc.gamma 
    crys = tbc.crys
    
    if ismissing(delta_q)
        delta_q =  get_dq(crys , tbc.eden)
    end
    
    return ewald_energy(crys, gamma, delta_q)

end

function ewald_energy(tbc::tb_crys_kspace, delta_q=missing)
    
    gamma = tbc.gamma 
    crys = tbc.crys
    
    if ismissing(delta_q)
        delta_q =  get_dq(crys , tbc.eden)
    end
    
    return ewald_energy(crys, gamma, delta_q)

end

function ewald_energy(crys::crystal, gamma, delta_q::Array{Float64,1})
    
    T = typeof(crys.coords[1,1])
    pot = zeros(T, crys.nat, crys.nat)

    for i = 1:crys.nat
        for j = 1:crys.nat
            pot[i,j] = gamma[i,j] * delta_q[i] * delta_q[j] 
        end
    end


    energy = 0.5*sum(pot)

#    println("ewald_energy ", energy, " " , delta_q, " ", gamma[1,1], " ", gamma[1,2], " ", gamma[2,1], " ", gamma[2,2])

    return energy, pot

end

function get_neutral_eden(tbc::tb_crys)

    return get_neutral_eden(tbc.crys, tbc.tb.nwan)

end

function get_neutral_eden(crys::crystal, nwan=missing)

    if ismissing(nwan)
        ind2orb, orb2ind, etotal, nval = orbital_index(crys)
        nwan = length(keys(ind2orb))
    end

    eden = zeros(nwan)
    counter = 0
    for (i, t) in enumerate(crys.types)
        at = atoms[t]
        z_ion = at.nval
        still_need = z_ion
        for o in at.orbitals
            if o == :s
                counter += 1
                if still_need <= 2.0 && still_need > 1e-5
                    eden[counter] = still_need/2.0
                    still_need = 0.0
                elseif still_need >= 2.0 && still_need > 1e-5
                    eden[counter] = 1.0
                    still_need = still_need - 2.0
                else
                    counter += 1
                end
            elseif o == :p
#                println("p")
                if still_need <= 6.0 && still_need > 1e-5
                    eden[counter+1] = (still_need/2.0)/3.0
                    eden[counter+2] = (still_need/2.0)/3.0
                    eden[counter+3] = (still_need/2.0)/3.0
                    still_need = 0.0
                    counter += 3
                elseif still_need >= 6.0 && still_need > 1e-5
                    eden[counter+1] = 1.0
                    eden[counter+2] = 1.0
                    eden[counter+3] = 1.0
                    still_need = still_need - 6.0
                    counter += 3
                else
                    counter += 3
                end
            elseif  o == :d
                if still_need <= 10.0 && still_need > 1e-5
                    eden[counter+1] = (still_need/2.0)/5.0
                    eden[counter+2] = (still_need/2.0)/5.0
                    eden[counter+3] = (still_need/2.0)/5.0
                    eden[counter+4] = (still_need/2.0)/5.0
                    eden[counter+5] = (still_need/2.0)/5.0
                    still_need = 0.0
                    counter += 5
                elseif still_need >= 10.0 && still_need > 1e-5
                    eden[counter+1] = 1.0
                    eden[counter+2] = 1.0
                    eden[counter+3] = 1.0
                    eden[counter+4] = 1.0
                    eden[counter+5] = 1.0
                    still_need = still_need - 10.0
                    counter += 5
                else
                    counter += 5
                end
            else
                println("bad orbital get_neutral_eden $o")
            end
        end
    end
    return eden

end

function get_dq(tbc::tb_crys_kspace)
    return get_dq(tbc.crys, tbc.eden)
end

function get_dq(tbc::tb_crys)
    return get_dq(tbc.crys, tbc.eden)
end


function get_dq(crys::crystal, chargeden::Array{Float64,1})


    e_den = zeros(Float64, crys.nat)
    z_ion = zeros(Float64, crys.nat)

    counter = 0
    for (i, t) in enumerate(crys.types)
        at = atoms[t]
        z_ion[i] = at.nval
        for o = 1:Int64(at.nwan/2)
            counter += 1
            e_den[i] += chargeden[counter]  * 2.0
        end
    end
    dq = -z_ion + e_den

    dq = dq .- sum(dq)/crys.nat #charge sum rule

#    println("e_den ", e_den)
#    println("z_ion ", z_ion)
#    println("dq ", dq)


    return dq

end

function get_h1(tbc::tb_crys)
    return get_h1(tbc, tbc.eden)
end


function get_h1(tbc::tb_crys, chargeden::Array{Float64,1})

    dq = get_dq(tbc.crys, chargeden)

    gamma = tbc.gamma
    
    epsilon = gamma * dq

    h1 = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    o1 = 1
    for i = 1:tbc.crys.nat
        at1 = atoms[tbc.crys.types[i]  ]
        nw1 = Int64(at1.nwan/2)
        o2 = 1
        for j = 1:tbc.crys.nat
            at2 = atoms[tbc.crys.types[j]]
            nw2 = Int64(at2.nwan/2)
            for c1 = o1:o1+nw1-1
                for c2 = o2:o2+nw2-1
                    h1[c1,c2] = 0.5 * (epsilon[i] + epsilon[j])
                end
            end
            o2 += nw2

        end
        o1 += nw1
    end

    return 0.5*(h1 + h1'), dq

end


function get_h1(tbc::tb_crys_kspace)
    return get_h1(tbc, tbc.eden)
end


function get_h1(tbc::tb_crys_kspace, chargeden::Array{Float64,1})

    dq = get_dq(tbc.crys, chargeden)

    gamma = tbc.gamma
    
    epsilon = gamma * dq

    h1 = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    o1 = 1
    for i = 1:tbc.crys.nat
        at1 = atoms[tbc.crys.types[i]  ]
        nw1 = Int64(at1.nwan/2)
        o2 = 1
        for j = 1:tbc.crys.nat
            at2 = atoms[tbc.crys.types[j]]
            nw2 = Int64(at2.nwan/2)
            for c1 = o1:o1+nw1-1
                for c2 = o2:o2+nw2-1
                    h1[c1,c2] = 0.5 * (epsilon[i] + epsilon[j])
                end
            end
            o2 += nw2

        end
        o1 += nw1
    end

    return 0.5*(h1 + h1'), dq

end

function get_energy_electron_density_kspace(tbcK::tb_crys_kspace; smearing = 0.01)

    bandenergy, eden, VECTS, VALS, efermi, error_flag = get_energy_electron_density_kspace(tbcK.tb, tbcK.nelec, smearing=smearing)

    if tbcK.scf
#        h1 = tbc.tb.h1
        echarge, pot = ewald_energy(tbcK)
    else
#        h1 = missing
        echarge = 0.0
    end

    etypes = types_energy(tbcK.crys)

    energy_smear = smearing_energy(VALS, tbcK.tb.kweights, efermi, smearing)
    println("CALC ENERGIES $etypes $echarge $bandenergy $energy_smear = ", bandenergy + etypes + echarge + energy_smear)

    return bandenergy + etypes + echarge + energy_smear, eden, VECTS, VALS, error_flag


end

function get_energy_electron_density_kspace(tb_k::tb_k, nelec; smearing = 0.01)

    temp = zeros(Complex{Float64}, tb_k.nwan, tb_k.nwan)
    denmat = zeros(Float64, tb_k.nwan, tb_k.nwan)

    VALS = zeros(Float64, tb_k.nk, tb_k.nwan)
    VALS0 = zeros(Float64, tb_k.nk, tb_k.nwan)
    VECTS = zeros(Complex{Float64}, tb_k.nk, tb_k.nwan, tb_k.nwan)
    SK = zeros(Complex{Float64}, tb_k.nk, tb_k.nwan, tb_k.nwan)

    error_flag = false

    for k in 1:tb_k.nk
        try
            vects, vals, hk, sk, vals0 = Hk(tb_k, tb_k.K[k,:])
            VALS[k,:] = vals
            VECTS[k,:,:] = vects
            SK[k,:,:] = sk
            VALS0[k,:] = vals0
        catch
            println("warning, get_energy_electron_density_kspace error")
            error_flag = true
        end
    end
    energy, efermi, occs = band_energy(VALS, tb_k.kweights, nelec, smearing, returnboth=true)
    energy0 = sum(occs .* VALS0 .* tb_k.kweights) / sum(tb_k.kweights) * 2.0

    for k in 1:tb_k.nk
        for a = 1:tb_k.nwan
            for i = 1:tb_k.nwan
                for j = 1:tb_k.nwan
                    temp[i,j] = VECTS[k,i,a]' * SK[k,i,j] * VECTS[k,j,a]
                end
            end
            temp = temp + conj(temp)
            denmat += 0.5 * occs[k,a] * real.(temp) * tb_k.kweights[k]
        end
    end
    electron_den = sum(denmat, dims=1) / sum(tb_k.kweights)
    return energy0, electron_den[:], VECTS, VALS, efermi, error_flag

end


end #end module
