##include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### QE specific 
module AtomicProj
"""
Scripts to analyze atomic projections
"""
#export load_xml
using LinearAlgebra
using DelimitedFiles

#include("Atomdata.jl")
using ..Atomdata:atoms
using JLD
#using ExXML
using XMLDict
using GZip

using FFTW
using ..DFToutMod
#using ..CrystalMod:crystal
#using ..Utility:arr2str
#using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..Utility:parse_str_ARR_complex
using ..DFToutMod:bandstructure
using ..TB:make_tb
using ..TB:make_tb_crys
using ..BandTools:band_energy
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..TB:write_tb_crys
#using ..TB:orbital_index
using ..TB:trim
using ..QE:loadXML_bs
using ..QE:loadXML
using ..TB:calc_energy
using ..TB:types_energy
using ..CrystalMod:crystal
using ..DFT:runSCF
using ..CrystalMod:get_grid
using ..TB:myfft
using ..TB:tb_indexes
using ..TB:symm_by_orbitals

using ..TB:tb_crys_kspace
using ..TB:make_tb_crys_kspace
using ..TB:make_tb_k
using ..TB:write_tb_crys_kspace


include("Commands.jl")



mutable struct proj_dat

    bs::bandstructure
    natwfc::Int64
    proj::Array{Complex{Float64}, 3}
    overlaps::Array{Complex{Float64}, 3}    
end

#print proj dat
Base.show(io::IO, d::proj_dat) = begin
    println(io,"projection data: nbnd = ", d.bs.nbnd, "; nkpts = ", d.bs.nks, "; natwfc = ", d.natwfc)
end   


function write_to_file(str, filename, directory="./")

    f=open("$directory/$filename", "w")
    write(f, str)
    close(f)

end

function makedict_proj(savedir)

    if isfile(savedir*"/atomic_proj.xml")
        filename=savedir*"/atomic_proj.xml"
    elseif isfile(savedir*"/atomic_proj.xml.gz")
        filename=savedir*"/atomic_proj.xml.gz"
    else
        println("error warning no atomic_proj.xml.gz or atomic_proj.xml")
        filename = missing
    end
    
    f = gzopen(filename, "r")
    fs = read(f, String)
    close(f)
    
    d = xml_dict(fs)

    return d
end

function make_projwfcx(prefix, tmpdir)
"""
makes the proj input file
"""
    c_dict = make_commands(1)
    julia_dir = c_dict["juliadir"]

    template_file=open("$julia_dir/template_inputs/template.proj")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    
    return temp
    
end


function run_projwfcx(projfile="proj.in"; directory="./", nprocs=1)
"""
run projwfc.x QE command
"""

    c_dict = make_commands(nprocs)
    proj = c_dict["proj"]
    command = `$proj $directory/$projfile `
    println("projwfc.x command")
    println(command)
    println()

    try
        println("Running projwfc.x")
        s = read(command, String)
        f = open(directory*"/"*projfile*".out", "w")
        write(f, s)
        close(f)
        
        println("Ran projwfc.x")
        println()
    catch
        println("Failed to run projwfc.x ")
        #println(s)
        return -1
    end

end

function makeOG(prefix, tmpdir )
"""
makes the OG file
"""
    template_file=open("template_inputs/template_og.in")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    
    return temp
    
end

function run_og(filename="og.in";  directory="./", nprocs=1)
"""
run open_grid.x command
"""
    c_dict = make_commands(nprocs)
    og = c_dict["og"]
    command = `$og $directory/$filename`
    println("opengrid command")
    println(command)
    try
        println("Running open_grid.x")

        s = read(command, String)
        f = open(directory*"/"*filename*".out", "w")
        write(f, s)
        close(f)
        
        println("Ran open_grid.x")
        return 0
    catch
        println("Failed to run open_grid.x")
        return -1
    end
        
end

function projwfx_workf(dft::dftout; directory="./", nprocs=1, freeze=true, writefile="projham.xml",writefilek="projham_K.xml", skip_og=true, skip_proj=true, shift_energy=true, cleanup=true, skip_nscf=true, localized_factor = 0.05, only_kspace=false, screening = 1.0)
"""
All steps to run wannier90 from converged QE SCF calculation

Steps:

1) open_grid scf calc 
2) run projwfc.x
2a) load output
3) make projected ham in k-space
4) do fourier transform
5) optionally write results to file
"""
    prefix=dft.prefix
    outdir=dft.outdir
    println()
    println("projwfc_workflow---------------------------------------------------------------------------")
    println()
    println("Step #0 nscf---------------------------------------------------------------------------------")
    s=""
    prefix_orig = deepcopy(prefix)

#    try

    olddir = "$directory/$prefix.save"
    nscfdir = "$directory/$prefix.nscf.save"

    function run_nscf()

        println("doing nscf")
        if isdir(nscfdir)
            println("removing original nscf directory $nscfdir")
            rm(nscfdir,recursive=true)
        end
        mkpath(nscfdir)
#        println("after mkdpath")

        files = readdir(olddir)
        for f in files
            if length(f) > 3 && f[end-2] == 'U' && f[end-1] == 'P' && f[end] == 'F'
#                println("found $f , copy")
                cp("$olddir/$f", "$nscfdir/$f")            
            end
        end

        try
            if isfile("$olddir/charge-density.dat")
                cp("$olddir/charge-density.dat",  "$nscfdir/charge-density.dat")
            end
            if isfile("$olddir/data-file-schema.xml")
                cp("$olddir/data-file-schema.xml", "$nscfdir/data-file-schema.xml" )
            end
            if isfile("$olddir/data-file-schema.xml.gz")
                cp("$olddir/data-file-schema.xml.gz", "$nscfdir/data-file-schema.xml.gz" )

                tounzip = "$nscfdir/data-file-schema.xml.gz"
                command = `gunzip $tounzip`
                s = read(command, String)
                println("gunzipped $tounzip")
            end
            if isfile("$olddir/charge-density.dat.gz")
                cp("$olddir/charge-density.dat.gz", "$nscfdir/charge-density.dat.gz" )

                tounzip = "$nscfdir/charge-density.dat.gz"
                command = `gunzip $tounzip`
                s = read(command, String)
                println("gunzipped $tounzip")

            end


        catch
            println("missing charge density or xml file, cannot run nscf")
        end
            
        crys = dft.crys
        tot_charge = dft.tot_charge

        if only_kspace
            calc = "nscf-sym"
        else
            calc = "nscf"
        end

        try
            dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=2, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge)
        catch
            println()
            println("first nscf failed, trying backup nscf with fewer extra bands")
            println()
            try
                dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=1, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge, use_backup=true)
            catch
                println("try 2")
                dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=-1, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge, use_backup=true)
            end
        end

        dft_nscf.energy = dft.energy
        dft_nscf.energy_smear = dft.energy_smear
        dft_nscf.forces = dft.forces
        dft_nscf.stress = dft.stress        
        
        prefix = prefix*".nscf"

        return dft_nscf

    end
    
    if !(skip_nscf) || !(isdir(nscfdir)) ||  ( !isfile(nscfdir*"/atomic_proj.xml") && !isfile(nscfdir*"/atomic_proj.xml.gz"))

        dft_nscf = run_nscf()

    else
        if (isdir(nscfdir))
            crys = dft.crys
            #            dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=true, nprocs=nprocs, skip=false, calculation="nscf")
            dft_nscf = loadXML("$directory/$prefix.nscf.save")            
            dft_nscf.energy = dft.energy
            dft_nscf.forces = dft.forces
            dft_nscf.stress = dft.stress        
            prefix = prefix*".nscf"
            println("loaded nscf aaaaaaaaaaaaaaaaaaaaaa")
        else
            crys = dft.crys
            dft_nscf = dft
        end



    end
        
#    catch
#        println("failed nscf")
#        println(s)
#    end

#    return


    newprefix = prefix
    
    useog = false
    if useog
        println("Step #1 og---------------------------------------------------------------------------------")
        println()

        newprefix=prefix*"_open"
        
        function og()
            og= makeOG(prefix, outdir)
            write_to_file(og, "og.in", directory)
            run_og("og.in", directory=directory, nprocs=nprocs)
        end
        
        if skip_og
            try
                dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
                println("using precomputed opengrid dir")
            catch
                println("precomputinga opengrid dir")
                og()
                dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
                
            end
        else
            og()
            dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
        
        end

#        kgrid = copy(dft.bandstruct.kgrid)
#        dft_nscf.bandstruct = dftopen_bs
#        dft_nscf.bandstruct.kgrid = kgrid

    else
        println("DON'T USE OG, since it works poorly")
        dftopen_bs = loadXML_bs("$directory/$newprefix.save")
    end

    

        
    println()
    println("Step #2 run projwfc.x----------------------------------------------------------------------")
    println()

    A = dft.crys.A             #needed to convert kpoints to crystal units
    a1=sum(A[1,:].^2)^0.5
    B = inv(A./a1)'
    
    function proj()
        projstr = make_projwfcx(newprefix, directory)
        write_to_file(projstr, "proj.in", directory)
        run_projwfcx("proj.in", directory=directory, nprocs=nprocs)
    end
    
    if skip_proj
        try
            p = loadXML_proj("$directory/$newprefix.save", B)
            println("using precomputed proj file")
        catch
            println("computing proj file")
            proj()
            p = loadXML_proj("$directory/$newprefix.save", B)
        end
    else
        proj()
        p = loadXML_proj("$directory/$newprefix.save", B)
        
    end

    println("P OVERLAPS ", size(p.overlaps))
    
    #check if p matches dft_nscf
    if sum(abs.(dft_nscf.bandstruct.eigs[1,:] - p.bs.eigs[1,:])) > 1e-5
        println("warning, projected eigs and dft_nscf eigs do not match ", sum(abs.(dft_nscf.bandstruct.eigs[1,:] - p.bs.eigs[1,:])))
        println("rerun nscf")
        prefix = deepcopy(prefix_orig)
        dft_nscf = run_nscf()
        println("rerun proj")
        proj()
        p = loadXML_proj("$directory/$newprefix.save", B)

    end


    println()
    println("Step #3 ham_k------------------------------------------------------------------------------")
    println()
    projection_warning = false
    if freeze
        band_froz = Int64(round(dft_nscf.bandstruct.nelec/2))+1
        en_froz = minimum(dft_nscf.bandstruct.eigs[:,band_froz])
        en_froz = max(en_froz, dft.bandstruct.efermi + 0.05)
        println("en_froz: ", en_froz, " band froz $band_froz")
        println("efermi energy ", dft.bandstruct.efermi)
        println("nk ", size(dft_nscf.bandstruct.kpts), "     ", size(dft_nscf.bandstruct.kweights))
        ham_k, EIG, Pmat, Nmat, VAL, projection_warning = AtomicProj.create_tb(p, dft_nscf, energy_froz=en_froz+.05, shift_energy=shift_energy);
        #A,B,C = AtomicProj.create_tb(p, dft, energy_froz=en_froz+.01); 
        #return A,B,C        
    else
        ham_k, EIG, Pmat, Nmat, VAL, projection_warning = AtomicProj.create_tb(p, dft_nscf);
    end
    println()
    println("Step #4 ham_r------------------------------------------------------------------------------")
    println()

    println("after create P OVERLAPS ", size(p.overlaps))


    #get tight binding
    tbck = prepare_ham_k(p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=true, localized_factor = localized_factor, screening=screening)

    println("after prepare P OVERLAPS ", size(p.overlaps))


    if !ismissing(writefilek)
        println("Step #5-0 write_tb_k-------")
        write_tb_crys_kspace("$directory/$writefilek", tbck);
    end

    if !only_kspace
        tbc = AtomicProj.get_ham_r(tbck)
        if !ismissing(writefile)
            println()
            println("Step #5 write_tb---------------------------------------------------------------------------")
            println()
            write_tb_crys("$directory/$writefile", tbc);
        end
    else
        tbc = missing
    end



#    tbc = AtomicProj.get_ham_r(p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=true, localized_factor = localized_factor);

    #tbc = AtomicProj.get_ham_r( p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=false);    


    #tbc = AtomicProj.get_ham_r_slow(p, dft_nscf, dft.bandstruct.kgrid ,ham_k, nonorth=true);    


#    trim(tbc.tb, 0.00005)
    
    println()
    println("Done proj----------------------------------------------------------------------------------")
    println()

    if cleanup
        println("clean: warning, removing wavefunctions to save space")

        newdir = "$directory/$newprefix.save"
        olddir = "$directory/$prefix.save"
        parentdir = "$directory"

        if (prefix_orig != prefix) && (newprefix != prefix_orig)
            origdir = "$directory/$prefix_orig.save"
            toclean = [newdir,olddir,origdir]
        else
            toclean = [newdir,olddir]
        end

        if (isdir(nscfdir))
            toclean = [toclean; nscfdir]
        end
        if (isdir(parentdir))
            toclean = [toclean; parentdir]
        end


        tot = 0
        for d in toclean
            files=read(`ls $d/`, String)
            for f in split(files, "\n")
                if occursin("wfc", f)
                    if tot <= 2
                        println("removing $d/$f")
                    elseif tot==3
                        println("removing $d/$f")
                        println("...")
                    end
                    cmd=`\rm $d/$f`
                    s = read(cmd, String)
                    tot += 1
                end
            end
        end
        println("clean-up done: $tot wfc files deleted")
    else
        println("don't clean wfcs")
    end        
    
    return tbc, tbck, projection_warning
    
end
    



function loadXML_proj(savedir, B=missing)
    
    d= makedict_proj(savedir)

    da = d["ATOMIC_PROJECTIONS"] #has all the real data in the xml

    nbnd = parse(Int, da["HEADER"]["NUMBER_OF_BANDS"][""])
    nk = parse(Int, da["HEADER"]["NUMBER_OF_K-POINTS"][""])    
    nspin = parse(Int, da["HEADER"]["NUMBER_OF_SPIN_COMPONENTS"][""])
    natwfc = parse(Int, da["HEADER"]["NUMBER_OF_ATOMIC_WFC"][""])
    nelec = parse(Float64, da["HEADER"]["NUMBER_OF_ELECTRONS"][""])
    efermi = parse(Float64, da["HEADER"]["FERMI_ENERGY"][""])                

    units_energy = da["HEADER"]["UNITS_FOR_ENERGY"][:UNITS]
    units_kpt = da["HEADER"]["UNITS_FOR_K-POINTS"][:UNITS]    

    println(units_energy, " " , units_kpt)
    println("nbands: ", nbnd," nk: ", nk)

    kpts = parse_str_ARR_float(da["K-POINTS"][""])
    kpts = reshape(kpts, 3,nk)'

    if !(ismissing(B))
        kpts = kpts * inv(B)
    end
    
    weights = parse_str_ARR_float(da["WEIGHT_OF_K-POINTS"][""])
    
    dae = da["EIGENVALUES"]
    
    eigs = zeros(nk,nbnd)
    for k in 1:nk
        eigs[k, :] = parse_str_ARR_float(dae["K-POINT."*string(k)]["EIG"][""])
        
    end

    dap = da["PROJECTIONS"]

    proj = zeros(Complex{Float64}, nk, natwfc, nbnd)
    for k in 1:nk
        t = dap["K-POINT."*string(k)]
        for a in 1:natwfc
            proj[k, a, :] =  parse_str_ARR_complex(t["ATMWFC."*string(a)][""])
        end
    end


    dao = da["OVERLAPS"]
    overlaps = zeros(Complex{Float64}, nk, natwfc, natwfc)

    for k in 1:nk
        t = dao["K-POINT."*string(k)]
        overlaps[k, :,:] = reshape(parse_str_ARR_complex(t["OVERLAP.1"][""]), natwfc, natwfc)'
    end

#    println(overlaps)

    bs = DFToutMod.makebs(nelec, efermi, kpts, weights, [0,0,0], eigs)
    
    return make_proj(bs, proj, overlaps)
    
end

function make_proj(bs, proj, overlaps)

    if bs.nks != size(proj)[1] || bs.nks != size(overlaps)[1]
        error("make_proj something wrong nks ", bs.nks," ",size(proj)[1]," ",size(overlaps)[1])
    end
    if bs.nbnd != size(proj)[3]
        error("make_proj something wrong nbnd ", bs.nbnd," ",size(proj)[3])
    end
    if size(proj)[2] != size(overlaps)[2] ||  size(proj)[2] !=	size(overlaps)[3]
        error("make_proj something wrong in natwfc ",  size(proj)[2]," ",size(overlaps)[2]," " ,size(overlaps)[3])
    end

    return proj_dat(bs, size(proj)[2], proj, overlaps)
    
end



#function create_temp(p::proj_dat, d::dftout; energy_froz=missing, nfroz=0, shift_energy=true)
##
#
#    return
#end

function create_tb(p::proj_dat, d::dftout; energy_froz=missing, nfroz=0, shift_energy=true)




    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)

    efermi_dft = d.bandstruct.efermi

    #decide semicore
    p2 = zeros(p.bs.nbnd)
    INDSEMI = zeros(Int64, p.bs.nks, nsemi)
    for k = 1:p.bs.nks
        p2[:] = real.(sum(p.proj[k,semicore,:] .* conj.(p.proj[k,semicore,:]) , dims=1))
        INDSEMI[k, :] = sortperm(p2, rev=true)[1:nsemi]
        if sum(p2) < nsemi - 0.5
            println("warning, identify semicore")
        end
        if k < 10
            println("p2 ", p2)
            println(k, " indsemi " , INDSEMI[k, :])
        end

    end
    println("INDSEMI ", INDSEMI)

    NBND = p.bs.nbnd - nsemi
    EIGS = zeros(p.bs.nks, NBND)
    PROJ = zeros(Complex{Float64}, p.bs.nks, nwan, NBND)

    for k = 1:p.bs.nks
        counter = 0
        for n = 1:p.bs.nbnd
            if !(n in INDSEMI[k, :])
                counter += 1
                EIGS[k,counter] = p.bs.eigs[k,n]
                PROJ[k,:,counter] = p.proj[k,wan,n]
                
            end
        end
#        if k < 20
#            println("$k EIGS ", EIGS[k,1:5])
#        end
        

#        if k < 2
#            println("ALLEIGS k $k ", p.bs.eigs[k,1:4])
#            println("EIGS    k $k ", EIGS[k,1:4])#
#
#        end
    end

    if shift_energy
        
        println("shifting eigenvalues to match dft atomization energy")
        ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
#        band_en = band_energy(d.bandstruct.eigs[:,nsemi+1:end], d.bandstruct.kweights, nval)
        band_en = band_energy(EIGS, d.bandstruct.kweights, nval)

#        println("d.crys")
#        println(d.crys)

        etypes = types_energy(d.crys)

        
        etot_dft = d.energy
        e_smear = d.energy_smear
        print("e_smear $e_smear")
 #       println("etot_dft $etot_dft etotal_atoms $etotal_atoms etypes $etypes")

#        println("atomization_energy ", etot_dft - etotal_atoms)

#        atomization_energy = etot_dft - etotal_atoms
        atomization_energy = etot_dft - etotal_atoms - etypes  - e_smear
#        println("subtract smearing energy")

        println("atomization_energy $atomization_energy")

        band_en = band_en 
        shift = (atomization_energy - band_en  )/nval

#        println("shift $shift")
#        println("band_en $band_en etot_dft $etot_dft atomization_energy $atomization_energy shift $shift nval $nval")
        
#        p.bs.eigs = p.bs.eigs .+ shift
        EIGS = EIGS .+ shift
        

        if !(ismissing(energy_froz))
            println("energy_froz before $energy_froz")
#            energy_froz = energy_froz + shift
            energy_froz = efermi_dft + shift
            

        end
        println("energy_froz after $energy_froz")
        println("sum shifted EIGS ", sum(EIGS))
    else
        println("no shift: match dft eigenvals")
    end
    
    #testing
#    nsemi = 0
 #   nwan = 2
 #   wan = [1,2]
 #   semicore = []
 #   nfroz = 0
    
    #nsemi = 0
    #nwan = 5
    #wan = [1,2,3,4,5]
    #semicore = []
    #nfroz = 0
    #---
    
#    println("semicore ", semicore, " nsemi ", nsemi)
#    println("wan ", wan, " nwan ", nwan)


    P = zeros(Complex{Float64}, NBND, NBND)  #p.bs.nbnd-nsemi,p.bs.nbnd-nsemi)

    Pmat = zeros(Float64, p.bs.nks, NBND)
    Nmat = zeros(Float64, p.bs.nks, nwan)

    Btilde = zeros(Complex{Float64}, nwan, NBND)
    B = zeros(Complex{Float64}, nwan, NBND)

    #    Pa = zeros(Complex{Float64}, nfroz, nfroz)
    #    Ba = zeros(Complex{Float64}, nwan, nfroz)
    #    Btildea = zeros(Complex{Float64}, nwan, nfroz)
    
    ham_dft = zeros(Float64, NBND)

    ham_k = zeros(Complex{Float64},  p.natwfc-nsemi, p.natwfc-nsemi,p.bs.nks)


    htemp = zeros(Complex{Float64}, nwan, nwan)
    #    htemp = zeros(Complex{Float64}, p.bs.nbnd-nsemi,p.bs.nbnd-nsemi)

#    max_c = minimum(p.bs.eigs[:,end])
    max_c = minimum(EIGS[:,end])
    min_c = max_c - 0.1

#    min_c = max(maximum(p.bs.eigs[:,wan])+.001, min_c) #ensure we don't cut off needed states

    min_c = max(maximum(EIGS[:,1:nwan])+.001, min_c) #ensure we don't cut off needed states
#    min_c = max(maximum(EIGS[:,1:nwan])+.001, min_c) #ensure we don't cut off needed states

    
#    println("min of highest eigenvals: ", max_c, " " , min_c)
    
    function cutoff(num, min_c, max_c)
        if num < min_c
            return 1.0
        elseif num > max_c
            return 0.0
        else
            t = (num - min_c)/(max_c-min_c)
            return 1.0 - 10.0 * t^3 + 15.0 *  t^4  - 6.0 * t^5
        end
    end
#    Palt = zeros(4,4)

    Palt_min = 1.0
    
    badk=0
    for k = 1:p.bs.nks
#    for k = 1:nks
        
        #this avoids breaking symmetry by cutting off a symmetrically equivalent pair/triplet, etc because we include a finite number of bands
#        energy_cutoff = p.bs.eigs[k,end]
        energy_cutoff = EIGS[k,end]
        max_ind = 0
#        for i=max(1,p.bs.nbnd-5):p.bs.nbnd
        for i=max(1,NBND-5):NBND
#            if p.bs.eigs[k,i] < energy_cutoff - 1e-3
            if EIGS[k,i] < energy_cutoff - 1e-3
                max_ind = i
            end
        end
#        if k == 1
#            println("wan $wan nsemi $nsemi max_ind $max_ind NBND $NBND")
#            println("EIGS ", EIGS[k,:])
#        end
        
#        B[:,1:max_ind-nsemi] = p.proj[k,wan, 1+nsemi:max_ind]
        B[:,1:max_ind] = PROJ[k,:, 1:max_ind]

        #smoothing
        for i in 1:max_ind
#            B[:,i] *= cutoff(p.bs.eigs[k,i+nsemi], min_c, max_c)
            B[:,i] *= cutoff(EIGS[k,i], min_c, max_c)
        end
        
        P[:,:] .= 0.0
#        P[1:max_ind-nsemi,1:max_ind-nsemi] = B[:,1:max_ind-nsemi]'*B[:,1:max_ind-nsemi]
        P[1:max_ind,1:max_ind] = B[:,1:max_ind]'*B[:,1:max_ind]
        

        Palt = B[:,1:max_ind]*B[:,1:max_ind]'


        
        for i in 1:nwan
            if real(Palt[i,i]) < Palt_min
                Palt_min = real(Palt[i,i])
            end
            #            println("Palt $i ", real(Palt[i,i]))
            if real(Palt[i,i]) < 0.80
                println("Warning, difficult to project atomic wfc ", i," ",  real(Palt[i,i]), " k ", k, " tr(P) = " , tr(real(P)))
            end
            if real(Palt[i,i]) < 0.45
                badk += 1
                println("bad k ", badk)
            end
            if real(Palt[i,i]) < 0.35
                badk += 1
                println("very bad k ", badk)
            end
            
        end
        
        #        println("P ", real(diag(P)))
        #        println("EIG ", (p.bs.eigs[k,1+nsemi:end]))

        val, vect = eigen( 0.5 * (P[1:max_ind,1:max_ind] + P[1:max_ind,1:max_ind]')     ) # 
#        if k == 1; println("P VAL ", val); end
        good_proj = (max_ind - nwan + 1 ) : max_ind
#        if k == 1; println("good proj", good_proj); end

        Btilde[:,1:max_ind] = vect[:,good_proj]'

        #        for (i, g) = enumerate(good_proj)
        #            Btilde[:,i] =  Btilde[:,i] ./ val[g]^0.5
        #        end
        
#        htemp[:,:] = Btilde[:,1:max_ind-nsemi] * Diagonal(p.bs.eigs[k,nsemi+1:max_ind])  * Btilde[:,1:max_ind-nsemi]'
#        htemp[:,:] = Btilde[:,1:max_ind-nsemi] * Diagonal(EIGS[k,1:max_ind])  * Btilde[:,1:max_ind-nsemi]'
        htemp[:,:] = Btilde[:,1:max_ind] * Diagonal(EIGS[k,1:max_ind])  * Btilde[:,1:max_ind]'

        neweigs, newvect = eigen((htemp+htemp')/2.0)

#        if k == 1
#            println( "NEW EIGS !!!!!!!!!!!!!!! ", neweigs)
#        end
        #        println("neweig ", neweigs)
        Pmat[k, :] = real(diag(P))
        Nmat[k, :] = neweigs
        
#        ham_k[:,:,k] = B[:,1:max_ind-nsemi] * Btilde[:,1:max_ind-nsemi]' * htemp * Btilde[:,1:max_ind-nsemi] *  B[:,1:max_ind-nsemi]'        
        ham_k[:,:,k] = B[:,1:max_ind] * Btilde[:,1:max_ind]' * htemp * Btilde[:,1:max_ind] *  B[:,1:max_ind]'        
        
        
        ham_k[:,:,k] = (ham_k[:,:,k]  + ham_k[:,:,k]')/2.0

#        if k == 1
#            println( "hamk ", eigvals(ham_k[:,:,k]))
#        end


#        if k < 10
#            println("k $k")
#            println(EIGS[k, 1:10])
#            println("calc xxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#            println(eigvals(ham_k[:,:,k]))
#            println()
#        end            

    end

    if badk / p.bs.nks > 0.15
        println("warning lots of bad kpoints ", badk, " of ", p.bs.nks)
        warn_badk=true
    else
        warn_badk=false
    end

    println("Min Palt: $Palt_min")

    if !(ismissing(energy_froz))
        energy_froz2 = energy_froz+0.25
        println("energy_froz: $energy_froz , $energy_froz2")
        for k = 1:p.bs.nks
            val_tbt, vect = eigen(ham_k[:,:,k] )
            val_tb = real(val_tbt)
            val_tb_new = deepcopy(val_tb)
#            val_pw = p.bs.eigs[k,nsemi+1:nsemi+nwan]
            val_pw = EIGS[k,1:nwan]
            
            order = Dict()
            
            score_mat = zeros(nwan, nwan)
            

            for n2 = 1:nwan
                cd_vect = real(vect[:,n2].*conj(vect[:,n2]))

                dist_min = 1000.0
                nmin = 0

                for n1 = 1:nwan
                    
#                    t = p.proj[k,wan, nsemi + n1] 
                    t = PROJ[k,:,n1]

                    cd_dft = real(t .* conj(t))
                
                    dist_en = (val_tb[n2] - val_pw[n1]).^2 * 10.0
                    dist_cd = sum((cd_dft - cd_vect).^2)
                    dist = dist_en + dist_cd + 0.1*(n1 - n2)^2
                    score_mat[n1,n2] = dist
#                    if dist < dist_min
#                        dist_min = dist
#                        nmin = n1
#                    end
                end
            end
            for n = 1:nwan
                s = sortperm(score_mat[n,:])
                for ss in s
                    if !(ss in keys(order))
                        order[ss] = n
                        break
                    end
                end
            end
#            if k < 10
#                for n = 1:nwan
#                    println("ORDER $n ", order[n])
#                end
#            end

            for n in 1:nwan
                if val_pw[order[n]] < energy_froz
                    val_tb_new[n] = val_pw[order[n]]
                elseif val_pw[order[n]] < energy_froz2
                    x=cutoff(val_tb[n], energy_froz, energy_froz2)
                    val_tb_new[n] = val_pw[order[n]] * x + val_tb[n]*(1.0-x)
                end
                
            end

#resym
            for n = 1:nwan-1
                c= [n]
                for n2 = n+1:nwan
                    if abs(val_tb_new[n] - val_tb_new[n2]) < 1e-4
                        push!(c, n2)
                    end
                end
                if length(c) > 1
                    t = sum(val_tb_new[c]) / length(c)
                    val_tb_new[c] .= t
                end
            end
                
            ham_k[:,:,k] = vect*Diagonal(val_tb_new)*vect'
            ham_k[:,:,k] = (ham_k[:,:,k]  + ham_k[:,:,k]')/2.0
            val_tb_new, vect = eigen(ham_k[:,:,k] )
            
            

        end
    end
    
#    for k = 1:10
#
#        println("$k ", p.bs.kpts[k,:])
#        val_tb_new, vect = eigen(ham_k[:,:,k] )
#        println("val new ", val_tb_new)
#        println("val dft ", EIGS[k,1:nwan])
#    end

    

    if nfroz >= 1 && (ismissing(energy_froz))
        println("nfroz: ", nfroz)
        for k = 1:p.bs.nks
            val, vect = eigen(ham_k[:,:,k] )
            val[1:nfroz] = EIGS[k,1:nfroz]
            if (abs(EIGS[k, nfroz+1] - EIGS[k, nfroz])< 1e-5) && nwan >= nfroz+1
                val[nfroz+1] = EIGS[k,nfroz+1]
                if (abs(EIGS[k, nfroz+2] - EIGS[k, nfroz])< 1e-5) && nwan >= nfroz+2
                    val[nfroz+2] = EIGS[k,nfroz+2]
                end
            end
            
                


            ham_k[:,:,k] = vect*Diagonal(val)*vect'
            ham_k[:,:,k] = (ham_k[:,:,k]  + ham_k[:,:,k]')/2.0

        end
    end

   
    #reshift!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    VAL = zeros(Float64, p.bs.nks, nwan)
    VECT = zeros(Complex{Float64}, p.bs.nks, nwan, nwan)
    for k = 1:p.bs.nks
        val, vect = eigen(ham_k[:,:,k])
        VAL[k,:] = val
        VECT[k,:,:] = vect
    end
        
    if shift_energy
        println("reshift")
#        println("size VAL ", size(VAL), " size kweights ", size(d.bandstruct.kweights))
        band_en = band_energy(VAL, d.bandstruct.kweights, nval)

        println("band_en_old " , band_en)

        shift = (atomization_energy - band_en)/nval
        VAL = VAL .+ shift
        for k = 1:p.bs.nks
            
            val = VAL[k,:]
            vect = VECT[k,:,:]
            ham_k[:,:,k] = vect*Diagonal(val)*vect'
            
        end
        println("done reshift")

        for k = 1:p.bs.nks
            val, vect = eigen((ham_k[:,:,k] + ham_k[:,:,k]')/2.0)
            VAL[k,:] = val
            VECT[k,:,:] = vect
        end

        band_en = band_energy(VAL, d.bandstruct.kweights, nval)

        println("band_en_new " , band_en, " " , band_en+etypes)


    end
    
    
#    writedlm("new.csv", real(ham_k))
    return ham_k, EIGS, Pmat, Nmat, VAL, warn_badk




end
    

function get_ham_r_slow(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 3}; nonorth=true, K=missing)
"""
DOESN'T USE FFTW. MOSTLY FOR DEBUGGING / missing fft libraries

Do Fourier transform k->R. grid is the size of the k-space grid. 
Optionally you can include the k space grid explictly, but normally 
K is already inside p. K must match ham_k
"""

    
    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)
    
#    grid = d.bandstruct.kgrid

    println("grid ", grid)
    println("ham_k size " , size(ham_k))
    println("nonorth: ", nonorth)
    nwan = size(ham_k)[2]

    if ismissing(K)
        K = p.bs.kpts
        nks = size(K)[1]
    else
        nks = size(K)[1]
    end
    
    if nonorth 
        #nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ham_kS = zeros(Complex{Float64},  nwan, nwan,nks)
        S = zeros(Complex{Float64}, nwan, nwan)
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        Sk2 = zeros(Complex{Float64},  nwan, nwan,nks)
        for k in 1:nks
            S[:,:] = p.overlaps[k, wan,wan]
            S = (S+S')/2.0
            Sk[:,:,k] =  S
        end
        for i=1:nwan #normalize overlaps
            sumi = sum(Sk[i,i,:])
            Sk[i,i,:] = Sk[i,i,:] / real(sumi) * Float64(nks)
        end
        for k in 1:nks
            S[:,:] = Sk[:,:,k]
            S = (S + S')/2.0
            Sk2[:,:,k] = sqrt(S)
            
            ham_kS[:,:,k] = Sk2[:,:,k]*ham_k[:,:,k]*Sk2[:,:,k]
            ham_kS[:,:,k] = (ham_kS[:,:,k]  + ham_kS[:,:,k]')/2.0
        
        end
        #end nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
        
    
    r=zeros(Float64,3)
    twopi_i = 1.0im*2.0*pi

    grid2 = [0,0,0]
    for i = 1:3
        if grid[i]%2 == 0
            grid2[i] = Int(grid[i]/2)
        else        
            grid2[i] = Int((grid[i]-1)/2)
        end
    end
    
    r_dict = Dict()


    ngrid2 = (grid2[1]*2+1)*(grid2[2]*2+1)*(grid2[3]*2+1)
    println("new real space grid: " , -grid2, " to " , grid2, " " , ngrid2)

    ind_arr = zeros(ngrid2,3)
    ham_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)

    if nonorth
        S_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)
    end
    
#    return 0

    sym = zero(Float64)
    exp_val = zero(Complex{Float64})
    c=0

    sym_mat = zeros(Float64, ngrid2)

    for r1 = -grid2[1]:(grid2[1])

        if (r1 + grid[1] == grid2[1]) ||  (r1 - grid[1] == -grid2[1])
            sym1 = 2.0
        else
            sym1 = 1.0
        end
            
        for r2 = -grid2[2]:(grid2[2])
            if (r2 + grid[2] == grid2[2]) ||  (r2 - grid[2] == -grid2[2])
                sym2 = 2.0
            else
                sym2 = 1.0
            end

            
            for r3 = -grid2[3]:(grid2[3])

                if (r3 + grid[3] == grid2[3]) ||  (r3 - grid[3] == -grid2[3])
                    sym3 = 2.0
                else
                    sym3 = 1.0
                end


#                println("rs $r1 $r2 $r3 $sym1 $sym2 $sym3 $sym")
                
                c+=1 
                sym_mat[c] = sym1 * sym2 * sym3
                r[:] = [r1 r2 r3]
                ind_arr[c,:] = r[:]
                r_dict[ind_arr[c,:]] = c
            end
        end
    end


    if nonorth
        for k in 1:nks
            for c = 1:ngrid2            
                exp_val = exp(twopi_i * (ind_arr[c,:]'*K[k,:])) / sym_mat[c]
                ham_r[:,:, c ] .+= ham_kS[:,:,k] .* exp_val
                S_r[:,:, c] .+=  Sk[:,:,k] .* exp_val

            end
        end
            
    else
        for k in 1:nks
            for c = 1:ngrid2            
                exp_val = exp(twopi_i * (ind_arr[c,:]'*K[k,:])) / sym_mat[c]
                ham_r[ :,:, c] .+= ham_k[:,:,k] .* exp_val
            end
        end
    end



    #turn ham_r into tight-binding object
    if nonorth
        tb = make_tb(ham_r/ Float64(nks), ind_arr, r_dict, S_r / Float64(nks))
    else
        tb = make_tb(ham_r/ Float64(nks), ind_arr, r_dict)
    end

#    renormalize_tb(d, tb)

    nelec=d.bandstruct.nelec - nsemi * 2.0

    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
    dftenergy=d.energy - etotal_atoms

    model = make_tb_crys(tb, d.crys, nelec, dftenergy)

    
    return model
    
end    
    

function prepare_ham_k(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 3}; nonorth=true, K=missing, localized_factor = 0.05, screening=1.0)

    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)
    
    println("grid ", grid)
    println("ham_k size " , size(ham_k))
    println("nonorth: ", nonorth)
    nwan = size(ham_k)[2]
    
    if ismissing(K)
        K = p.bs.kpts
        nks = size(K)[1]
    else
        nks = size(K)[1]
    end
    
    #fix mysterious sign error K-space
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
    
    OVERLAPS = deepcopy(p.overlaps)
    println("wan ", wan)
    println("so ", size(OVERLAPS))
    OVERLAPS = OVERLAPS[:,wan, wan]
    println("so2 ", size(OVERLAPS))
    for n = 1:nwan
        a1,t1,orb = ind2orb[n]
        if orb == :px || orb == :py
            println("fix mysterious sign issue K-space $n, $a1, $t1, $orb")
            ham_k[n,:,:] = -1.0 * ham_k[n,:,:]
            ham_k[:,n,:] = -1.0 * ham_k[:,n,:]

            if nonorth
                OVERLAPS[:,n,:] = -1.0*OVERLAPS[:,n,:]
                OVERLAPS[:,:,n] = -1.0*OVERLAPS[:,:,n]
                
            end
        elseif orb == :dxz || orb == :dyz
            println("fix mysterious sign issue K-space D ORBITALS $n, $a1, $t1, $orb")
            ham_k[n,:,:] = -1.0 * ham_k[n,:,:]
            ham_k[:,n,:] = -1.0 * ham_k[:,n,:]

            if nonorth
                OVERLAPS[:,n,:] = -1.0*OVERLAPS[:,n,:]
                OVERLAPS[:,:,n] = -1.0*OVERLAPS[:,:,n]
        
            end

        end
    end




    if nonorth 
        #nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ham_kS = zeros(Complex{Float64},  nwan, nwan,nks)
        S = zeros(Complex{Float64}, nwan, nwan)
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        Sk2 = zeros(Complex{Float64},  nwan, nwan,nks)
            

        
        for k in 1:nks
#            S[:,:] = p.overlaps[k, wan,wan]
            S[:,:] = OVERLAPS[k, :,:]

            S = (S+S')/2.0           #+ I(nwan) * localized_factor
            Sk[:,:,k] =  S
        end

#        println("fourteen")
#        println(Sk[:,:,14])

        wtot = sum(p.bs.kweights)
        SUMI = zeros(nwan)
        for i=1:nwan #normalize overlaps
            SUMI[i] = sum(Sk[i,i,:] .* p.bs.kweights) / wtot
        end
        println("old SUMI ", SUMI)
        SUMI = symm_by_orbitals(d.crys, SUMI) #this handles the fact the k-grid can break orbital symmetry
        println("new SUMI ", SUMI)
        for i=1:nwan #normalize overlaps
            Sk[i,i,:] = Sk[i,i,:] / SUMI[i]
        end


#        for i=1:nwan #normalize overlaps
 #           sumi = sum(Sk[i,i,:] .* p.bs.kweights)
#            println("normalize $i ", sumi/sum(p.bs.kweights) )
#            Sk[i,i,:] = Sk[i,i,:] / real(sumi) * sum(p.bs.kweights)
#        end

#        println("fourteen2")
#        println(Sk[:,:,14])

        for k in 1:nks
            S = Sk[:,:,k]
            S = (S + S')/2.0
            Sk[:,:,k] = S
#            vals, vects = eigen(S)
#            Sprime = vects*diagm(vals .+ localized_factor)*vects'
#            Sk[:,:,k] = Sprime
        end            

#        Sk = Sk ./ (1.0 + localized_factor)

        
#        for i=1:nwan #normalize overlaps
#            sumi = sum(Sk[i,i,:])
#            Sk[i,i,:] = Sk[i,i,:] / real(sumi) * Float64(nks)
#        end

#        for i=1:nwan #normalize overlaps
#            sumi = sum(Sk[i,i,:] .* p.bs.kweights)
#            Sk[i,i,:] = Sk[i,i,:] / real(sumi) * sum(p.bs.kweights)
#        end

        
        if localized_factor > 1e-5
            for k = 1:nks
                Sk[:,:,k] = Sk[:,:,k] * (1.0 - localized_factor) + I(nwan) * localized_factor 
            end
        end

#        println("fourteen3")
#        println(Sk[:,:,14])
        
        
        for k in 1:nks
            S[:,:] = Sk[:,:,k]
            S = (S+S')/2.0
            Sk2[:,:,k] = sqrt(S)
            
            ham_kS[:,:,k] = Sk2[:,:,k]*ham_k[:,:,k]*Sk2[:,:,k]
            ham_kS[:,:,k] = (ham_kS[:,:,k]  + ham_kS[:,:,k]')/2.0
        
        end

        tbk = make_tb_k(ham_kS, K, d.bandstruct.kweights, Sk, grid=grid, nonorth=true)
        tbck = make_tb_crys_kspace(tbk, d.crys, nval, d.energy - etotal_atoms, scf=false, screening=screening)

        return tbck

        #end nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        for i = 1:nks
            Sk[i,:,:] = I(nwan)
        end

        tbk = make_tb_k(ham_k, K, d.bandstruct.kweights, Sk, grid=grid, nonorth=false)
        tbck = make_tb_crys_kspace(tbk, d.crys, nval, d.energy - etotal_atoms, scf=false)
        
        return tbck

    end

    ###########


end



#function get_ham_r(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 3}; nonorth=true, K=missing, localized_factor = 0.05)
function get_ham_r(tbck::tb_crys_kspace)
"""
Do Fourier transform k->R. grid is the size of the k-space grid. 
Optionally you can include the k space grid explictly, but normally 
K is already inside p. K must match ham_k
"""

    nonorth = tbck.tb.nonorth
    ham_k = tbck.tb.Hk
    Sk = tbck.tb.Sk
    K = tbck.tb.K
    grid = tbck.tb.grid
    nwan = tbck.tb.nwan

    if nonorth
        ham_r, S_r, r_dict, ind_arr = myfft(tbck.crys, nonorth, grid, K,ham_k, Sk)
    else
        ham_r, r_dict, ind_arr =      myfft(tbck.crys, nonorth, grid, K,ham_k, missing)
    end             

    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbck.crys)

#=
    #fix mysterious sign issue


    for n = 1:nwan
        a1,t1,orb = ind2orb[n]
        if orb == :px || orb == :py
            println("fix mysterious sign issue R-space $n, $a1, $t1, $orb")
            ham_r[n,:,:] = -1.0 * ham_r[n,:,:]
            ham_r[:,n,:] = -1.0 * ham_r[:,n,:]

            if nonorth
                S_r[:,n,:] = -1.0*S_r[:,n,:]
                S_r[n,:,:] = -1.0*S_r[n,:,:]
                
            end
        elseif orb == :dxz || orb == :dyz
            println("fix mysterious sign issue R-space D ORBITALS $n, $a1, $t1, $orb")
            ham_r[n,:,:] = -1.0 * ham_r[n,:,:]
            ham_r[:,n,:] = -1.0 * ham_r[:,n,:]

            if nonorth
                S_r[:,n,:] = -1.0*S_r[:,n,:]
                S_r[n,:,:] = -1.0*S_r[n,:,:]
        
            end

        end
    end
=#
    #turn ham_r into tight-binding object
    if nonorth
        tb = make_tb(ham_r, ind_arr, r_dict, S_r )        
    else
        S_r = zeros(size(ham_r))

        println(size(S_r))

        r = r_dict[[0,0,0]]
        for i =  1:nwan
            S_r[i,i,r] = 1.0
        end
        tb = make_tb(ham_r, ind_arr, r_dict, S_r)        
    end

#    renormalize_tb(d, tb)
#    nelec=.bandstruct.nelec ##- nsemi * 2.0
    nelec = nval

    dftenergy=tbck.dftenergy
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbck.crys)

    model = make_tb_crys(tb, tbck.crys, nelec, dftenergy )
    
    return model
    
end    
    








##function run_proj(inputstr, outputstr, nprocs=1, directory="./")
##"""
##run command
##"""
##    directory=rstrip(directory, '/')
##    
##    #get commandline
##    c_dict = make_commands(nprocs)
##    qe = c_dict["qe"]
##    
##
##    command = `$qe $directory/$inputstr`
##    println("actual command")
##    println(command)
##    s=""
##    try
##        println("Running DFT")
##
##        s = read(command, String)
##        
##        println("Writing output")
##        f = open(directory*"/"*outputstr, "w")
##        write(f, s)
##        close(f)
##        
##        println("Ran DFT")
##        return 0
##    catch
##        println("Failed to run DFT")
##        println(s)
##        return -1
##    end
##        
##end
##


end
