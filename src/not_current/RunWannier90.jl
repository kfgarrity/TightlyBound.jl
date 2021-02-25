#Not currently used for anything

#=

###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module W90
"""
Scripts to run Wannier90
"""
#export load_xml
using LinearAlgebra

using XMLDict
using ..DFToutMod
using ..CrystalMod:crystal
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..QE:loadXML, loadXML_bs

include("Atomdata.jl")
include("Commands.jl")

function write_to_file(str, filename, directory="./")

    f=open("$directory/$filename", "w")
    write(f, str)
    close(f)

end

#function test(dft::dftout)
#    return dft
#end

function qe2w90_workf(dft::dftout; seedname="w90", directory="./", nprocs=1)
"""
All steps to run wannier90 from converged QE SCF calculation

Steps:

1) open_grid scf calc 
2) wannier90 preprocessing
3) pw2wan
4) wannier90
"""

    prefix=dft.prefix
    outdir=dft.outdir

    
    println()
    println("qe2w90_workflow-------------------------------------------------------------------")
    println()
    println("Step #1 og------------------------------------------------------------------------")
    println()
    
    #1 og
    og= makeOG(prefix, outdir)
    write_to_file(og, "og.in", directory)
    run_og("og.in", directory=directory, nprocs=nprocs)
    
    newprefix=prefix*"_open"

    dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs

    kgrid = copy(dft.bandstruct.kgrid)
    dft.bandstruct = dftopen_bs
    dft.bandstruct.kgrid = kgrid
    
    println()
    println("Step #2 w90 pp--------------------------------------------------------------------")
    println()    
    #2 w90 pp

    println("PREFIX: $prefix")
    println("OUTDIR: $outdir")
    
    win = makeWIN(dft) #needs the unfolded k-points and band structure !!!!!

    write_to_file(win, "$seedname.win", directory)

    run_wannier90x(seedname, preprocess=true, directory=directory) #preprocess
    
    println()
    println("Step #3 pw2wan--------------------------------------------------------------------")
    println()
    #3 pw2wan
    pw2wan = makePW2WAN(newprefix, outdir, seedname)

    write_to_file(pw2wan, "qe.pw2wan", directory)
    
    run_pw2wannier90x("qe.pw2wan"; directory=directory, nprocs=nprocs) 

    println()
    println("Step #4 w90 main run--------------------------------------------------------------")
    println()
    #4 w90 for real
    run_wannier90x(seedname, preprocess=false, directory=directory) #normal mode

    println()    
    println("Done!!!!!!!!!!!!!!!!--------------------------------------------------------------")
    println()

    
end
                      


function run_wannier90x(seedname; preprocess = False, directory="./")
"""
run wannier90.x command
"""
    c_dict = make_commands(1)
    w90 = c_dict["wannier90"]


    if preprocess == true
        command=`$w90 -pp $directory/$seedname `
    else
        command=`$w90  $directory/$seedname `
    end

    println("wannier90 command")
    println(command)
    
    try
        println("Running wannier90.x")

        s = read(command, String)
        println("wannier90.x stdout:")
        println(s)
        println()
        
#        f = open(directory*"/"*seedname*".wout", "w")
#        write(f, s)
#        close(f)
        if preprocess
            println("Ran wannier90.x -pp")
        else
            println("Ran wannier90.x")
        end            
        return 0
    catch
        println("Failed to run wannier90.x")
        return -1
    end
        
end

function run_pw2wannier90x(pw2wanfile; directory="./", nprocs=1)
"""
run pw2wannier90.x command
"""

    c_dict = make_commands(nprocs)
    qe = c_dict["pw2wan"]
    command = `$qe $directory/$pw2wanfile `
    println("pw2wan command")
    println(command)
    println()
    
    try
        println("Running pw2wannier90.x")

        s = read(command, String)
        f = open(directory*"/"*pw2wanfile*".out", "w")
        write(f, s)
        close(f)
        
        println("Ran pw2wannier90.x")
        return 0
    catch
        println("Failed to run pw2wannier90.x")
        return -1
    end
        
end

function run_og(filename="og.in";  directory="./", nprocs=1)
"""
run pw2wannier90.x command
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

function makeOG(prefix, tmpdir )
"""
makes the pw2wan file
"""
    template_file=open("template_inputs/template_og.in")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    
    return temp
    
end


function makePW2WAN(prefix, tmpdir, seedname; spinpol=false)
"""
makes the pw2wan file
"""
    template_file=open("template_inputs/template.pw2wan")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    temp = replace(temp, "SEEDNAME" => "$tmpdir/$seedname")
    other = ""
    if spinpol
        println("makepw2wan no spinpol yet")
    end
    temp = replace(temp, "JULIAOTHER" => other)            
    
    return temp
    
end
    

function makeWIN(dft::dftout)
"""
Make inputfile for W90 calculation
"""
    convert_ryd_ev = 13.605693122

    
    template_file=open("template_inputs/template.win")
    temp = read(template_file, String)
    close(template_file)
    
    crys=dft.crys
    

    temp = replace(temp, "JULIANAT" => crys.nat)

    nbandsemi = 0
    nwan = 0        

    for t in crys.types
        nbandsemi += atoms[t].nsemicore
        nwan += atoms[t].nwan

    end

   
    nbandsemi=convert(Int64, nbandsemi/2)
    nwan=convert(Int64, nwan/2)
    
    nbndtot = dft.bandstruct.nbnd - nbandsemi

    if nbndtot < nwan
        exit("error, dft calculation needs more bands to accomadate all wannier functions", nbandsemi,nwan,dft.bandstruct.nbnd)
    end
    
    settypes = []
    for t in crys.types
        if !(t in settypes)
            push!(settypes, t)
        end
    end

    proj = ""
    for t in settypes
        proj *= t*":" 
        for o in atoms[t].orbitals
            proj *= string(o)*","
        end
        proj = strip(proj, ',')
        proj=proj*"\n"
    end        
    proj = strip(proj, '\n')
        
    
    temp = replace(temp, "JULIANUMWANN" => nwan)
    temp = replace(temp, "JULIANUMBANDS" => nbndtot)

    if nbandsemi > 0 
        temp = replace(temp, "JULIAEXCLUDE" =>  "1-$nbandsemi")
    end
        
    temp = replace(temp,"JULIAPROJECTIONS" => proj)



    
    eigs = vec(dft.bandstruct.eigs)*convert_ryd_ev
    fermi = dft.bandstruct.efermi*convert_ryd_ev

    ind=findall(x->x > fermi, eigs)

    vol = abs(det(crys.A))
    
    if (crys.nat == 1 &&  vol > 10*10*9.99) || (crys.nat == 2 &&  vol > 20*10*9.99)
        println("Frozen +0.5 low density system")
        frozmax = minimum(eigs[ind]) + 0.5        #to deal with spread out systems with highly localized bands. we don't want continuum states frozen, that is bad.
    else
        frozmax = minimum(eigs[ind]) + 1.5        
    end

    # there cannot ever be more frozen bands than WFs
    if nwan+nbandsemi < dft.bandstruct.nbnd
        frozmax_max = minimum(dft.bandstruct.eigs[:,nbandsemi+nwan+1])-0.1
        frozmax = min(frozmax_max, frozmax)
    end
    
    temp = replace(temp, "JULIAFROZMAX" => frozmax)

    
    
    temp = replace(temp, "JULIACELL\n" => arr2str(crys.A))

    t=""
    for i = 1:crys.nat
        t=t*crys.types[i]*" "
        t=t*arr2str(crys.coords[i,:])
        if i != crys.nat
            t=t*"\n"
        end
    end
    temp = replace(temp, "JULIAATOMS" => t)

    
    temp = replace(temp, "JULIAKGRID\n" => arr2str(dft.bandstruct.kgrid))
    temp = replace(temp, "JULIAKPOINTS\n" => arr2str(dft.bandstruct.kpts))  

    return temp
    
                   
    
end

    
end #end W90

=#
