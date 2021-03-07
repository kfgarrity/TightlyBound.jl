"""
    module ManageDatabase

Module for reading `coefs` from files and making database as needed for calculations
"""
module ManageDatabase

#using FileIO
#using JLD2

using ..CrystalMod:crystal
using ..CalcTB:coefs

using ..TightlyBound:DATSDIR1
using ..TightlyBound:DATSDIR2
using ..CalcTB:read_coefs

datdir1 = DATSDIR1 #primary directory
datdir2 = DATSDIR2 #backup directory

database_list = Set()

database_cached = Dict()
database_cached["SCF"] = true
database_cached["scf"] = true

"""
    function prepare_database(c::crystal)

Get ready database of precalculated `coefs` for `crystal`
"""
function prepare_database(c::crystal)
#    println("prepare c ", c.types)
    prepare_database(c.types)
end

"""
    function prepare_database(at_list)
"""
function prepare_database(at_list)
    
    println("prepare atoms ", at_list)
    at_list = Symbol.(at_list)
    #    println("database_list ", database_list)
    s = Set(at_list)
    for s1 in s
        for s2 in s
            if ! ( Set([s1,s2]) in database_list)
                add_to_database(Set([s1,s2]))
            end
        end
    end
end

"""
    function add_to_database(s::Set)

Load elements or twobody terms from precalcuated `coefs` from files.
"""
function add_to_database(s::Set)#

    println("add_to_database $s")

    at_arr = collect(s)
    if length(s) == 1
        a1 = Symbol(at_arr[1])
        if !haskey(database_cached , (a1, a1))
#            f = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_el.$a1.xml"
            f =  "$datdir1/coef.el.2bdy.$a1.xml.gz"
            f2 = "$datdir2/coef.el.2bdy.$a1.xml.gz"
            if isfile(f) || isfile(f*".gz")
                try
                    dat = read_coefs(f)
                    database_cached[(a1, a1)] = dat
                    println("added to cache ", (a1, a1))
                catch
                    println("warning - error loading $f")
                end
            elseif isfile(f2) || isfile(f2*".gz")
                try
                    dat = read_coefs(f2)
                    database_cached[(a1, a1)] = dat
                    println("added to cache ", (a1, a1))
                catch
                    println("warning - error loading $f2")
                end
            else
                println("warning - no file for database $a1 ")
            end
        end

        if !haskey(database_cached , (a1, a1, a1))
#            f = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_el.$a1.xml"
            f =  "$datdir1/coef.el.3bdy.$a1.xml.gz"
            f2 = "$datdir2/coef.el.3bdy.$a1.xml.gz"
            if isfile(f) || isfile(f*".gz")
                try
#                    jldopen(f)                    
                    dat = read_coefs(f)
                    database_cached[(a1, a1, a1)] = dat
                    println("added to cache ", (a1, a1, a1))
                catch
                    println("warning - error loading $f")
                end
            elseif isfile(f2) || isfile(f2*".gz")
                try
#                    jldopen(f)                    
                    dat = read_coefs(f2)
                    database_cached[(a1, a1, a1)] = dat
                    println("added to cache ", (a1, a1, a1))
                catch
                    println("warning - error loading $f2")
                end
            else
                println("warning, no file for 3bdy database $a1")
            end
        end
           
 
    elseif length(s) == 2

        a1 = Symbol(at_arr[1])
        a2 = Symbol(at_arr[2])
        if !haskey(database_cached , (a1, a2))
            
#            fab = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a1.$a2.xml"
#            fba = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a2.$a1.xml"

            fab =  "$datdir1/coef.el.2bdy.$a1.$a2.xml.gz"
            fba =  "$datdir1/coef.el.2bdy.$a2.$a1.xml.gz"

            fab2 = "$datdir2/coef.el.2bdy.$a1.$a2.xml.gz"
            fba2 = "$datdir2/coef.el.2bdy.$a2.$a1.xml.gz"

            
            if isfile(fab) || isfile(fab*".gz")
                f = fab
            elseif isfile(fba) || isfile(fba*".gz")
                f = fba
            elseif isfile(fab2) || isfile(fab2*".gz")
                f = fab2
            elseif isfile(fba2) || isfile(fba2*".gz")
                f = fba2
            else
                f = missing
                println("warning - binary file missing ")
            end

            if !ismissing(f)
                try

                    dat = read_coefs(f)


                    database_cached[(a1, a2)] = dat
                    database_cached[(a2, a1)] = dat
                    println("added to cache ", (a1, a2), " twobody ")
                    
                catch
                    println("warning - error loading binary file $f")
                end
            end
#################

#            fab = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a1.$a2.xml"
#            fba = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a2.$a1.xml"

            fab =  "$datdir1/coef.el.3bdy.$a1.$a2.xml.gz"
            fba =  "$datdir1/coef.el.3bdy.$a2.$a1.xml.gz"

            fab2 = "$datdir2/coef.el.3bdy.$a1.$a2.xml.gz"
            fba2 = "$datdir2/coef.el.3bdy.$a2.$a1.xml.gz"
            
            if isfile(fab) || isfile(fab*".gz")
                f = fab
            elseif isfile(fba) || isfile(fba*".gz")
                f = fba
            elseif isfile(fab2) || isfile(fab2*".gz")
                f = fab2
            elseif isfile(fba2) || isfile(fba2*".gz")
                f = fba2

            else
                f = missing
                println("warning - binary file missing ")
            end

            if !ismissing(f)
                try
#                    jldopen(f)
                    dat = read_coefs(f)

                    database_cached[(a1,a1,a2)] = dat
                    database_cached[(a1,a2,a1)] = dat
                    database_cached[(a2,a1,a1)] = dat
                    
                    database_cached[(a1,a2,a2)] = dat
                    database_cached[(a2,a1,a2)] = dat
                    database_cached[(a2,a2,a1)] = dat
                    println("added to cache ", (a1, a2), " threebody ")
                catch
                    println("warning - error loading binary 3body $f ")
                end
            end
            
        end
    else
        println("warning - not setup to load ternary ", s)
    end

    
    push!(database_list, s)

    
end





end #end module
