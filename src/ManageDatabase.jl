module ManageDatabase

#using FileIO
#using JLD2

using ..CrystalMod:crystal
using ..CalcTB:coefs

using ..TightlyBound:DATSDIR
using ..CalcTB:read_coefs

defaultdatdir = DATSDIR

database_list = Set()

database_cached = Dict()
database_cached["SCF"] = true
database_cached["scf"] = true


function prepare_database(c::crystal)
    println("prepare c ", c.types)
    prepare_database(c.types)
end

function prepare_database(at_list)

    println("prepare at ", at_list)
    println("database_list ", database_list)
    s = Set(at_list)
    for s1 in s
        for s2 in s
            if ! ( Set([s1,s2]) in database_list)
                add_to_database(Set([s1,s2]))
            end
        end
    end
end


function add_to_database(s::Set)#

    println("add_to_database $s")

    at_arr = collect(s)
    if length(s) == 1
        a1 = at_arr[1]
        if !haskey(database_cached , (a1, a1))
            f = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_el.$a1.xml"
            if isfile(f) || isfile(f*".gz")
                try
                    dat = read_coefs(f)
                    database_cached[(a1, a1)] = dat
                    println("added to cache ", (a1, a1))
                catch
                    println("warning - error loading $f")
                end
            else
                println("warning - no file for database $a1 ")
            end
        end

        if !haskey(database_cached , (a1, a1, a1))
            f = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_el.$a1.xml"
            if isfile(f) || isfile(f*".gz")
                try
#                    jldopen(f)                    
                    dat = read_coefs(f)
                    database_cached[(a1, a1, a1)] = dat
                    println("added to cache ", (a1, a1, a1))
                catch
                    println("warning - error loading $f")
                end
            else
                println("warning, no file for 3bdy database $a1")
            end
        end
           
 
    elseif length(s) == 2

        a1 = at_arr[1]
        a2 = at_arr[2]
        if !haskey(database_cached , (a1, a2))
            
            fab = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a1.$a2.xml"
            fba = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a2.$a1.xml"
            
            if isfile(fab) || isfile(fab*".gz")
                f = fab
            elseif isfile(fba) || isfile(fba*".gz")
                f = fba
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

            fab = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a1.$a2.xml"
            fba = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a2.$a1.xml"
            
            if isfile(fab) || isfile(fab*".gz")
                f = fab
            elseif isfile(fba) || isfile(fba*".gz")
                f = fba
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
