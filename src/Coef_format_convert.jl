

function get_data_info_oldversion(at_set, dim)

    # very  minor changes, see differs. This is soley for format conversion from old coefs to new format, which has extra 3body terms.

    n_2body = 5
    n_2body_onsite = 4
    n_2body_S = 7
    n_3body = 5 #differs. is 7 in the main version
    n_3body_same = 5
    n_3body_onsite = 4
    n_3body_onsite_same = 2


    data_info = Dict{Tuple, Array{Int64,1}}()
    orbs = []
    if dim == 2 #2body

        at_list = [i for i in at_set]
#        println(at_list)
        if length(at_list) == 1
            at_list = [at_list[1], at_list[1]]
        end
        sort!(at_list)
#        println(at_list)

        orbs1 = atoms[at_list[1]].orbitals
        orbs2 = atoms[at_list[2]].orbitals

        at1 = at_list[1]
        at2 = at_list[2]

        if at1 == at2
            same_at = true
        else
            same_at = false
        end

#        orbs = []

        #2body part
        function get2bdy(n, symb)
            tot=0
            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
#                    push!(orbs, (o1, o2, symb))
                    if o1 == :s && o2 == :s
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n
                        tot += n
                        
                    elseif (o1 == :s && o2 == :p ) || (o1 == :p && o2 == :s )
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n
                        tot += n
#                        if same_at
#                            data_info[(o2, o1, symb)] = data_info[(o1, o2, symb)]
#                        end
                        
                    elseif (o1 == :p && o2 == :p )
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n*2
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n*2
                        tot += n*2

#                    elseif (o1 == :p && o2 == :p )
#                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n*2
#                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n*2
#                        tot += n*2

                    elseif (o1 == :s && o2 == :d ) || (o1 == :d && o2 == :s )
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n
                        tot += n

                    elseif (o1 == :p && o2 == :d ) || (o1 == :d && o2 == :p )
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n*2
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n*2
                        tot += n*2

                    elseif (o1 == :d && o2 == :d ) 
                        data_info[(at1, o1, at2, o2, symb)] = tot+1:tot+n*3
                        data_info[(at2, o2, at1, o1, symb)] = tot+1:tot+n*3
                        tot += n*3


                        
                    end
                end
            end
            return tot
        end

        totH = get2bdy(n_2body, :H)
        totS = get2bdy(n_2body_S, :S)
#        println("totH $totH totS $totS")
        #onsite part
        function getonsite(atX,orbsX, tot, n)
            for o1 in orbsX
                for o2 in orbsX
                    if (o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d)
                        continue
                    end

#                    push!(orbs, (o1, o2, :O))
                    if o1 == :s && o2 == :s
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n
#                        println("data_info" , (atX, o1, o2, :O), tot+1:tot+n)
                        tot += n
                    elseif (o1 == :s && o2 == :p )
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n
                        data_info[(atX, o2, o1, :O)] = data_info[(atX, o1, o2, :O)]
                        tot += n
                    elseif (o1 == :p && o2 == :p )
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n*2
                        tot += n*2

                    elseif o1 == :s && o2 == :d
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n
                        data_info[(atX, o2, o1, :O)] = data_info[(atX, o1, o2, :O)]
                        tot += n

                    elseif o1 == :p && o2 == :d
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n
                        data_info[(atX, o2, o1, :O)] = data_info[(atX, o1, o2, :O)]
                        tot += n

                    elseif o1 == :d && o2 == :d
                        data_info[(atX, o1, o2, :O)] = tot+1:tot+n*2
                        data_info[(atX, o2, o1, :O)] = data_info[(atX, o1, o2, :O)]
                        tot += n*2

                    end
                end
            end
            return tot
        end

        if same_at #true onsite terms
           for o in orbs1
#                println("true onsite ", o)
                data_info[(at1, o, :A)] = [totH+1]
                totH += 1
            end
        end


        totHO = getonsite(at1, orbs1, totH, n_2body_onsite)

        if !(same_at) #need reverse if not same atom
            totHO = getonsite(at2, orbs2, totHO, n_2body_onsite)
        end

    elseif dim == 3 #3body
        
        totS = 0 #no 3body overlap terms

        at_list = [i for i in at_set]
        sort!(at_list)
        
        if length(at_list) == 1


            #permutations are trivial
            perm_ij = [[at_list[1], at_list[1], at_list[1]]]
            perm_on = [[at_list[1], at_list[1], at_list[1]]]
        elseif length(at_list) == 2


            #unique permutations
            perm_ij = [[at_list[1], at_list[1], at_list[2]] ,
                       [at_list[2], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[2]] ]

            perm_on = [[at_list[1], at_list[2], at_list[2]] ,
                        [at_list[1], at_list[1], at_list[2]] ,
                        [at_list[2], at_list[1], at_list[1]] ,
                        [at_list[2], at_list[1], at_list[2]] ]
            
        elseif length(at_list) == 3


            #all permutations exist hij
            perm_ij = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[1], at_list[3], at_list[2]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[2], at_list[3], at_list[1]] ,
                       [at_list[3], at_list[1], at_list[2]] ,
                       [at_list[3], at_list[2], at_list[1]] ]

            #onsite can flip last 2 atoms
            perm_on = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[3], at_list[1], at_list[2]]]

            
        else
            println("ERROR  get_data_info $dim $at_set $at_list")
        end
        



        function get3bdy(n, symb, start, at1, at2, at3)
            tot=start

            orbs1 = atoms[at1].orbitals
            orbs2 = atoms[at2].orbitals

            if at1 == at2
                same_at = true
            else
                same_at = false
            end

            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
                    
#                    push!(orbs, (o1, o2, symb))

                    if same_at
                        data_info[(at1, o1, at2, o2, at3,  symb)] = collect(tot+1:tot+n)
                        data_info[(at2, o2, at1, o1, at3,  symb)] = collect(tot+1:tot+n)
                    else                                                  #[1 2 3 4 5 6 7 8  9  10 11 12 13 14 15 16 17 18]
                        data_info[(at1, o1, at2, o2, at3,  symb)] = collect(tot+1:tot+n)
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  13 16 18 14 17 15]' #switch 2 4 and 3 6
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  ]' #switch 2 4 and 3 6
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 3 2 4 5 7 6 8 9]' #switch 2 4 and 3 6
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 3 2 4 5 7 6]' #switch 2 4 and 3 6
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 3 2  4 6 5]' #switch 2 4 and 3 6
#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1 3 2 4 6 5  7 9 8 ]' #switch 2 4 and 3 6

#                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1, 3, 2, 4, 5, 7, 6  ] #switch 2 4 and 3 6

                        data_info[(at2, o2, at1, o1, at3,  symb)] = tot .+ [1, 3, 2, 4, 5  ] #differs
                    end
                    
                    
                    tot += n

#                    if same_at
#                        data_info[(o2, o1, symb)] = data_info[(o1, o2, symb)]
#                    end
                        
                    
                end
            end
            return tot
        end

#        if at_list[2] == at_list[3]
#            same_at_on = true
#        else
#            same_at_on = false
 #       end
        

        function get3bdy_onsite(n, same_at,symb, start, at1, at2, at3)
#            if at2 == at3  #|| at1 == at2 || at1 == at3
#                same_at = true
#            else
#                same_at = false
#            end

            orbs1 = atoms[at1].orbitals

#            println("get3bdy_onsite $at1 $at2 $at3 $n")

            tot=start
            for o1 in orbs1
#                data_info[(at1, o1,at2, at3,  symb)] = collect(tot+1:tot+n)
#                data_info[(at1, o1,at3, at2,  symb)] = collect(tot+1:tot+n)

#                push!(orbs, (at1, o1,at2, at3,  symb))
#                push!(orbs, (at1, o1,at3, at2,  symb))

                if same_at
                    data_info[(at1, o1,at2, at3,  symb)] = collect(tot+1:tot+n)
                    data_info[(at1, o1,at3, at2,  symb)] = collect(tot+1:tot+n)
                else
                    data_info[(at1, o1,at2, at3,  symb)] = collect(tot+1:tot+n)
#                    data_info[(at1, o1,at3, at2,  symb)] = collect(tot+1:tot+n)
                    data_info[(at1, o1,at3, at2,  symb)] = tot .+ [1, 3, 2, 4]

#                    data_info[(at1, o1,at2, at3,  symb)] = collect(tot+1:tot+n)
#                    data_info[(at1, o1,at3, at2,  symb)] = tot .+ [1 3 2 4]'
                end
                tot += n                               #       1 2 3 4 5 6 7 8

            end
            return tot
        end
        
        tot_size = 0
        for p in perm_ij
            if p[1] == p[2]
                tot_size = get3bdy(n_3body_same, :H, tot_size, p[1], p[2], p[3])
            else
                tot_size = get3bdy(n_3body, :H, tot_size, p[1], p[2], p[3])
            end
        end

        for p in perm_on
            if  (p[1] == p[2] ||  p[2] == p[3] || p[1] == p[3])
                tot_size = get3bdy_onsite(n_3body_onsite_same,true, :O, tot_size, p[1], p[2], p[3]) #all diff
            else
                tot_size = get3bdy_onsite(n_3body_onsite,false, :O, tot_size, p[1], p[2], p[3]) #
            end
        end
    
        totHO = tot_size

    else
        println("error, only 2 or 3 body terms, you gave me : ", at_list)
    end
    return totHO ,totS, data_info, orbs

end            


function fix_format_change(datH, totHnew, dim, at_list, data_info)

    #this converts from old data_info to new data_info, setting the extra terms to zero.

    totH,totS, data_info_old, orbs = get_data_info_oldversion(at_list, dim)
    
    datHnew = zeros(totHnew)

    for k in keys(data_info)
        dnew = data_info[k]
        dold = data_info_old[k]

        n=length(dold)
        datHnew[dnew[1:n] ] = datH[dold]  #the extra terms are left as zero
    end
    
    return datHnew

end
