
""" Return the open reading frames c||responding to a given HXB2 index. """
function index2fram(i)
   frames = [] 
    if( ( 790<=i<=2292) || (5041<=i<=5619) || (8379<=i<=8469) || (8797<=i<=9417) )
        push!(frames, 1)
    end
    if( (5831<=i<=6045) || (6062<=i<=6310) || (8379<=i<=8653) )
        push!(frames, 2)
    end
    if( (2253<=i<=5096) || (5559<=i<=5850) || (5970<=i<=6045) || (6225<=i<=8795) )
        push!(frames, 3)
    end
    return frames
end;
""" Return the amino acid character c||responding to the input codon. """
function codon2aa(c; noq=false)
    # If all nucleotides are missing, return gap
    if(c[1]=="-" && c[2]=="-" && c[3]=="-")   return "-"
    # Else if some nucleotides are missing, return "?"
    elseif(c[1]=="-" || c[2]=="-" || c[3]=="-")
        if(noq) return "-" 
        else return "?" 
        end

    elseif( c[1] in ["W", "S", "M", "K", "R", "Y"] || c[2] in ["W", "S", "M", "K", "R", "Y"] ) return "X"

    elseif( c[1]=="T" )
        if( c[2]=="T" )
            if( c[3] in ["T", "C", "Y"] )       return "F"
            elseif( c[3] in ["A", "G", "R"] )   return "L"
            else                                return "X" 
            end

        elseif( c[2]=="C" )                     return "S"
        elseif( c[2]=="A" )
            if( c[3] in ["T", "C", "Y"] )       return "Y"
            elseif( c[3] in ["A", "G", "R"] )   return "*"
            else                                return "X" 
            end
            
        elseif( c[2]=="G" ) 
            if( c[3] in ["T", "C", "Y"] )       return "C"
            elseif  c[3]=="A"                   return "*"
            elseif  c[3]=="G"                   return "W"
            else                                return "X" 
            end
        else                                    return "X" 
        end 
            
    elseif( c[1]=="C" )
        if( c[2]=="T" )
            return "L"
        elseif( c[2]=="C" )                     return "P"
        elseif( c[2]=="A" ) 
            if( c[3] in ["T", "C", "Y"] )       return "H"
            elseif( c[3] in ["A", "G", "R"] )   return "Q"
            else                                return "X" 
            end
        elseif( c[2]=="G" )                     return "R"
        else                                    return "X" 
        end
        
    elseif(c[1]=="A")
        if( c[2]=="T" )
            if( c[3] in ["T", "C", "Y"] )       return "I"
            elseif( c[3] in ["A", "M", "W"] )   return "I"
            elseif( c[3]=="G" )                 return "M"
            else                                return "X" 
            end
        elseif( c[2]=="C" )                     return "T"
        elseif( c[2]=="A" )
            if( c[3] in ["T", "C", "Y"] )       return "N"
            elseif( c[3] in ["A", "G", "R"] )   return "K"
            else                                return "X" 
            end
        elseif( c[2]=="G" )
            if( c[3] in ["T", "C", "Y"] )       return "S"
            elseif( c[3] in ["A", "G", "R"] )   return "R"
            else                                return "X" 
            end
        else                                    return "X"
        end
            
    elseif( c[1]=="G" ) 
        if(c[2]=="T")                           return "V"
        elseif( c[2]=="C" )                     return "A"
        elseif( c[2]=="A" )
            if( c[3] in ["T", "C", "Y"] )       return "D"
            elseif( c[3] in ["A", "G", "R"] )   return "E"
            else                                return "X"
            end
        elseif(c[2]=="G")                       return "G"
        else                                    return "X"
        end
    else                                        return "X" 
    end
end;





""" get nonsynonymous mutation. This code requires only TF_seq.
"""
function get_nonsynonymous(i, mut, len_sites, site_DNA_norepeat, TF_seq)
    #len_sites = length(site_DNA_norepeat)
    i_HXB2 = site_DNA_norepeat[i]
    frames = index2fram(i_HXB2)
    count_nonsyno = 0

    for fr in frames
        pos = Int((i_HXB2-fr)%3)+1
        codon_indices = []
        if(i_HXB2 % 3 == 0)
            if(fr==3) 
             codon_indices = collect(i:(i+2))
            end
            if(fr==2)
             codon_indices = collect((i-1):(i+1))
            end
            if(fr==1)
             codon_indices = collect((i-2):i)
            end
        end
        
        if(i_HXB2 % 3 == 1)
            if(fr==3) 
             codon_indices = collect((i-1):(i+1))
            end
            if(fr==2)
             codon_indices = collect((i-2):i)
            end
            if(fr==1)
             codon_indices = collect(i:(i+2))
            end
        end
        
        if(i_HXB2 % 3 == 2)
            if(fr==3) 
             codon_indices = collect((i-2):i)
            end
            if(fr==2)
             codon_indices = collect(i:(i+2))
            end
            if(fr==1)
             codon_indices = collect((i-1):(i+1))
            end
        end        
        
        if( 0<codon_indices[1] && codon_indices[end]<=len_sites )
            TF_codon = TF_seq[codon_indices]
            mut_codon = copy(TF_codon)

            mut_codon[pos] = mut
            TF_AA = codon2aa(TF_codon)
            mut_AA = codon2aa(mut_codon)
            if(TF_AA != mut_AA)
                count_nonsyno += 1
            end
        else
            println("mutant at site $i_HXB2, mutation $mut in coon that does not terminate in alignment, assuming syn")
        end
    
    end
    return count_nonsyno
end;




""" Get N-linked glycosylation mutation 
"""
function get_glycosylation(i, mut, len_sites, site_DNA_norepeat, TF_seq)
    notX = ["P", "*", "?", "-"]
    SorT = ["S", "T"]
    codon_shifts = [-6, -3, 3, 6]
    shift = 2
    count_glycan = 0 # Positive => creation of Glycan, Negative => Deletion
    i_HXB2 = site_DNA_norepeat[i]
    if(6225<=i_HXB2 && i_HXB2 <= 8795 )
        #frames = index2fram(i_HXB2)
        #ONLY for Frames3!!!
	frames = [3]
	for fr in frames
            pos = Int((i_HXB2-fr)%3)+1
            codon_indices = []
            if(fr==1) # pos = 3
                codon_indices = (i-2+shift):(i+shift)
            end
            if(fr==2) # pos = 2
                codon_indices = (i-3+shift):(i-1+shift)
            end
            if(fr==3) # pos = 1
                codon_indices = (i-4+shift):(i-2+shift)
            end
            if( 0<codon_indices[1] && codon_indices[end]<=len_sites )
                TF_codon = TF_seq[codon_indices]
                mut_codon = copy(TF_codon)

                mut_codon[pos] = mut
                TF_AA = codon2aa(TF_codon)
                mut_AA = codon2aa(mut_codon)

                if(TF_AA != mut_AA)
                    #------- The case for INTRODUCING NEW Glycan mutation.  ------#
                    # The case for [, , !N, X, S/T] ->  [, , N, X, S/T]                
                    if((mut_AA == "N") && (TF_AA != "N") 
                            && codon_indices[end]+codon_shifts[4] <= len_sites )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[3]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[4]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA ∉ notX) && (TF_2_AA in SorT) )
                            count_glycan +=1
                        end
                    end
                    # The case for [, N, !X, S/T, ] ->  [, N, X, S/T, ]
                    if( (mut_AA ∉ notX) && (TF_AA ∈ notX)
                            && (codon_indices[1]  +codon_shifts[2] > 0)
                            && (codon_indices[end]+codon_shifts[3] <= len_sites) )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[2]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[3]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA == "N") && (TF_2_AA in SorT) )
                            count_glycan +=1
                        end
                    end
                    # The case for [N, X, !S/T, ,] ->  [N, !X, S/T, ,]
                    if( (mut_AA ∈ SorT) && (TF_AA ∉ SorT)
                            && (codon_indices[1]+codon_shifts[1] > 0) )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[1]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[2]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA == "N") && (TF_2_AA ∉ notX) )
                            count_glycan +=1
                        end
                    end

                    #------- The case for REMOVING EXISTED Glycan mutation.  ------#
                    # The case for [, , N, X, S/T] ->  [, , !N, X, S/T]                
                    if((TF_AA == "N") && (mut_AA != "N") 
                            && codon_indices[end]+codon_shifts[4] <= len_sites )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[3]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[4]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA ∉ notX) && (TF_2_AA in SorT) )
                            count_glycan -=1
                        end
                    end
                    # The case for [, N, X, S/T, ] -> [, N, !X, S/T, ]
                    if( (TF_AA ∉ notX) && (mut_AA ∈ notX)
                            && (codon_indices[1]  +codon_shifts[2] > 0)
                            && (codon_indices[end]+codon_shifts[3] <= len_sites) )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[2]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[3]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA == "N") && (TF_2_AA in SorT) )
                            count_glycan -=1
                        end
                    end
                    # The case for [N, !X, S/T, ,] -> [N, X, !S/T, ,]
                    if( (TF_AA ∈ SorT) && (mut_AA ∉ SorT)
                            && (codon_indices[1]+codon_shifts[1] > 0) )
                        TF_1_codon = TF_seq[collect(codon_indices) .+ codon_shifts[1]] # collect is need to add a codon_shift
                        TF_2_codon = TF_seq[collect(codon_indices) .+ codon_shifts[2]] # collect is need to add a codon_shift
                        TF_1_AA = codon2aa(TF_1_codon)
                        TF_2_AA = codon2aa(TF_2_codon)
                        if( (TF_1_AA == "N") && (TF_2_AA ∉ notX) )
                            count_glycan -=1
                        end
                    end
                end

            else
                println("mutant at site $i_HXB2, mutation $mut in coon that does not terminate in alignment, assuming syn")
            end
        end

    end
    return count_glycan
end;


""" Obtaining date in which the first mutaitons are observed. 
	0 => TF 
	>0 => first observed date of mutation
	-1 => neverobserved, i.e., monomorphic.
"""
function get_polymorphic(q, L, NUC, date_in, MSA, TF_seq_in)
    n_seq_sample = length(date_in)+1
    polymorphic = [-1 for i in 1:L for _ in 1:q ]; # this list should contain date of the first mutation is observed. 
    # get polymorphic: 0=>TF, >0=>mutation.
    n_count = 1
    for i in 1:L
        for a in NUC
            notyet_found_mutation = true
            n=2
            while(n<=n_seq_sample)
                if( notyet_found_mutation && (MSA[n][i]==a) )
                    if(TF_seq_in[i] != a)
                        polymorphic[n_count] = date_in[n-1]
                    end
                    if(TF_seq_in[i] == a)
                        polymorphic[n_count] = 0 
                    end
                    notyet_found_mutation = false
                    n = n_seq_sample # this condition is sufficient to break out the loop. 
                end
                n += 1
            end        
            n_count += 1
        end
    end
    return polymorphic
end;

""" Normalize selection coefficients such that S_i(mutant) <- S_i(mutant) - S_i(TF)
"""
function normalize_selection(fname_MPL, TF_seq_in)
    S_MPL = readdlm(fname_MPL)[:,1]
    S_MPL_TF = []
    for i in 1:L
        a = TF_seq_in[i]
        a_num = NUC2NUM(a)+1
        s_in = S_MPL[km(i, a_num, q)] 
        for _ in 1:q
            push!(S_MPL_TF, s_in)
        end
    end
    return S_MPL - S_MPL_TF
end;

