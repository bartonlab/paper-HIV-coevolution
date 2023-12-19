idx_HXB2_V1 = collect(6615:6693)
idx_HXB2_V2 = collect(6694:6761)
idx_HXB2_V3 = collect(7047:7124)
idx_HXB2_V4 = collect(7377:7479)
idx_HXB2_V5 = collect(7622:7708)
idx_HXB2_LD = collect(6901:7002)
idx_HXB2_MPER = collect(8560:8714)
idx_HXB2_CD4BS = collect(6783:6858)
# --------------------------------- CH848 Specific --------------------------------- #
idx_HXB2_half = collect(2000:9719)
idx_HXB2_ENV = collect(6225:8795)
idx_HXB2_POL = collect(2085:5096)
idx_HXB2_VIF = collect(5041:5619)
idx_HXB2_NEF = collect(8797:9417)
idx_HXB2_VPU = collect(6062:6310)
idx_HXB2_VPR = collect(5559:5850)
idx_HXB2_REV1 = collect(5970:6045)
idx_HXB2_REV2 = collect(8379:8653)
idx_HXB2_TAT1 = collect(5831:6045)
idx_HXB2_TAT2 = collect(8379:8469);
# ------------------------------- Glycan Related Analysis ------------------------------- #
notX = ["P", "*", "?", "-"]
SorT = ["S", "T"]




struct HXB2
    idx_nuc::Union{Number, Missing, Nothing}
    idx_AA::Union{String, Missing, Nothing}
    gene::Union{String, Missing, Nothing}
end

function get_num_nuc(x)
    x = string(x)
    num_nuc = -1
    if !has_letters(x)
        num_nuc = parse(Int, x)
    else
        this_nuc = match(r"\d+", x).match
        num_nuc = parse(Int, this_nuc)
    end
    return num_nuc
end;

idx_n2a(iL, iR) = Int(floor((iR+1-iL)/3)+1)
idx_n2a(iL, iR, fr) = Int(floor((iR+3-fr-iL)/3)+1)
function map_numNUC_to_numAA(x, fr_set)
    i_hxb2 = get_num_nuc(x) 
    gen_set, i_AA_set = ["" for _ in 1:3], ["" for _ in 1:3]
    for i_fr in 1:length(fr_set)
        fr = fr_set[i_fr]
        i_AA, gen = "", ""
        if(fr==1)
            if( 790 ≤ i_hxb2 ≤ 2292) i_AA = idx_n2a( 790, i_hxb2, fr); gen = "gag" end
            if(5041 ≤ i_hxb2 ≤ 5619) i_AA = idx_n2a(5041, i_hxb2, fr); gen = "vif" end
            if(8379 ≤ i_hxb2 ≤ 8424) i_AA = idx_n2a(8379, i_hxb2, fr); gen = "tat" end
            if(8797 ≤ i_hxb2 ≤ 9417) i_AA = idx_n2a(8797, i_hxb2, fr); gen = "nef" end      
        end
        if(fr==2)
            if(5831 ≤ i_hxb2 ≤ 6045) i_AA = idx_n2a(5831, i_hxb2, fr); gen = "tat" end
            if(6062 ≤ i_hxb2 ≤ 6310) i_AA = idx_n2a(6062, i_hxb2, fr); gen = "vpu" end
            if(8379 ≤ i_hxb2 ≤ 8653) i_AA = idx_n2a(8379, i_hxb2, fr); gen = "rev" end
        end
        if(fr==3)
            if(2085 ≤ i_hxb2 ≤ 5096) i_AA = idx_n2a(2085, i_hxb2, fr); gen = "pol" end
            if(5559 ≤ i_hxb2 ≤ 5850) i_AA = idx_n2a(5559, i_hxb2, fr); gen = "vpr" end
            if(5970 ≤ i_hxb2 ≤ 6045) i_AA = idx_n2a(5970, i_hxb2, fr); gen = "rev" end
            if(6225 ≤ i_hxb2 ≤ 8795) i_AA = idx_n2a(6225, i_hxb2, fr); gen = "env" end
        end
        gen_set[fr] = gen; i_AA_set[fr] = string(i_AA) 
    end
    return i_hxb2, i_AA_set, gen_set
end

function has_letters(str)
    # Check if the string contains any letters (A-Z or a-z)
    return occursin(r"[A-Za-z]", str)
end

function search_HXB2(hxb2_set, target_idx)
    for obj in hxb2_set
        if obj.idx_nuc == target_idx
            return obj.idx_AA, obj.gene
        end
    end
    return nothing, nothing  # Return nothing if target_a1 is not found
end

function get_total_index_variable_region(csv_raw_in)
    len_csv_raw_in = length(csv_raw_in.End_AA)
    csv_raw_in_index_set = []
    for n in 1:len_csv_raw_in
        for x in csv_raw_in.Start_AA[n]:csv_raw_in.End_AA[n]
            push!(csv_raw_in_index_set, x)
        end
    end
    csv_raw_in_index_set = copy(unique(csv_raw_in_index_set))
    return csv_raw_in_index_set
end

function print_matrix(matrix)
    rows, cols = size(matrix)
    for i in 1:rows
        for j in 1:cols
            @printf("%d ", matrix[i,j])
        end
        println()
    end
end

function get_when_mutation_occored(date_num, observed_nuc, a_WT)
    @assert length(date_num) == length(observed_nuc)
    n_entry = length(observed_nuc)
    a_MT_set = []; date_mut_set = []
    for n in 1:n_entry
        a_MT = observed_nuc[n]
        if(a_MT != a_WT && a_MT ∉ a_MT_set)
            push!(a_MT_set, a_MT)
            push!(date_mut_set, date_num[n])
        end
    end
    return a_MT_set, date_mut_set
end;


function replace_nothing_with_missing(arr)
    return [x === nothing ? missing : x for x in arr]
end

#@show mutant_types_set_AA;
""" Split a string argument arg that cosist in arg="string1Xstring2" to ["string1", "X", "string2"]
""" 
function split_mutant_string(s::String)
    m = match(r"(.*?)(\d+)(.*)", s)
    if m === nothing
        return []
    else
        return [m.captures...]
    end
end;

function my_replace_nothing(str_in)
    str_in_out = str_in
    if(occursin(r"nothing", str_in))
        str_in_out = replace(str_in, "nothing" => "0")
    end
    return str_in_out
end;

function check_mutant_is_in_reported(mut_in_origin, reported_mutant_data)
    flag_muched = false
    if(mut_in_origin != "" && mut_in_origin != "NA")
        mut_in_set = split(mut_in_origin, "|")
        for mut_in in mut_in_set
            s1,i,s2 = split_mutant_string(string(mut_in))
            s1_splited = split(s1, "/")
            s2_splited = split(s2, "/");
            for ref_in in reported_mutant_data
                ref_s1, ref_i, ref_s2 = split_mutant_string(ref_in)
                if(ref_s1!="X" && ref_s2!="X")
                    ref_eff = ref_s1*ref_i*ref_s2
                    for s1_eff in s1_splited
                        for s2_eff in s2_splited
                            s_eff = s1_eff * i * s2_eff
                            if(s_eff == ref_eff)
                                flag_muched = true
                            end
                        end
                    end
                end
                #     
                if(ref_s1=="X" && ref_s2!="X")
                    ref_eff = ref_i*ref_s2
                    for s2_eff in s2_splited
                        s_eff = i*s2_eff
                        if(s_eff == ref_eff)
                            flag_muched = true
                        end
                    end
                end
                #
                if(ref_s1!="X" && ref_s2=="X")
                    ref_eff = ref_s1*ref_i
                    for s1_eff in s1_splited
                        s_eff = s1_eff*i
                        if(s_eff == ref_eff)
                            flag_muched = true
                        end
                    end
                end
            end
        end
    end
    return flag_muched
end;

function get_true_false_variable_region(vec_AA_idx, csv_raw_in)
    vec_csv_raw_out = []
    for i in 1:length(vec_AA_idx)
        i_AA = vec_AA_idx[i]
        csv_raw_in_set = get_total_index_variable_region(csv_raw_in)
        if(string(i_AA) ∈ string.(csv_raw_in_set)) 
            push!(vec_csv_raw_out, true)
        else
            push!(vec_csv_raw_out, false)
        end
    end
    return vec_csv_raw_out
end;


function get_possible_glycan(i_mut, i_raw, date_mut_set, NUC, a_MT_set, n_poly_idx_max,  collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
#    idx_hxb2 = parse(Int, csv_index_and_TF.HXB2[i_raw])
    idx_hxb2 = parse(Int, match(r"(\d+)", csv_index_and_TF.HXB2[i_raw]).match)
    
    n_plus_glycan, n_minus_glycan, n_glycan_shift = zeros(Int,3), zeros(Int,3), zeros(Int,3)
    #if( 6225<= idx_hxb2 <= 8795)
        date_mut = date_mut_set[i_mut]
        data_after_mut = copy(data_num[collected_time .== date_mut, :])
        # Note, the ensemble should have only having mutation at idx_poly
        nuc_at_idx_poly = copy(data_after_mut[:, idx_poly])
        # this set shouldn't be empby because this site is polymorphic
        data_after_mut = copy(data_after_mut[nuc_at_idx_poly .== a_MT_set[i_mut] , :]) 
        data_before_mut = []
        if(date_mut>collected_time_unique[1])                
            n_time = collect(1:n_time_max)[collected_time_unique .== date_mut]
            date_earlier = collected_time_unique[n_time .- 1]
            data_before_mut = copy(data_num[collected_time .== date_earlier, :])
        else
            data_before_mut = copy(data_after_mut) # mutation took place at the earliest time
        end

        n_after = size(data_after_mut,1)
        n_before = size(data_before_mut,1)
        data_after_mut_extend = [x for i in 1:n_after , x in seq_TF]
        data_before_mut_extend = [x for i in 1:n_before , x in seq_TF]
        for i in 1:n_before 
            data_before_mut_extend[i, idx_only_poly] = [NUC[a+1] for a in data_before_mut[i, :]]
        end
        for i in 1:n_after 
            data_after_mut_extend[i, idx_only_poly] = [NUC[a+1] for a in data_after_mut[i, :]]
        end

        # Don't consider the edges
        if( 7<=i_raw && (i_raw + i_raw%3 + 2*3) <= n_poly_idx_max) 
        
            for i_fr in 1:length(this_frame_set)
                codon_location_set = [] # that should contains the 5 types of sites. 
                push!(codon_location_set, )
                frame_temp = this_frame_set[i_fr];
                aa_set_before = []
                aa_set_after = []
                codon_location = collect(1:3) # this index should extend 1 to to 2652
                if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                for k in -2:1:2 push!(codon_location_set, codon_location .+ 3*k) end     
            
                if(minimum(codon_location)>0 && maximum(codon_location)<=length(seq_TF))
                    for n in 1:n_before
                        aa_set_temp = []
                        for k in 1:5
                            if(minimum(codon_location_set[k])>0 && maximum((codon_location_set[k]))<=length(seq_TF))
                                x = join(data_before_mut_extend[n, codon_location_set[k]])
                                y = haskey(NUC2AA, x) ? NUC2AA[x] : "-"
                                push!(aa_set_temp, y)
                            end
                        end
                        push!(aa_set_before, join(aa_set_temp)) # Suppose we are only interested in the emergence of specific mutations, not statistical properties. 
                    end
                    for n in 1:n_after
                        aa_set_temp = []
                        for k in 1:5
                            if(minimum(codon_location_set[k])>0 && maximum((codon_location_set[k]))<=length(seq_TF))
                                x = join(data_after_mut_extend[n, codon_location_set[k]])
                                y = haskey(NUC2AA, x) ? NUC2AA[x] : "-"
                                push!(aa_set_temp, y)
                            end
                        end
                        push!(aa_set_after, join(aa_set_temp))
                    end
                    aa_set_before = unique(aa_set_before)
                    aa_set_after = unique(aa_set_after)

                    for x in aa_set_after
                        x_sp = split(x, "")
                        for y in aa_set_before
                            y_sp = split(y, "")
                            (n_glycan_shift_plus_temp, n_glycan_shift_minus_temp) = check_glycan_shield_hole_shift(x_sp, y_sp)
                            if(n_glycan_shift_plus_temp) n_plus_glycan[frame_temp] += 1 end
                            if(n_glycan_shift_minus_temp) n_minus_glycan[frame_temp] += 1 end
                            if(n_glycan_shift_minus_temp * n_glycan_shift_plus_temp) n_glycan_shift[frame_temp] += 1 end
                        end
                    end
                end
            end
        end # end for the if of i_raw range     
    #end
    return (n_plus_glycan, n_minus_glycan, n_glycan_shift)
end;

""" check_glycan_shield_hole_shift(x_sp, y_sp)
"""
function check_glycan_shield_hole_shift(x_sp, y_sp)
    idx_mut = 3 # The codon windows have a mutation at the middle.
    n_glycan_shift_plus_temp, n_glycan_shift_minus_temp = false, false
    #------------- Detect + N-linked Glycan ----------  
    #Here y and x are before and after mutation, respectively.
    # !N,!X,S/T -> N,!X,S/T 
    if(x_sp[idx_mut] == "N" && y_sp[idx_mut] != "N")
        if((y_sp[idx_mut+1] ∉ notX && y_sp[idx_mut+2] ∈ SorT) || 
           (x_sp[idx_mut+1] ∉ notX && x_sp[idx_mut+2] ∈ SorT)) 
            n_glycan_shift_plus_temp = true
        end
    end
    # N,X,S/T -> N,!X,S/T
    if(x_sp[idx_mut] ∉ notX && y_sp[idx_mut] ∈ notX)
        if( (y_sp[idx_mut-1] == "N" && y_sp[idx_mut+1] ∈ SorT) || 
            (x_sp[idx_mut-1] == "N" && x_sp[idx_mut+1] ∈ SorT) )
            n_glycan_shift_plus_temp = true
        end
    end

    # N,!X,!S/T -> N,!X,S/T 
    if(x_sp[idx_mut] ∈ SorT && y_sp[idx_mut] ∉ SorT )
        if( (y_sp[idx_mut-2] == "N" && y_sp[idx_mut-1] ∉ notX) ||
            (x_sp[idx_mut-2] == "N" && x_sp[idx_mut-1] ∉ notX))                            
            n_glycan_shift_plus_temp = true
        end
    end

    #------------- Detect - N-linked Glycan ----------  
    # !N,!X,S/T -> N,!X,S/T 
    if(y_sp[idx_mut] == "N" && x_sp[idx_mut] != "N")
        if((y_sp[idx_mut+1] ∉ notX && y_sp[idx_mut+2] ∈ SorT)||
           (x_sp[idx_mut+1] ∉ notX && x_sp[idx_mut+2] ∈ SorT))
            n_glycan_shift_minus_temp = true
        end
    end
    # N,X,S/T -> N,!X,S/T
    if(y_sp[idx_mut] ∉ notX && x_sp[idx_mut] ∈ notX)
        if( (y_sp[idx_mut-1] == "N" && y_sp[idx_mut+1] ∈ SorT)||
            (x_sp[idx_mut-1] == "N" && x_sp[idx_mut+1] ∈ SorT))
            n_glycan_shift_minus_temp = true
        end
    end
    # N,!X,!S/T -> N,!X,S/T 
    if(y_sp[idx_mut] ∈ SorT && x_sp[idx_mut] ∉ SorT )
        if( (y_sp[idx_mut-2] == "N" && y_sp[idx_mut-1] ∉ notX)||
            (x_sp[idx_mut-2] == "N" && x_sp[idx_mut-1] ∉ notX))
            n_glycan_shift_minus_temp = true
        end
    end
    return (n_glycan_shift_plus_temp, n_glycan_shift_minus_temp)
end;


"""
    Make_combination_of_mutations_with_genetic_background_w_glycan(csv_index_and_TF, seq_num_raw, hxb2csv)

TBW
"""
function Make_combination_of_mutations_with_genetic_background_w_glycan(csv_index_and_TF, seq_num_raw, hxb2csv)
    poly_idx = csv_index_and_TF.polymorphic
    seq_TF = string.(csv_index_and_TF.TF);
    L_TF = length(csv_index_and_TF.TF)# note this length can be different from the length of TF seq. in different csv files. 
    collected_time = seq_num_raw[:, 1]
    # Note mapping is 0=>'-', 1=>'A', .., 4=>'T', which are the same as my convention. 
    data_num = copy(seq_num_raw[:, 3:end]);
    collected_time_unique = sort(unique(collected_time), rev=false) # ensure that this is assending order. 
    n_time_max = length(collected_time_unique)
    idx_only_poly = csv_index_and_TF.polymorphic .!= "NA"
    mutant_types_set_nuc = [];
    
    # Look at only the site is !="NA"
    mutant_hxb2 = []
    mutant_nuc = []
    mutant_date_found = []
    
    mutant_gene = [[] for _ in 1:3]
    mutant_types_set_nuc = [[] for _ in 1:3]
    mutant_types_set_nuc_simple = [[] for _ in 1:3]
    mutant_types_set_AA = [[] for _ in 1:3]
    
    plus_glycan_set = [[] for _ in 1:3];
    minus_glycan_set = [[] for _ in 1:3];
    shifted_glycan_set = [[] for _ in 1:3]
    n_poly_idx_max = length(poly_idx) 
    for i_raw in 1:n_poly_idx_max
        if(poly_idx[i_raw] != "NA")
            # Find when first mutation occared 
            idx_hxb2 = csv_index_and_TF.HXB2[i_raw]
            a_WT = csv_index_and_TF.TF[i_raw]
            #@printf("\ni_raw=%d hxb2=%s a_WT=%s\n", i_raw, idx_hxb2, a_WT)
            idx_hxb2_num = parse(Int, match(r"\d+", idx_hxb2).match)
            this_frame_set, this_gene_set = index2frame(idx_hxb2_num)
            idx_poly = parse(Int, poly_idx[i_raw])+1         
            observed_nuc = copy(data_num[:, idx_poly])

            (a_MT_set, date_mut_set) = get_when_mutation_occored(collected_time, observed_nuc, a_WT) # Collect sequence that isn't same as a_WT at the subjected site
            perm = sortperm(a_MT_set) #ordering as -, A, C, G, T
            date_mut_set = copy(date_mut_set[perm])
            a_MT_set = copy(a_MT_set[perm])
            for i_mut in 1:length(date_mut_set) # this set shouldn't be empty set because we are looking at polymorpic site.
                a_mut = NUC[a_MT_set[i_mut]+1]
                (date_found, this_gene_set_out, mutation_tot_string_nuc_set_out, mutation_tot_string_nuc_simple_set_out, mutation_tot_string_AA_set_out) = get_possible_mutation_AA_and_NUC(
                    i_mut, i_raw, idx_hxb2, date_mut_set, NUC, a_MT_set, a_WT, n_poly_idx_max, 
                    collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, 
                    seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
                
                push!(mutant_hxb2, idx_hxb2)
                push!(mutant_nuc, a_mut)
                push!(mutant_date_found, date_found)
                for i_fr in 1:3 
                    #@printf("fr:%d gen:%s mut:%s mut_simple:%s mut_AA:%s\n", 
                    #i_fr, this_gene_set_out[i_fr], mutation_tot_string_nuc_set_out[i_fr], mutation_tot_string_nuc_simple_set_out[i_fr], mutation_tot_string_AA_set_out[i_fr]) 
                    push!(mutant_gene[i_fr], this_gene_set_out[i_fr], mutation_tot_string_nuc_simple_set_out[i_fr])
                    push!(mutant_types_set_nuc[i_fr], mutation_tot_string_nuc_set_out[i_fr]) 
                    push!(mutant_types_set_nuc_simple[i_fr], mutation_tot_string_nuc_simple_set_out[i_fr])
                    push!(mutant_types_set_AA[i_fr], mutation_tot_string_AA_set_out[i_fr])
                end

                (n_plus_glycan_set_out, n_minus_glycan_set_out, n_glycan_shift_set_out) = get_possible_glycan(i_mut, i_raw, date_mut_set, NUC, a_MT_set, n_poly_idx_max, collected_time, 
                    collected_time_unique, idx_poly, idx_only_poly, n_time_max, seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
                
                for i_fr in 1:3 
                    push!(plus_glycan_set[i_fr], n_plus_glycan_set_out[i_fr])
                    push!(minus_glycan_set[i_fr], n_minus_glycan_set_out[i_fr])
                    push!(shifted_glycan_set[i_fr], n_glycan_shift_set_out[i_fr])
                end
                
            end # eif(!(i_raw%3>0 && .... x)) # to avoid bounding error
        end # if(poly_idx[i_raw] != "NA")
    end #for of i_raw
    return (mutant_hxb2, mutant_nuc, mutant_date_found, mutant_gene, 
            mutant_types_set_AA, mutant_types_set_nuc, mutant_types_set_nuc_simple, plus_glycan_set, minus_glycan_set, shifted_glycan_set)
end;
function filtering_mutations_AA(mutant_types_set_AA)
    mutant_types_set_AA_filtered = []
    for x in mutant_types_set_AA
        if(x!="")
            letters_left_set = []; i_nuc_set = []; letters_right_set = []
            x_splited = split(x, "|")
            for y in x_splited
                match_obj = match(r"(.*?)(\d+)(.*)", y)
                if(match_obj != nothing) # if match_obj is only nothing, then return the empty lists which gives "" at the end.
                    letters_left = split(match_obj.captures[1], '/')
                    i_nuc_temp = match_obj.captures[2]
                    letters_right = split(match_obj.captures[3], '/')
                    [push!(letters_left_set, z) for z in letters_left if z != ""]
                    [push!(letters_right_set, z) for z in letters_right if z != ""]
                    push!(i_nuc_set, i_nuc_temp)
                end
            end
            out_str = join(sort(unique(letters_left_set)), "/") * join(unique(i_nuc_set)) * join(sort(unique(letters_right_set)), "/")
            push!(mutant_types_set_AA_filtered, out_str)
        else
            push!(mutant_types_set_AA_filtered, "")
        end
    end
    return mutant_types_set_AA_filtered
end;

""" This fucntion returns amino acid for each frame. Input files are only hxb2 index of mutaions, and the csv file including the TF's nuc and polymorphic indicies.  
"""
function get_TF_AA(csv_index_and_TF, mutant_hxb2)
    seq_TF = copy(csv_index_and_TF.TF);
    mutant_types_set_TF_AA = [[] for _ in 1:3]
    poly_idx = csv_index_and_TF.polymorphic
    n_poly_idx_max = length(poly_idx) 
    for i_raw in 1:n_poly_idx_max
        if(poly_idx[i_raw] != "NA")
            # Find when first mutation occared 
            idx_hxb2 = csv_index_and_TF.HXB2[i_raw]
            this_frame_set, this_gene_set = index2frame(idx_hxb2)
            a_WT = csv_index_and_TF.TF[i_raw]
            for i_fr in 1:3
                frame_temp = this_frame_set[i_fr];
                num_nuc, i_AA_set, gene_check = map_numNUC_to_numAA(idx_hxb2, frame_temp)
                i_AA = i_AA_set[frame_temp]
                #@assert(this_gene_set[i_fr] == gene_check[frame_temp]) ;# this should be same

                codon_location = collect(1:3) # this index should extend 1 to to 2652
                if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                if(minimum(codon_location)>0 && maximum(codon_location)<=length(seq_TF))
                    x_TF = join(seq_TF[codon_location])
                    AA_TF = haskey(NUC2AA, x_TF) ? NUC2AA[x_TF] : "-" 
                    push!(mutant_types_set_TF_AA[i_fr], AA_TF)
                else
                    push!(mutant_types_set_TF_AA[i_fr], "NA")
                end
            end
        else
            for i_fr in 1:3
                push!(mutant_types_set_TF_AA[i_fr], "NA")
            end
        end
    end
    
    mutant_types_set_TF_AA_out = [[] for _ in 1:3]        
    for i_hxb2 in mutant_hxb2
        idx = i_hxb2 .== csv_index_and_TF.HXB2
        if(count(idx)>0)
            for i_fr in 1:3
                push!(mutant_types_set_TF_AA_out[i_fr], mutant_types_set_TF_AA[i_fr][idx][1])
            end
        end
    end
    return mutant_types_set_TF_AA_out
end;

""" By giving any sequences, gives possible number of glycan holdes, shields and shifts
        get_glycan_plus_minus_shift_statistics(seq_TF)
"""
function get_glycan_plus_minus_shift_statistics(csv_index_and_TF)
    n_non_syn, n_N_add, n_N_rem, n_N_sht = 0, 0, 0, 0
    i_fr = 3
    seq_TF = copy(csv_index_and_TF.TF)
    for i_raw in 1:length(seq_TF)
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        #this_frame_set, this_gene_set = index2frame(idx_hxb2)
        #frame_temp = this_frame_set[end];
        frame_temp = 3
        codon_location_set = [] # that should contains the 5 types of sites. 
        push!(codon_location_set, )
        a_TF = seq_TF[i_raw]
        # Consider only non_gap sites
        if(a_TF != "-")
            codon_location = collect(1:3)
            if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
            if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
            if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end

            for k in -2:1:2 push!(codon_location_set, codon_location .+ 3*k) end
            aa_TF = []
            for k in 1:5
                if(minimum(codon_location_set[k])>0 && maximum((codon_location_set[k]))<=length(seq_TF))
                    x = join(seq_TF[codon_location_set[k]])
                    y = haskey(NUC2AA, x) ? NUC2AA[x] : "-"
                    push!(aa_TF, y)
                end
            end

            nuc_TF = seq_TF[i_raw]
            for nuc_MT in NUC 
                if(nuc_MT != nuc_TF)
                    seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT; aa_MT = []
                    for k in 1:5
                        if(minimum(codon_location_set[k])>0 && maximum((codon_location_set[k]))<=length(seq_MT))
                            x = join(seq_MT[codon_location_set[k]])
                            y = haskey(NUC2AA, x) ? NUC2AA[x] : "-"
                            push!(aa_MT, y)
                        end
                    end
                    #@printf("MT:%s, TF:%s\n", join(aa_MT), join(aa_TF))
                    if(join(aa_MT) != join(aa_TF))
                        n_non_syn += 1    
                        (bool_N_plus, bool_N_minus) = check_glycan_shield_hole_shift(aa_MT, aa_TF);
                        if(bool_N_plus) n_N_add += 1 end
                        if(bool_N_minus) n_N_rem += 1 end
                        if(bool_N_plus * bool_N_minus) n_N_sht += 1 end
                    end
                end 
            end
        end
    end
    return (n_non_syn, n_N_add, n_N_rem, n_N_sht)
end;

""" This function split the mutaion consiting XNY, where X and Y are amino acid letters and N is a digit to [X, N, Y]
"""
function extract_parts(s)
    results = []
    if(occursin(r"\d", s) && s!="")  # Check if the string contains a digit
        left_letters = match(r"^\D*", s).match  # Letters before the digit
        digit = match(r"\d+", s).match  # The digit
        right_letters = match(r"\D*$", s).match  # Letters after the digit

        push!(results, left_letters)
        push!(results, digit)
        push!(results, right_letters)
    else
        for _ in 1:3 push!(results, "NA" ) end  # No digit found, add false for all parts
        
    end
    return results
end;

""" This function replace the multiple occurance of letters as wildtype to the AA of the TF sequence.  
"""
# Example usage
function replacing_redundant_AA_by_TF(mutation_in, mutation_TF)
    strings = copy(mutation_in)
    to_be_replaced = copy(mutation_TF);
    output_letters = []
    for n in 1:length(strings)
        x = strings[n]

        if(x != "" && mutation_TF[n]!="NA")
            extracted_x = extract_parts(x)
            # check are there multiple letters at the left 
            splited_left_letters = split(extracted_x[1], "/")
            if(length(splited_left_letters)>1 || splited_left_letters[1]=="?")
                extracted_x[1] = to_be_replaced[n]
            end
            push!(output_letters, join(extracted_x) )
        else
            push!(output_letters, "" )
        end
        
    end
    return output_letters
end;


function get_possible_mutation_AA_and_NUC(i_mut, i_raw, idx_hxb2, date_mut_set, NUC, a_MT_set, a_WT, n_poly_idx_max, 
            collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
    date_found = date_mut_set[i_mut]
    a_mut = NUC[a_MT_set[i_mut]+1]
    this_gene_set_out = ["NA" for _ in 1:3]
    mutation_tot_string_nuc_set_out = ["NA" for _ in 1:3]
    mutation_tot_string_nuc_simple_set_out = ["NA" for _ in 1:3]
    mutation_tot_string_AA_set_out = ["NA" for _ in 1:3]

    date_mut = date_mut_set[i_mut]
    data_after_mut = copy(data_num[collected_time .== date_mut, :])
    # Note, the ensemble should have only having mutation at idx_poly
    nuc_at_idx_poly = copy(data_after_mut[:, idx_poly])
    # this set shouldn't be empby because this site is polymorphic
    data_after_mut = copy(data_after_mut[nuc_at_idx_poly .== a_MT_set[i_mut] , :]) 
    data_before_mut = []
    if(date_mut>collected_time_unique[1])                
        n_time = collect(1:n_time_max)[collected_time_unique .== date_mut]
        date_earlier = collected_time_unique[n_time .- 1]
        data_before_mut = copy(data_num[collected_time .== date_earlier, :])
    else
        data_before_mut = copy(data_after_mut) # mutation took place at the earliest time
    end
    n_after = size(data_after_mut,1)
    n_before = size(data_before_mut,1)
    data_after_mut_extend = [x for i in 1:n_after , x in seq_TF]
    data_before_mut_extend = [x for i in 1:n_before , x in seq_TF]
    for i in 1:n_before 
        data_before_mut_extend[i, idx_only_poly] = [NUC[a+1] for a in data_before_mut[i, :]]
    end
    for i in 1:n_after 
        data_after_mut_extend[i, idx_only_poly] = [NUC[a+1] for a in data_after_mut[i, :]]
    end
    # -- From here need to change the code and consider the difference of frames. --#

    this_gene = ""; mutation_tot_string_nuc = ""
    mutation_tot_string_nuc_simple = ""; mutation_tot_string_AA = ""        
    for i_fr in 1:length(this_frame_set)
        frame_temp = this_frame_set[i_fr];
        this_gene_set_out[frame_temp] = this_gene_set[i_fr]        
        num_nuc, i_AA_set, gene_check = map_numNUC_to_numAA(idx_hxb2, frame_temp)
        i_AA = i_AA_set[frame_temp]
        #@assert(this_gene_set[i_fr] == gene_check[frame_temp]) ;# this should be same
        
        codon_set_tot_after = []
        codon_set_tot_before = []
        codon_location = collect(1:3) # this index should extend 1 to to 2652
        if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
        if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
        if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
        
        if(minimum(codon_location)>0 && maximum(codon_location)<=length(seq_TF))
            [push!(codon_set_tot_before, x) for x in unique([join(data_before_mut_extend[n, codon_location]) for n in 1:n_before])]
            [push!(codon_set_tot_after, x) for x in unique([join(data_after_mut_extend[n, codon_location]) for n in 1:n_after])]
            codon_set_tot_before = copy(unique(codon_set_tot_before))
            codon_set_tot_after = copy(unique(codon_set_tot_after))
            """
            @printf("cdn_bfr:%s\n", join(codon_set_tot_before, ","))
            @printf("cdn_afr:%s\n", join(codon_set_tot_after, ","))
            """
            AA_set_after = [haskey(NUC2AA,x) ? NUC2AA[x] : "-" for x in codon_set_tot_after]
            AA_set_before = [haskey(NUC2AA,x) ? NUC2AA[x] : "-" for x in codon_set_tot_before]
            #@printf("fr:%d AA_before:%s AA_after:%s \n", i_fr, join(AA_set_before, "|"), join(AA_set_after, "|"))
            
            # --------------------- Making transition matrix that can change by a single substitution. ----------------------
            # ----> this process also consider the transition between mutations that observed at the same time. 
            len_codon_before = length(codon_set_tot_before)
            len_codon_after = length(codon_set_tot_after)
            paired_before = zeros(Int, len_codon_before)
            paired_after = zeros(Int, len_codon_after)
            table_paired = falses(len_codon_before, len_codon_after)
            for i_before in 1:len_codon_before
                x = codon_set_tot_before[i_before]
                x_split = split(x, "")
                for i_after in 1:len_codon_after
                    y = codon_set_tot_after[i_after]
                    y_split = split(y, "")
                    hamming_dist = sum( 1 .- kr.(x_split, y_split))
                    if(hamming_dist == 1)
                        # keep trucking paired/non-paired mutations
                        paired_before[i_before] += 1
                        paired_after[i_after] += 1
                        table_paired[i_before, i_after] = true
                    end
                end
            end
            # This table will be used if cannot make the pair of mutation in between two time steps. 
            table_paired_after = falses(len_codon_after, len_codon_after)
            for i_after in 1:len_codon_after
                x = codon_set_tot_after[i_after]
                x_split = split(x, "")
                for i_after_v2 in 1:len_codon_after
                    y = codon_set_tot_after[i_after_v2]
                    y_split = split(y, "")
                    if(i_after != i_after_v2)
                        hamming_dist = sum( 1 .- kr.(x_split, y_split))
                        if(hamming_dist == 1) table_paired_after[i_after, i_after_v2] = true end
                    end
                end
            end
            """
            @printf("Detected Codon: before>%s & after>%s\n", join(codon_set_tot_before, "/"), join(codon_set_tot_after, "/"))
            @printf("Detected AA: before>%s & after>%s\n", join(AA_set_before, "/"), join(AA_set_after, "/"))            
            @printf("Possible Transition (row: before, col: after): \n")
            print_matrix(table_paired)
            @printf("Possible Transition (After): \n")
            print_matrix(table_paired_after)
            """
            # ------------------------------- Adjoint paired/unpaired mutations ----------------------------------#
            mutation_tot_string_AA = ""
            mutation_tot_string_nuc = ""
            mutation_tot_string_nuc_simple = join(unique(data_before_mut_extend[:, i_raw]), "/") * idx_hxb2 * a_mut
            if(size(codon_set_tot_after, 1)>0 || size(codon_set_tot_before,1)>0)
                # --- Making the paires that has single nucleotide sustitution and regarding also the genetic background ----        
                # -------- Make the format to write out columns in the CSV file ----------#
                nuc_codon_idex = csv_index_and_TF.HXB2[codon_location]
                nuc_mut_string = join([nuc_codon_idex[1], nuc_codon_idex[3]], "-")
                # ------------------------ Formatting the mutations ----------------------#
                mutations_nuc_tot = []
                mutations_AA_tot = []
                for i_after in 1:len_codon_after
                    # Making a pair that the mutant (after) has pair(s) with a single substitution.
                    if(count(table_paired[:, i_after])>0)
                        after_str_nuc =  nuc_mut_string * codon_set_tot_after[i_after]
                        after_str_AA = string(i_AA) * AA_set_after[i_after]
                        before_str_nuc = join(codon_set_tot_before[table_paired[:, i_after]], "/")
                        before_str_AA = join(unique(AA_set_before[table_paired[:, i_after]]), "/")
                        push!(mutations_nuc_tot, before_str_nuc * after_str_nuc)
                        push!(mutations_AA_tot, before_str_AA * after_str_AA)
                    end
                end
                
                # ----- Treatments for mutations couldn't make a pair with a single substituions. ------#
                if(length(mutations_nuc_tot)==0)
                    for i_after in 1:len_codon_after
                        if(count(table_paired[:, i_after])==0) # if table_paired didn't much, use table_paired_after !
                            after_str_nuc =  nuc_mut_string * codon_set_tot_after[i_after]
                            after_str_AA = string(i_AA) * AA_set_after[i_after]                            
                            # --- Considering pairs that observed at the same times. --- #
                            after_str_nuc_bepaired = join(codon_set_tot_after[table_paired_after[i_after, :]], "/")
                            after_str_AA_bepaired = join(unique(AA_set_after[table_paired_after[i_after, :]]), "/")
                            
                            
                            AA_set_before_rem_gap = filter(x -> x != "-", AA_set_before)
                            codon_set_tot_before_rem_gap = filter(x -> x != "---", codon_set_tot_before)
                            tf_codon = join(seq_TF[codon_location])
                            if(after_str_nuc_bepaired != "")
                                push!(mutations_nuc_tot, after_str_nuc_bepaired * after_str_nuc)
                                push!(mutations_AA_tot, after_str_AA_bepaired * after_str_AA)
                            # --- If variant is just gap, then accept any mutations --- #
                            elseif(AA_set_after[i_after] == "-")
                                before_str_nuc = join(codon_set_tot_before_rem_gap, "/")
                                before_str_AA = join(unique(AA_set_before_rem_gap), "/")
                                push!(mutations_nuc_tot, before_str_nuc * after_str_nuc)
                                push!(mutations_AA_tot, before_str_AA * after_str_AA)
                            # --- If TF codons were still presence at the previous time, then assume mutation happend from T/F.  --- #                             
                            elseif( tf_codon ∈ codon_set_tot_before)
                                push!(mutations_nuc_tot, tf_codon * after_str_nuc)
                                a_TF = haskey(NUC2AA, tf_codon) ? NUC2AA[tf_codon] : "-"
                                push!(mutations_AA_tot, a_TF * after_str_AA)                                
                            # --- If there current is not gap but is gap at previous time, then assume gap  --- #
                            elseif(AA_set_after[i_after] != "-" && ("-" ∈ AA_set_before_rem_gap || [""]==AA_set_before))
                                push!(mutations_nuc_tot, "---" * after_str_nuc)
                                push!(mutations_AA_tot, "-" * after_str_AA)
                            else
                                push!(mutations_nuc_tot, "?" * after_str_nuc)
                                push!(mutations_AA_tot, "?" * after_str_AA)
                            end
                        end
                    end
                end
                
                mutation_tot_string_nuc = join(mutations_nuc_tot, "|")
                mutation_tot_string_AA = join(mutations_AA_tot, "|")
            end

            mutation_tot_string_nuc_set_out[frame_temp] = mutation_tot_string_nuc
            mutation_tot_string_nuc_simple_set_out[frame_temp] = mutation_tot_string_nuc_simple
            mutation_tot_string_AA_set_out[frame_temp] = mutation_tot_string_AA

            """
            @printf("fr:%d gene:%s poly_idx:%d HXB2_idx:%s nuc:%s detected:%dd\n", frame_temp, this_gene_set_out[frame_temp], i_raw, idx_hxb2, a_mut, date_found)
            @printf("GENE:%s AA:%s NUC:%s\n", mutation_tot_string_nuc_simple, mutation_tot_string_AA, mutation_tot_string_nuc)
            """
            end # end of frame.

        end
    
    return (date_found, this_gene_set_out, mutation_tot_string_nuc_set_out, 
        mutation_tot_string_nuc_simple_set_out, mutation_tot_string_AA_set_out) 
end;

function is_numeric_string(s)
    ismissing(s) && return false
    s == "NA" && return false
    try
        parse(Float64, s)
        return true
    catch e
        return false
    end
end;
function custom_mode(arr)
    # Filter out "NA"
    filtered_arr = filter(x -> x != "NA", arr)
    
    # If the filtered array is empty, return "NA"
    if isempty(filtered_arr)
        return "NA"
    end

    # Otherwise, return the mode
    return mode(filtered_arr)
end

function filter_empty_strings(arr)
    # If all elements are empty strings, return the array as-is
    if all(x -> x == "", arr)
        return arr
    end

    # Otherwise, filter out empty strings
    return filter(x -> x != "", arr)
end


function get_csv_RMs(f_dir, csv_index_and_TF, s_MPL_RMs, s_SL_RMs, dx_RMs, fname_key_RMs_and_Human)
    n_max_RMs = length(fname_key_RMs_and_Human);
    indices_marginalized = [i for i in 1:length(csv_index_and_TF.polymorphic_marginalized) if is_numeric_string(csv_index_and_TF.polymorphic_marginalized[i])]
    HXB2_detected = csv_index_and_TF.HXB2[indices_marginalized] # this can be used as an unique key. 

    # Temporaly store all elements to be easier the frequent access to find specific sites. 
    TF_set = []; HXB2_set = []; nucleotide_set = []
    consensus_set = []; edge_gap_set = []; exposed_set = []
    flanking_set = []; glycan_set = []; s_MPL_set = []
    s_SL_set = []; nonsynonymous_set = [] 
    for k in 1:length(fname_key_RMs_and_Human)
        f_key = fname_key_RMs_and_Human[k]
        fname_csv = f_dir *f_key* "-3-poly.csv"
        csv_temp = DataFrame(CSV.File(fname_csv));
        push!(HXB2_set, csv_temp.HXB2_index)
        push!(nucleotide_set, csv_temp.nucleotide)
        push!(TF_set, csv_temp.TF)
        push!(consensus_set, csv_temp.consensus)
        push!(edge_gap_set, csv_temp.edge_gap)
        push!(exposed_set, csv_temp.exposed)
        push!(flanking_set, csv_temp.flanking)
        push!(glycan_set, csv_temp.glycan) 
        push!(s_MPL_set, csv_temp.s_MPL)
        push!(s_SL_set, csv_temp.s_SL) 
        push!(nonsynonymous_set, csv_temp.nonsynonymous)
    end;

    # Identify what mutations were observed. 
    nuc_detected = []
    nuc_detected_RMs = []
    for x in HXB2_detected
        # ----- This is for Human -----#
        temp_vec = []
        if(x ∈ HXB2_set[1])
            [push!(temp_vec, y) for y in nucleotide_set[1][HXB2_set[1] .== x]]
        end
        if(x ∉ HXB2_set[1])
            push!(temp_vec, "")
        end
        push!(nuc_detected, join(unique(temp_vec), "/"))

        temp_vec = []
        for k in 2:n_max_RMs            
            if(x ∈ HXB2_set[k])
                [push!(temp_vec, y) for y in nucleotide_set[k][HXB2_set[k] .== x]]
            end
            if(x ∉ HXB2_set[k])
                push!(temp_vec, "")
            end        
        end
        push!(nuc_detected_RMs, join(unique(temp_vec), "/"))
    end

    polymorphic_index_global = []
    nucleotide_global = []
    HXB2_index_global = []
    nonsynonymous_global = []
    nonsynonymous_global_RMs = []

    TF_global = []; TF_global_RMs = []
    consensus_global = []; consensus_global_RMs = []
    exposed_global = []; exposed_global_RMs = []
    flanking_global = []; flanking_global_RMs = []
    glycan_global = []; glycan_global_RMs = []
    edge_gap_global = []; edge_gap_global_RMs = []
    s_MPL_global = [[] for k in 1:n_max_RMs]
    s_SL_global = [[] for k in 1:n_max_RMs]
    s_MPL_marginalized = []; s_SL_marginalized = []; dx_marginalized = []
    n_count = 0
    for n in 1:length(HXB2_detected)    
        nuc_detected_split = split(nuc_detected[n], "/")
        nuc_detected_RMs_split = split(nuc_detected_RMs[n], "/")
        for a in NUC
            n_count += 1
            if(a ∈ nuc_detected_split || a ∈ nuc_detected_RMs_split)
                push!(polymorphic_index_global, n-1) # this should be placed at a deeper level.
                push!(nucleotide_global, a)
                push!(HXB2_index_global, HXB2_detected[n])

                idx1 = HXB2_set[1] .== HXB2_detected[n]
                idx2 = nucleotide_set[1] .== a
                idx = idx1 .* idx2

                if(count(idx)>0)
                    [push!(TF_global, y) for y in TF_set[1][idx]]
                    [push!(nonsynonymous_global, y) for y in nonsynonymous_set[1][idx]]
                    [push!(consensus_global, y) for y in consensus_set[1][idx]]
                    [push!(edge_gap_global, y) for y in edge_gap_set[1][idx]]
                    [push!(exposed_global, y) for y in exposed_set[1][idx]]
                    [push!(flanking_global, y) for y in flanking_set[1][idx]]
                    [push!(glycan_global, y) for y in glycan_set[1][idx]]
                    [push!(s_MPL_global[1], y) for y in s_MPL_set[1][idx]]
                    [push!(s_SL_global[1], y) for y in s_SL_set[1][idx]]
                end
                if(count(idx)==0)
                    push!(TF_global, "NA")
                    push!(nonsynonymous_global, "NA")
                    push!(consensus_global, "NA")
                    push!(edge_gap_global, "NA")
                    push!(exposed_global, "NA")
                    push!(flanking_global, "NA")
                    push!(glycan_global, "NA")
                    push!(s_MPL_global[1], "NA")
                    push!(s_SL_global[1], "NA")
                end

                nonsynonymous_temp = []; consensus_temp = []
                edge_gap_temp = []; exposed_temp = []; TF_temp = []
                flanking_temp = []; glycan_temp = []
                for k in 2:n_max_RMs
                    idx1 = HXB2_set[k] .== HXB2_detected[n]
                    idx2 = nucleotide_set[k] .== a
                    idx = idx1 .* idx2
                    if(count(idx)>0)
                        [push!(TF_temp, y) for y in TF_set[k][idx]]
                        [push!(nonsynonymous_temp, y) for y in nonsynonymous_set[k][idx]]
                        [push!(consensus_temp, y) for y in consensus_set[k][idx]]
                        [push!(edge_gap_temp, y) for y in edge_gap_set[k][idx]]
                        [push!(exposed_temp, y) for y in exposed_set[k][idx]]
                        [push!(flanking_temp, y) for y in flanking_set[k][idx]]
                        [push!(glycan_temp, y) for y in glycan_set[k][idx]]
                        #
                        [push!(s_MPL_global[k], y) for y in s_MPL_set[k][idx]]
                        [push!(s_SL_global[k], y) for y in s_SL_set[k][idx]]
                    end
                    if(count(idx)==0)
                        push!(TF_temp, "NA")
                        push!(nonsynonymous_temp, "NA")
                        push!(consensus_temp, "NA")
                        push!(edge_gap_temp, "NA")
                        push!(exposed_temp, "NA")
                        push!(flanking_temp, "NA")
                        push!(glycan_temp, "NA")
                        #
                        push!(s_MPL_global[k], "NA")
                        push!(s_SL_global[k], "NA")                    
                    end
                end
                push!(TF_global_RMs, custom_mode(TF_temp))
                push!(nonsynonymous_global_RMs, custom_mode(nonsynonymous_temp))
                push!(consensus_global_RMs, custom_mode(consensus_temp))
                push!(edge_gap_global_RMs, custom_mode(edge_gap_temp))
                push!(exposed_global_RMs, custom_mode(exposed_temp))
                push!(flanking_global_RMs, custom_mode(flanking_temp))
                push!(glycan_global_RMs, custom_mode(glycan_temp))

                push!(s_MPL_marginalized, s_MPL_RMs[n_count])
                push!(s_SL_marginalized, s_SL_RMs[n_count])
                push!(dx_marginalized, dx_RMs[n_count])            

             end # (a ∈ nuc_detected_split || a ∈ nuc_detected_RMs_split)
        end # end NUC
    end
    return (polymorphic_index_global, HXB2_index_global, nucleotide_global, 
    TF_global, TF_global_RMs, nonsynonymous_global, nonsynonymous_global_RMs, 
    edge_gap_global, edge_gap_global_RMs, exposed_global, exposed_global_RMs, 
    flanking_global, flanking_global_RMs, glycan_global, glycan_global_RMs, consensus_global, consensus_global_RMs,   
    dx_marginalized, s_MPL_marginalized, s_SL_marginalized, s_MPL_global, s_SL_global)
end;


function get_jointed_RMs_for_CSV(csv_selection, csv_index_and_TF, fname_key_RMs_and_Human, fname_key_RMs)
    fr = 3 
    vec_gene_idx = []
    vec_nuc_idx = []; vec_AA_idx = []; vec_TF = []
    for x in csv_selection.HXB2
        num_nuc = get_num_nuc(x)
        frame, gene = index2frame(num_nuc)
        num_nuc, num_AA, this_gene = map_numNUC_to_numAA(x, frame)
        push!(vec_nuc_idx, num_nuc)
        push!(vec_AA_idx, num_AA[fr])
        push!(vec_gene_idx, this_gene[fr])
    end;

    @printf("gene: %s\n", join(unique(vec_gene_idx), ", "));
    #@printf("size: %d %d %d\n", size(HXB2_set,1), size(vec_nuc_idx,1), size(vec_AA_idx,1)) 
    @printf("max: %d, %d\n", maximum(skipmissing(vec_nuc_idx)), maximum(skipmissing(vec_nuc_idx))) 
    @printf("min: %d, %d\n", minimum(skipmissing(vec_nuc_idx)), minimum(skipmissing(vec_nuc_idx)));

    vec_CH103_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_CH103binding)
    vec_CH235_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_CH235binding)
    vec_CD4contact_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_CD4contact)
    vec_V1 = get_true_false_variable_region(vec_AA_idx, csv_raw_V1_HV)
    vec_V2 = get_true_false_variable_region(vec_AA_idx, csv_raw_V2_HV)
    vec_V3 = get_true_false_variable_region(vec_AA_idx, csv_raw_V3_HV)
    vec_V4 = get_true_false_variable_region(vec_AA_idx, csv_raw_V4_HV)
    vec_V5 = get_true_false_variable_region(vec_AA_idx, csv_raw_V5_HV)
    vec_LoopD = get_true_false_variable_region(vec_AA_idx, csv_raw_LoopD);

    #csv_temp.HXB2_index <=> csv_selection.HXB2
    mutant_types_set_nuc_simple = [[] for _ in 1:length(csv_selection.HXB2)]
    mutant_types_set_AA_filtered = [[] for _ in 1:length(csv_selection.HXB2)]
    Nlinked_plus_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    Nlinked_minus_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    Nlinked_shift_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    for k in 1:length(fname_key_RMs)
        key_RM = fname_key_RMs[k]
        csv_temp = DataFrame(CSV.File("../out/" *key_RM* "-3-poly_sorted_with_mutant.csv"));
        for i in 1:length(csv_selection.HXB2)
            site_hxb2 = csv_selection.HXB2[i]
            nuc_this = csv_selection.nuc[i]
            idx1 = csv_temp.HXB2_index .== site_hxb2
            idx2 = csv_temp.nucleotide .== nuc_this
            idx = idx1 .* idx2 
            if(count(idx)>0)
                [push!(mutant_types_set_AA_filtered[i], convert(String, x)) for x in csv_temp.mutants_AA_fr3[idx]]
                [push!(mutant_types_set_nuc_simple[i], convert(String, x)) for x in csv_temp.mutants_nuc_filtered[idx]]
            else
                push!(mutant_types_set_AA_filtered[i], "")
                push!(mutant_types_set_nuc_simple[i], "")
            end
        end
    end
    for k in 1:length(fname_key_RMs_and_Human)
        key_RM = fname_key_RMs_and_Human[k]
        csv_temp = DataFrame(CSV.File("../out/" *key_RM* "-3-poly_sorted_with_mutant.csv"));
        @printf("Processing %s, done.\n", key_RM) 
        for i in 1:length(csv_selection.HXB2)
            site_hxb2 = csv_selection.HXB2[i]
            nuc_this = csv_selection.nuc[i]
            idx1 = csv_temp.HXB2_index .== site_hxb2
            idx2 = csv_temp.nucleotide .== nuc_this
            idx = idx1 .* idx2 
            if(count(idx)>0)
                push!(Nlinked_plus_set[k],  csv_temp.N_linked_glycan_plus_fr3[idx][1])
                push!(Nlinked_minus_set[k], csv_temp.N_linked_glycan_minus_fr3[idx][1])
                push!(Nlinked_shift_set[k], csv_temp.N_linked_glycan_shift_fr3[idx][1])
            else
                push!(Nlinked_plus_set[k], "")
                push!(Nlinked_minus_set[k], "")
                push!(Nlinked_shift_set[k], "")
            end
        end
    end

    mutant_types_set_AA_filtered = copy([unique(x) for x in mutant_types_set_AA_filtered])
    mutant_types_set_nuc_simple = copy([unique(x) for x in mutant_types_set_nuc_simple]);

    mutant_types_set_AA_filtered_temp = []
    for x in mutant_types_set_AA_filtered
        letters_left_set = []
        letters_right_set = []
        i_nuc_set = []
        w = ""
        count_flag = false
        for y in x
            if(y != "")
                count_flag = true
                match_obj = match(r"(.*?)(\d+)(.*)", y)
                letters_left = split(match_obj.captures[1], '/')
                i_nuc_temp = parse(Int, match_obj.captures[2])
                letters_right = split(match_obj.captures[3], '/')
                [push!(letters_left_set, z) for z in letters_left]
                [push!(letters_right_set, z) for z in letters_right]
                push!(i_nuc_set, i_nuc_temp)
            end
        end
        if(count_flag)
            w = join(sort(unique(letters_left_set)), "/") * join(unique(i_nuc_set)) * join(sort(unique(letters_right_set)), "/")
        else
           w = "" 
        end
        push!(mutant_types_set_AA_filtered_temp, w) 
    end

    mutant_types_set_nuc_simple_temp = []
    for x in mutant_types_set_nuc_simple
        letters_left_set = []
        letters_right_set = []
        i_nuc_set = []
        w = ""
        count_flag = false
        for y in x        
            if(y != "")
                count_flag = true
                match_obj = match(r"(.*?)(\d+)(.*)", y)
                letters_left = split(match_obj.captures[1], '/')
                i_nuc_temp = parse(Int, match_obj.captures[2])
                letter_right = match_obj.captures[3]
                [push!(letters_left_set, z) for z in letters_left]
                push!(letters_right_set, letter_right)
                push!(i_nuc_set, i_nuc_temp)
            end
        end
        if(count_flag)
            w = join(sort(unique(letters_left_set)), "/") * join(unique(i_nuc_set)) * join(sort(unique(letters_right_set)))
        else
           w = "" 
        end
        push!(mutant_types_set_nuc_simple_temp, w) 
    end
    mutant_types_set_nuc_simple = copy(mutant_types_set_nuc_simple_temp);
    mutant_types_set_AA_filtered = copy(mutant_types_set_AA_filtered_temp);

    # --- Check whether mutations have been reported or not. ---# 
    reported_mutant_data = copy(String.(csv_raw_resist_mut_CH103[:, 1]))
    bool_resist_mut_CH103 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_resist_mut_CH235[:, 1]))
    bool_resist_mut_CH235 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_resist_strain_specific_Abs_CH505[:, 1]))
    bool_resist_strain_specific_Abs_CH505 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_common_mut_SHIV_CH505[:, 1]))
    bool_common_mut_SHIV_CH505 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    reversion_true_false = csv_selection.consensus_RMs .== csv_selection.nuc .!= csv_selection.TF_RMs;
    
    return (vec_gene_idx, vec_nuc_idx, vec_AA_idx, vec_TF, 
    vec_CH103_binding, vec_CH235_binding, vec_CD4contact_binding, 
    vec_V1, vec_V2, vec_V3, vec_V4, vec_V5, vec_LoopD, 
    mutant_types_set_nuc_simple, mutant_types_set_AA_filtered, 
    Nlinked_plus_set, Nlinked_minus_set, Nlinked_shift_set,
    bool_resist_mut_CH103, bool_resist_mut_CH235, bool_resist_strain_specific_Abs_CH505, 
    bool_common_mut_SHIV_CH505, reversion_true_false) 
end;


function get_jointed_RMs_for_CSV_SHIV848(csv_selection, csv_index_and_TF, fname_key_RMs_and_Human, fname_key_RMs)
    fr = 3 
    vec_gene_idx = []
    vec_nuc_idx = []; vec_AA_idx = []; vec_TF = []
    for x in csv_selection.HXB2
        num_nuc = get_num_nuc(x)
        frame, gene = index2frame(num_nuc)
        num_nuc, num_AA, this_gene = map_numNUC_to_numAA(x, frame)
        push!(vec_nuc_idx, num_nuc)
        push!(vec_AA_idx, num_AA[fr])
        push!(vec_gene_idx, this_gene[fr])
    end;

    @printf("gene: %s\n", join(unique(vec_gene_idx), ", "));
    #@printf("size: %d %d %d\n", size(HXB2_set,1), size(vec_nuc_idx,1), size(vec_AA_idx,1)) 
    @printf("max: %d, %d\n", maximum(skipmissing(vec_nuc_idx)), maximum(skipmissing(vec_nuc_idx))) 
    @printf("min: %d, %d\n", minimum(skipmissing(vec_nuc_idx)), minimum(skipmissing(vec_nuc_idx)));

    vec_DH270_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_DH270_binding)
    vec_DH272_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_DH272_binding)
    vec_DH475_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_DH475_binding)
    vec_CD4contact_binding = get_true_false_variable_region(vec_AA_idx, csv_raw_CD4contact)
    vec_V1 = get_true_false_variable_region(vec_AA_idx, csv_raw_V1_HV)
    vec_V2 = get_true_false_variable_region(vec_AA_idx, csv_raw_V2_HV)
    vec_V3 = get_true_false_variable_region(vec_AA_idx, csv_raw_V3_HV)
    vec_V4 = get_true_false_variable_region(vec_AA_idx, csv_raw_V4_HV)
    vec_V5 = get_true_false_variable_region(vec_AA_idx, csv_raw_V5_HV)
    vec_LoopD = get_true_false_variable_region(vec_AA_idx, csv_raw_LoopD);

    #csv_temp.HXB2_index <=> csv_selection.HXB2
    mutant_types_set_nuc_simple = [[] for _ in 1:length(csv_selection.HXB2)]
    mutant_types_set_AA_filtered = [[] for _ in 1:length(csv_selection.HXB2)]
    Nlinked_plus_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    Nlinked_minus_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    Nlinked_shift_set = [[] for _ in 1:length(fname_key_RMs_and_Human)]
    for k in 1:length(fname_key_RMs)
        key_RM = fname_key_RMs[k]
        csv_temp = DataFrame(CSV.File("../out/" *key_RM* "-3-poly_sorted_with_mutant.csv"));
        for i in 1:length(csv_selection.HXB2)
            site_hxb2 = csv_selection.HXB2[i]
            nuc_this = csv_selection.nuc[i]
            idx1 = csv_temp.HXB2_index .== site_hxb2
            idx2 = csv_temp.nucleotide .== nuc_this
            idx = idx1 .* idx2 
            if(count(idx)>0)
                [push!(mutant_types_set_AA_filtered[i], convert(String, x)) for x in csv_temp.mutants_AA_fr3[idx]]
                [push!(mutant_types_set_nuc_simple[i], convert(String, x)) for x in csv_temp.mutants_nuc_filtered[idx]]
            else
                push!(mutant_types_set_AA_filtered[i], "")
                push!(mutant_types_set_nuc_simple[i], "")
            end
        end
    end
    for k in 1:length(fname_key_RMs_and_Human)
        key_RM = fname_key_RMs_and_Human[k]
        csv_temp = DataFrame(CSV.File("../out/" *key_RM* "-3-poly_sorted_with_mutant.csv"));
        @printf("Processing %s, done.\n", key_RM) 
        for i in 1:length(csv_selection.HXB2)
            site_hxb2 = csv_selection.HXB2[i]
            nuc_this = csv_selection.nuc[i]
            idx1 = csv_temp.HXB2_index .== site_hxb2
            idx2 = csv_temp.nucleotide .== nuc_this
            idx = idx1 .* idx2 
            if(count(idx)>0)
                push!(Nlinked_plus_set[k],  csv_temp.N_linked_glycan_plus_fr3[idx][1])
                push!(Nlinked_minus_set[k], csv_temp.N_linked_glycan_minus_fr3[idx][1])
                push!(Nlinked_shift_set[k], csv_temp.N_linked_glycan_shift_fr3[idx][1])
            else
                push!(Nlinked_plus_set[k], "")
                push!(Nlinked_minus_set[k], "")
                push!(Nlinked_shift_set[k], "")
            end
        end
    end

    mutant_types_set_AA_filtered = copy([unique(x) for x in mutant_types_set_AA_filtered])
    mutant_types_set_nuc_simple = copy([unique(x) for x in mutant_types_set_nuc_simple]);

    mutant_types_set_AA_filtered_temp = []
    for x in mutant_types_set_AA_filtered
        letters_left_set = []
        letters_right_set = []
        i_nuc_set = []
        w = ""
        count_flag = false
        for y in x
            if(y != "")
                count_flag = true
                match_obj = match(r"(.*?)(\d+)(.*)", y)
                letters_left = split(match_obj.captures[1], '/')
                i_nuc_temp = parse(Int, match_obj.captures[2])
                letters_right = split(match_obj.captures[3], '/')
                [push!(letters_left_set, z) for z in letters_left]
                [push!(letters_right_set, z) for z in letters_right]
                push!(i_nuc_set, i_nuc_temp)
            end
        end
        if(count_flag)
            w = join(sort(unique(letters_left_set)), "/") * join(unique(i_nuc_set)) * join(sort(unique(letters_right_set)), "/")
        else
           w = "" 
        end
        push!(mutant_types_set_AA_filtered_temp, w) 
    end

    mutant_types_set_nuc_simple_temp = []
    for x in mutant_types_set_nuc_simple
        letters_left_set = []
        letters_right_set = []
        i_nuc_set = []
        w = ""
        count_flag = false
        for y in x        
            if(y != "")
                count_flag = true
                match_obj = match(r"(.*?)(\d+)(.*)", y)
                letters_left = split(match_obj.captures[1], '/')
                i_nuc_temp = parse(Int, match_obj.captures[2])
                letter_right = match_obj.captures[3]
                [push!(letters_left_set, z) for z in letters_left]
                push!(letters_right_set, letter_right)
                push!(i_nuc_set, i_nuc_temp)
            end
        end
        if(count_flag)
            w = join(sort(unique(letters_left_set)), "/") * join(unique(i_nuc_set)) * join(sort(unique(letters_right_set)))
        else
           w = "" 
        end
        push!(mutant_types_set_nuc_simple_temp, w) 
    end
    mutant_types_set_nuc_simple = copy(mutant_types_set_nuc_simple_temp);
    mutant_types_set_AA_filtered = copy(mutant_types_set_AA_filtered_temp);

    # --- Check whether mutations have been reported or not. ---# 
    reported_mutant_data = copy(String.(csv_raw_resist_mut_DH270[:, 1]))
    bool_resist_mut_DH270 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_resist_mut_DH272[:, 1]))
    bool_resist_mut_DH272 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_resist_mut_DH475[:, 1]))
    bool_resist_mut_DH475 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_resist_strain_specific_Abs_CH848[:, 1]))
    bool_resist_strain_specific_Abs_CH848 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    #
    reported_mutant_data = copy(String.(csv_raw_common_mut_SHIV_CH848[:, 1]))
    bool_common_mut_SHIV_CH848 = [check_mutant_is_in_reported(my_replace_nothing(string(x)), reported_mutant_data) for x in mutant_types_set_AA_filtered]
    reversion_true_false = csv_selection.consensus_RMs .== csv_selection.nuc .!= csv_selection.TF_RMs;
    
    return (vec_gene_idx, vec_nuc_idx, vec_AA_idx, vec_TF, 
    vec_DH270_binding, vec_DH272_binding, vec_DH475_binding, vec_CD4contact_binding, 
    vec_V1, vec_V2, vec_V3, vec_V4, vec_V5, vec_LoopD, 
    mutant_types_set_nuc_simple, mutant_types_set_AA_filtered, 
    Nlinked_plus_set, Nlinked_minus_set, Nlinked_shift_set,
    bool_resist_mut_DH270, bool_resist_mut_DH272, bool_resist_mut_DH475, bool_resist_strain_specific_Abs_CH848, 
    bool_common_mut_SHIV_CH848, reversion_true_false) 
end;

function check_N_linked_glycan(aa_set_temp)
    notX = ["P", "*", "?", "-"]
    SorT = ["S", "T"]
    flag_N_linked_glycan = false
    if(aa_set_temp[1] == "N")
        if(aa_set_temp[3] ∈ SorT)
            if(aa_set_temp[2] ∉ notX)
                flag_N_linked_glycan = true
            end
        end
    end
    return flag_N_linked_glycan
end;       

""" get_x_fold(vec_in, idx_sel, csv_index_and_TF, idx_type; reversion=false)
"""
function get_x_fold(vec_in, idx_sel, csv_index_and_TF, idx_type; n_null=0, N_null=0)
    # Franction for null or hypothetical mutations
    if(n_null==0)
        n_null = count( [x ∈ idx_type for x in extract_integer.(csv_index_and_TF.HXB2[csv_index_and_TF.TF .!= "-", 1]) ] );
    end
    if(N_null==0)
        N_null = length(csv_index_and_TF[csv_index_and_TF.TF .!= "-", 1])
    end
    α_null = n_null / N_null
    
    N_sel = count(vec_in)
    n_sel = count(vec_in[idx_sel])
    α_sel, x_fold = 0, 0
    if(N_sel > 0) α_sel = n_sel / N_sel end 
    if(α_null > 0) x_fold = α_sel / α_null end
    
    #@printf("N_sel:%d\tn_sel:%d\tN_null:%d\tn_null:%d", N_sel, n_sel, N_null, n_null)
    return (n_sel, n_null, N_sel, x_fold)
end;


""" get_x_fold(vec_in, idx_sel, csv_index_and_TF, idx_type; reversion=false)
"""
function get_x_fold(n_sel, N_sel, n_null, N_null)
    # Franction for null or hypothetical mutations
    α_null = n_null / N_null
    α_sel, x_fold = 0, 0
    if(N_sel > 0) α_sel = n_sel / N_sel end 
    if(α_null > 0) x_fold = α_sel / α_null end
    #@printf("N_sel:%d\tn_sel:%d\tN_null:%d\tn_null:%d", N_sel, n_sel, N_null, n_null)
    return (n_sel, n_null, N_sel, x_fold)
end;


""" Count the number of nonsynonymous mutations and nonsynonymous mutations restricted to a specific region or types.
    get_num_of_nonsyn(csv_index_and_TF, idx_type)
"""
function get_num_of_nonsyn(csv_index_and_TF, idx_type)
    n_nsyn, n_nsyn_restricted = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    for i_raw in 1:length(seq_TF)
        nuc_TF = seq_TF[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            for nuc_MT in NUC
                if(nuc_MT != nuc_TF)
                    seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT        
                    flag_nsyn, flag_nsyn_restricted = false, false
                    for i_fr in 1:3
                        codon_location = collect(1:3)
                        if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                        if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                        if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                        codon_TF = join(seq_TF[codon_location])
                        aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                        codon_MT = join(seq_MT[codon_location])
                        aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                        if(aa_MT != aa_TF) 
                            flag_nsyn = true
                            if(idx_hxb2 ∈ idx_type)
                                flag_nsyn_restricted = true
                            end
                        end
                    end
                    if(flag_nsyn) n_nsyn += 1 end
                    if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
                end
            end
        end
    end
    return (n_nsyn, n_nsyn_restricted)
end;

# idx_type is not necessary for this function. --> need to fix this!
function get_num_of_nonsyn_reversion(csv_index_and_TF)
    n_nsyn, n_nsyn_restricted = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    seq_consensus = copy(csv_index_and_TF.consensus)
    for i_raw in 1:length(seq_TF)
        nuc_TF = seq_TF[i_raw]
        nuc_consensus = seq_consensus[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            for nuc_MT in NUC
                if(nuc_MT != nuc_TF) 
                    seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT
                    flag_nsyn, flag_nsyn_restricted = false, false
                    for i_fr in 1:3
                        codon_location = collect(1:3)
                        if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                        if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                        if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                        #codon_location
                        codon_TF = join(seq_TF[codon_location])
                        aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                        codon_MT = join(seq_MT[codon_location])
                        aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                        if(aa_MT != aa_TF) 
                            flag_nsyn = true
                            if(nuc_consensus == nuc_MT)
                                flag_nsyn_restricted = true
                            end
                        end
                    end
                    if(flag_nsyn) n_nsyn += 1 end
                    if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
                end
            end
        end
    end 
    return (n_nsyn, n_nsyn_restricted)
end;



function get_n_sel_and_N_sel(idx_type_in, csv_raw_in, csv_index_and_TF, idx_type)
    # The following computation is similar to the n_null and N_null computation
    n_nsyn, n_nsyn_restricted = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    L_TF = length(csv_index_and_TF.HXB2)
    i_eff_restricted_max = count(idx_significant)
    for i_eff in 1:count(idx_type_in)
        i_raw = collect(1:L_TF)[csv_index_and_TF.HXB2 .== csv_raw_in.HXB2_index[idx_type_in][i_eff]][1] # This line is clitical to get the HXB2 index in TF seq.
        nuc_MT = csv_raw_in.nucleotide[idx_type_in][i_eff]
        nuc_TF = seq_TF[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            if(nuc_MT != nuc_TF)
                seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT        
                flag_nsyn, flag_nsyn_restricted = false, false
                for i_fr in 1:3
                    codon_location = collect(1:3)
                    if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                    if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                    if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                    codon_TF = join(seq_TF[codon_location])
                    aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                    codon_MT = join(seq_MT[codon_location])
                    aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                    if(aa_MT != aa_TF) 
                        flag_nsyn = true
                        if(i_eff <= i_eff_restricted_max)
                            flag_nsyn_restricted = true
                        end
                    end
                end
                if(flag_nsyn) n_nsyn += 1 end
                if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
            end
        end
    end
    return (n_nsyn, n_nsyn_restricted)
end;

"""
    get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF)
"""
function get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF)
    n_nsyn, n_nsyn_restricted = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    seq_consensus = copy(csv_index_and_TF.consensus)
    L_TF = length(csv_index_and_TF.HXB2)
    i_eff_restricted_max = count(idx_significant)
    for i_eff in 1:length(csv_raw_in.HXB2_index)
        i_raw = collect(1:L_TF)[csv_index_and_TF.HXB2 .== csv_raw_in.HXB2_index[i_eff]][1] # This line is clitical to get the HXB2 index in TF seq.
        nuc_consensus = seq_consensus[i_raw]
        nuc_MT = csv_raw_in.nucleotide[i_eff]
        nuc_TF = seq_TF[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            if(nuc_MT != nuc_TF)
                seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT        
                flag_nsyn, flag_nsyn_restricted = false, false
                for i_fr in 1:3
                    codon_location = collect(1:3)
                    if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                    if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                    if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end
                    codon_TF = join(seq_TF[codon_location])
                    aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                    codon_MT = join(seq_MT[codon_location])
                    aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                    if((aa_MT != aa_TF) && (nuc_consensus == nuc_MT))  
                        flag_nsyn = true
                        #if(i_eff <= i_eff_restricted_max)
                        if(n_nsyn_restricted <= i_eff_restricted_max)
                            flag_nsyn_restricted = true
                        end
                    end
                end
                if(flag_nsyn) n_nsyn += 1 end
                if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
            end
        end
    end 
    return (n_nsyn, n_nsyn_restricted)
end;

function filter_nuc_mut(mutant_types_set_nuc_simple)
    mutant_nuc_simple_filtered = []    
    for i in 1:length(mutant_types_set_nuc_simple[2])
        x1 = mutant_types_set_nuc_simple[1][i]
        x2 = mutant_types_set_nuc_simple[2][i]
        x3 = mutant_types_set_nuc_simple[3][i]
        temp_set = []
        for x in [x1,x2,x3]
            if(x != "NA" && x != "")
                push!(temp_set, x)
            end
        end
        if(length(temp_set)==0) push!(temp_set, "NA") end 

        push!(mutant_nuc_simple_filtered, unique(temp_set)[1])
    end;
    #@assert count(length.(mutant_nuc_simple_filtered) .!= 1) == 0;
    return mutant_nuc_simple_filtered
end;
