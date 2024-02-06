idx_HXB2_V1 = collect(6615:6693)
idx_HXB2_V2 = collect(6693:6812)
idx_HXB2_V3 = collect(7110:7217)
idx_HXB2_V4 = collect(7377:7478)
idx_HXB2_V5 = collect(7602:7634)
idx_HXB2_LD = collect(7047:7073)
idx_HXB2_MPER = collect(8202:8273)
idx_HXB2_CD4BS = [
    collect(6594:6605); 
    collect(6695:6697);
    collect(6810:6812);
    collect(6816:6818);
    collect(7059:7073);
    collect(7314:7334);
    collect(7344:7346);
    collect(7497:7520);
    collect(7587:7607);
    collect(7629:7631); 
    collect(7635:7655) ]

# ------- Protein Sequence ------- #
idx_HXB2_Pro_V1 = collect(131:156)
idx_HXB2_Pro_V2 = collect(157:196)
idx_HXB2_Pro_V3 = collect(296:331)
idx_HXB2_Pro_V4 = collect(385:418)
idx_HXB2_Pro_V5 = collect(460:470)
idx_HXB2_Pro_LD = collect(275:283)
idx_HXB2_Pro_MPER = collect(660:683)
idx_HXB2_Pro_CD4BS = [
    collect(124:127); 
    collect(157:158);
    collect(196:196);
    collect(198:198);
    collect(279:283);
    collect(364:370);
    collect(374:374);
    collect(425:432);
    collect(455:461);
    collect(469:469); 
    collect(471:477) ]


# --------------------------------- CH848 Specific --------------------------------- #
idx_HXB2_entire_gene = collect(1:9719)
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

""" This function returns the presence and absence of glycan, i.e., glycan shield and glycan hole. 
    Consider the full genetic background comparing the set of consective 5 codons between before and after the mutaions. 
    However, this method consider all possible combination between the before and after the mutations; the possible pair of sequence should be considered to get more accurate glycan dynamics. 
"""
unction get_possible_glycan_naive(i_mut, i_raw, date_mut_set, NUC, a_MT_set, n_poly_idx_max,  collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
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
                    @show aa_set_before
                    @show aa_set_after

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

""" This function returns the presence and absence of glycan, i.e., glycan shield and glycan hole.
    Unlike the previous version, named get_possible_glycan_naive, this function considers the possible pair of sequence to get more accurate glycan dynamics.
    To this end, we import the paired sequences from the mutaion detection function.
"""
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
                    @show aa_set_before
                    @show aa_set_after

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
                
                # --- Get the possible mutation at the site fully considering the genetic backgroud.  --- #
                (date_found, this_gene_set_out, mutation_tot_string_nuc_set_out, mutation_tot_string_nuc_simple_set_out, mutation_tot_string_AA_set_out) = get_possible_mutation_AA_and_NUC(
                    i_mut, i_raw, idx_hxb2, date_mut_set, NUC, a_MT_set, a_WT, n_poly_idx_max, 
                    collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, 
                    seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
                
                # --- Alternative method to get the possible mutation using the only TF's backgrond. This (classical) method does not accurately detect mutations. --- #
                (date_found_naive, this_gene_set_out_naive, mutation_tot_string_nuc_set_out_naive, 
                    mutation_tot_string_nuc_simple_set_out_naive, mutation_tot_string_AA_set_out_naive) = get_possible_mutation_AA_and_NUC_naive(
                    i_mut, i_raw, idx_hxb2, date_mut_set, NUC, a_MT_set, a_WT, n_poly_idx_max, 
                    collected_time, collected_time_unique, idx_poly, idx_only_poly, n_time_max, 
                    seq_TF, csv_index_and_TF, data_num, this_frame_set, this_gene_set)
                @assert date_found == date_found_naive
                @assert this_gene_set_out == this_gene_set_out_naive

                push!(mutant_hxb2, idx_hxb2)
                push!(mutant_nuc, a_mut)
                push!(mutant_date_found, date_found)
                for i_fr in 1:3 
                    #@printf("fr:%d gen:%s mut:%s mut_simple:%s mut_AA:%s\n", 
                    #i_fr, this_gene_set_out[i_fr], mutation_tot_string_nuc_set_out[i_fr], mutation_tot_string_nuc_simple_set_out[i_fr], mutation_tot_string_AA_set_out[i_fr]) 
                    push!(mutant_gene[i_fr], this_gene_set_out[i_fr])
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
            idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw]) 
            this_frame_set, this_gene_set = index2frame(idx_hxb2)
            #a_WT = csv_index_and_TF.TF[i_raw]
            for i_fr in 1:3
               if(i_fr ∈ this_frame_set)                    
                    #frame_temp = this_frame_set[i_fr];
                    codon_location = collect(1:3) # this index should extend 1 to to 2652
                    if(i_raw%3 == (i_fr+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                    if(i_raw%3 == (i_fr+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                    if(i_raw%3 == i_fr%3) codon_location = i_raw .+ collect(-2:1:0) end
                    if(minimum(codon_location)>0 && maximum(codon_location)<=length(seq_TF))
                        x_TF = join(seq_TF[codon_location])
                        AA_TF = haskey(NUC2AA, x_TF) ? NUC2AA[x_TF] : "-" 
                        push!(mutant_types_set_TF_AA[i_fr], AA_TF)
                    else
                        push!(mutant_types_set_TF_AA[i_fr], "NA")
                    end
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
    frame_temp = 3 # consider only 3rd frame where glycan appeares.
    seq_TF = copy(csv_index_and_TF.TF)
    for i_raw in 1:length(seq_TF)
        codon_location_set = [] # that should contains the 5 types of sites. 
        push!(codon_location_set, )
        a_TF = seq_TF[i_raw]
        # Consider only non_gap sites
        if(a_TF != "-")
            codon_location = collect(1:3)
            if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
            if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
            if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end

            if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
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
                        if(length(aa_MT) == length(aa_TF) == 5) # check if the length of aa is 5 to get rid of the edge effect.
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
        
        if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
            [push!(codon_set_tot_before, x) for x in unique([join(data_before_mut_extend[n, codon_location]) for n in 1:n_before])]
            [push!(codon_set_tot_after, x) for x in unique([join(data_after_mut_extend[n, codon_location]) for n in 1:n_after])]
            codon_set_tot_before_original = copy([join(data_before_mut_extend[n, codon_location]) for n in 1:n_before])
            codon_set_tot_after_original = copy([join(data_after_mut_extend[n, codon_location]) for n in 1:n_after])
            codon_set_tot_before= copy(unique(codon_set_tot_before))
            codon_set_tot_after = copy(unique(codon_set_tot_after))
            
            # The following frequency of the codon will be used if there are multiple options. 
            num_of_codon_before = []; num_of_codon_after = []
            [push!(num_of_codon_before, count(x .== codon_set_tot_before_original)) for x in codon_set_tot_before]
            [push!(num_of_codon_after, count(x .== codon_set_tot_after_original)) for x in codon_set_tot_after]
            
            """
            if(extract_integer(idx_hxb2)==7607)
                @printf("cdn_bfr:%s\n", join(codon_set_tot_before, ","))
                @printf("cdn_afr:%s\n", join(codon_set_tot_after, ","))
            end
            """
        
            AA_set_after = [haskey(NUC2AA,x) ? NUC2AA[x] : "-" for x in codon_set_tot_after]
            AA_set_before = [haskey(NUC2AA,x) ? NUC2AA[x] : "-" for x in codon_set_tot_before]
        
            """
            if(extract_integer(idx_hxb2)==7607)
            flag_show = (AA_set_before != AA_set_after) && (i_raw%3 == frame_temp%3 )
            if(flag_show)
                @printf("fr:%d AA_before:%s AA_after:%s \n\n", i_fr, join(AA_set_before, "|"), join(AA_set_after, "|"))
            end
            """
            
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
            #if(extract_integer(idx_hxb2)==7607)
            if(flag_show)
                @printf("TF codon:%s\n", join(seq_TF[codon_location]))
                @printf("Detected Codon: before>%s & after>%s\n", join(codon_set_tot_before, "/"), join(codon_set_tot_after, "/"))
                @printf("Num of Codon: before>%s & after>%s\n", join(string.(num_of_codon_before), ","), join(string.(num_of_codon_after), ","))
                @printf("Detected AA: before>%s & after>%s\n", join(AA_set_before, "/"), join(AA_set_after, "/"))            
                @printf("Possible Transition (row: before, col: after): \n")
                print_matrix(table_paired)
                @printf("Possible Transition (After): \n")
                print_matrix(table_paired_after)
            end
            """
            # ------------------------------- Adjoint paired/unpaired mutations ----------------------------------#
            mutation_tot_string_AA = ""
            mutation_tot_string_nuc = ""
            mutation_tot_string_nuc_simple = join(unique(data_before_mut_extend[:, i_raw]), "/") * idx_hxb2 * a_mut
            if(size(codon_set_tot_after, 1)>0 || size(codon_set_tot_before,1)>0)
                # --- Making the paires that has single nucleotide while considerng the genetic background ----        
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
                                codon_set_tot_before_rem_gap_copy = copy(codon_set_tot_before_rem_gap)
                                codon_set_tot_before_rem_gap_copy[codon_set_tot_before_rem_gap_copy .== ""] .= "---"
                                AA_set_before_rem_gap_copy = copy(AA_set_before_rem_gap)
                                AA_set_before_rem_gap_copy[AA_set_before_rem_gap_copy .== ""] .= "-"
                                before_str_nuc = join(codon_set_tot_before_rem_gap_copy, "/")
                                before_str_AA = join(unique(AA_set_before_rem_gap_copy), "/")
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
            #if(extract_integer(idx_hxb2)==7607)
            if(flag_show)
                @printf("fr:%d gene:%s poly_idx:%d HXB2_idx:%s nuc:%s detected:%dd\n", frame_temp, this_gene_set_out[frame_temp], i_raw, idx_hxb2, a_mut, date_found)
                @printf("GENE:%s AA:%s NUC:%s\n\n", mutation_tot_string_nuc_simple, mutation_tot_string_AA, mutation_tot_string_nuc)
            end
            #"""
            end # end of frame.

        end
    
    return (date_found, this_gene_set_out, mutation_tot_string_nuc_set_out, 
        mutation_tot_string_nuc_simple_set_out, mutation_tot_string_AA_set_out) 
end;

""" This function returns the possible mutations comparing with TF, there is no information of backgroud sequence. 
"""
function get_possible_mutation_AA_and_NUC_naive(i_mut, i_raw, idx_hxb2, date_mut_set, NUC, a_MT_set, a_WT, n_poly_idx_max, 
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
    
    n_after = size(data_after_mut,1)
    data_after_mut_extend = [x for i in 1:n_after , x in seq_TF]
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
        
        codon_set_tot_after = []
        codon_location = collect(1:3) # this index should extend 1 to to 2652
        if(i_raw%3 == (frame_temp+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
        if(i_raw%3 == (frame_temp+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
        if(i_raw%3 == frame_temp%3) codon_location = i_raw .+ collect(-2:1:0) end

        if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
            codon_tf = seq_TF[codon_location]
            AA_tf = haskey(NUC2AA, join(codon_tf)) ? NUC2AA[join(codon_tf)] : "-"
            [push!(codon_set_tot_after, x) for x in unique([join(data_after_mut_extend[n, codon_location]) for n in 1:n_after])]
            codon_set_tot_after = copy(unique(codon_set_tot_after))
            len_codon_after = length(codon_set_tot_after)
            AA_set_after = [haskey(NUC2AA,x) ? NUC2AA[x] : "-" for x in codon_set_tot_after]
        
            nuc_codon_idex = csv_index_and_TF.HXB2[codon_location]
            nuc_mut_string = join([nuc_codon_idex[1], nuc_codon_idex[3]], "-")        
            mutation_tot_string_nuc_set_out[frame_temp] = join(codon_tf) * nuc_mut_string * join(codon_set_tot_after, "/") 
            mutation_tot_string_nuc_simple_set_out[frame_temp] = a_WT * idx_hxb2 * a_mut
            mutation_tot_string_AA_set_out[frame_temp] = AA_tf * i_AA * join(AA_set_after, "/") 
        end
    end # end of frame.

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
    #if(α_null > 0) x_fold = α_sel / α_null end
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

""" get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_type)

"""
function get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_type)
    n_nsyn_restricted, n_nsyn = 0, 0
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
                        if(i_fr ∈ this_frame_set)
                            codon_location = collect(1:3)
                            if(i_raw%3 == (i_fr+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                            if(i_raw%3 == (i_fr+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                            if(i_raw%3 == i_fr%3) codon_location = i_raw .+ collect(-2:1:0) end
                            if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
                                codon_TF = join(seq_TF[codon_location])
                                aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                                codon_MT = join(seq_MT[codon_location])
                                aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                                if(aa_MT != aa_TF) 
                                    flag_nsyn = true
                                    if(idx_hxb2 ∈ idx_HXB2_type)
                                        flag_nsyn_restricted = true
                                    end
                                end
                            end
                        end
                    end
                    if(flag_nsyn) n_nsyn += 1 end
                    if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
                end
            end
        end
    end
    return (n_nsyn_restricted, n_nsyn)
end

# idx_type is not necessary for this function. --> need to fix this!
function get_num_of_nonsyn_reversion(csv_index_and_TF)
    n_nsyn_restricted, n_nsyn = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    seq_consensus = copy(csv_index_and_TF.consensus)
    for i_raw in 1:length(seq_TF)
        nuc_TF = seq_TF[i_raw]
        nuc_consensus = seq_consensus[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            for nuc_MT in NUC
                if(nuc_MT != nuc_TF && nuc_MT != "-") 
                    seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT
                    flag_nsyn_restricted, flag_nsyn = false, false
                    for i_fr in 1:3
                        if(i_fr ∈ this_frame_set)
                            codon_location = collect(1:3)
                            if(i_raw%3 == (i_fr+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                            if(i_raw%3 == (i_fr+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                            if(i_raw%3 == i_fr%3) codon_location = i_raw .+ collect(-2:1:0) end
                            if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
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
                        end
                    end
                    if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
                    if(flag_nsyn) n_nsyn += 1 end
                end
            end
        end
    end 
    return (n_nsyn_restricted, n_nsyn)
end;

function get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in, idx_significant)
    # The following computation is similar to the n_null and N_null computation
    n_nsyn_restricted = 0
    seq_TF = copy(csv_index_and_TF.TF);
    L_TF = length(csv_index_and_TF.HXB2)
    i_eff_restricted_max = count(idx_significant)
    for i_eff in 1:count(idx_significant)
        i_raw = collect(1:L_TF)[csv_index_and_TF.HXB2 .== csv_raw_in.HXB2_index[i_eff]][1] # This line is clitical to get the HXB2 index in TF seq.
        nuc_MT = csv_raw_in.nucleotide[i_eff]
        nuc_TF = seq_TF[i_raw]
        idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            if(nuc_MT != nuc_TF)
                seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT        
                flag_nsyn_restricted = false
                for i_fr in 1:3
                    if(i_fr ∈ this_frame_set)
                        codon_location = collect(1:3)
                        if(i_raw%3 == (i_fr+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                        if(i_raw%3 == (i_fr+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                        if(i_raw%3 == i_fr%3) codon_location = i_raw .+ collect(-2:1:0) end
                        if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
                            codon_TF = join(seq_TF[codon_location])
                            aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                            codon_MT = join(seq_MT[codon_location])
                            aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                            if(aa_MT != aa_TF) 
                                #if(idx_hxb2 ∈ idx_type_in)
                                if(idx_type_in[i_eff])
                                    flag_nsyn_restricted = true
                                end
                            end
                        end
                    end
                end
                if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
            end
        end
    end
    return n_nsyn_restricted
end;

""" get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF)
"""
function get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF, idx_significant)
    n_nsyn_restricted, n_nsyn = 0, 0
    seq_TF = copy(csv_index_and_TF.TF);
    seq_consensus = copy(csv_index_and_TF.consensus)
    L_TF = length(csv_index_and_TF.HXB2)
    i_eff_restricted_max = count(idx_significant)
    for i_eff in 1:count(idx_significant)
        i_raw = collect(1:L_TF)[csv_index_and_TF.HXB2 .== csv_raw_in.HXB2_index[i_eff]][1] # This line is clitical to get the HXB2 index in TF seq.
        nuc_consensus = seq_consensus[i_raw]
        nuc_MT = csv_raw_in.nucleotide[i_eff]
        nuc_TF = seq_TF[i_raw]
        #idx_hxb2 = extract_integer(csv_index_and_TF.HXB2[i_raw])
        idx_hxb2 = extract_integer(csv_raw_in.HXB2_index[i_eff])
        this_frame_set, this_gene_set = index2frame(idx_hxb2)
        if(nuc_TF != "-")
            if(nuc_MT != nuc_TF)
                seq_MT = copy(seq_TF); seq_MT[i_raw] = nuc_MT        
                flag_nsyn_restricted, flag_nsyn = false, false
                for i_fr in 1:3
                    if(i_fr ∈ this_frame_set)
                        codon_location = collect(1:3)
                        if(i_raw%3 == (i_fr+1)%3) codon_location = i_raw .+ collect( 0:1:2) end
                        if(i_raw%3 == (i_fr+2)%3) codon_location = i_raw .+ collect(-1:1:1) end
                        if(i_raw%3 == i_fr%3) codon_location = i_raw .+ collect(-2:1:0) end
                        if(0 < codon_location[1] && codon_location[end] <= length(seq_TF)) 
                            codon_TF = join(seq_TF[codon_location])
                            aa_TF = haskey(NUC2AA, codon_TF) ? NUC2AA[codon_TF] : "-"
                            codon_MT = join(seq_MT[codon_location])
                            aa_MT = haskey(NUC2AA, codon_MT) ? NUC2AA[codon_MT] : "-"
                            if(aa_MT != aa_TF)  
                                flag_nsyn = true
                                #if(nuc_consensus == nuc_MT)
                                if(csv_raw_in.reversion[i_eff])
                                    flag_nsyn_restricted = true
                                end
                            end
                        end
                    end
                end
                if(flag_nsyn_restricted) n_nsyn_restricted += 1 end
                if(flag_nsyn) n_nsyn += 1 end
            end
        end
    end 
    return (n_nsyn_restricted, n_nsyn) 
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

""" log_factorial_stirling(n)
"""
function log_factorial_stirling(n)
    if n <=0
        return 0.0
    else
        return (n * log(n) - n + log(2π * n) / 2) / log(10)  # Stiring's approximation.
    end
end

""" log_binomial_coefficient(n, k)
"""
function log_binomial_coefficient(n, k)
    return log_factorial_stirling(n) - log_factorial_stirling(k) - log_factorial_stirling(n - k)
end

""" Fishers_exact_test(A, B, C, D)
"""
function fishers_exact_test(A, B, C, D)
    total = A + B + C + D
    #@printf("A+B:%d\tA:%d\tC+D:%d\tC:%d\ttotal:%d\tA+D:%d\n", A, B, C, D, total, A+D)
    log_p_value = log_binomial_coefficient(A + B, A) + 
                  log_binomial_coefficient(C + D, C) - 
                  log_binomial_coefficient(total, A + C)
    return log_p_value
end;

""" get_enrichment_and_pvalues(csv_raw_in, csv_index_and_TF, s_threshold_set)
    # This function is used for the enrichment calculation and Fisher's exact test. 
    # The input is the csv file of the mutations and the csv file of the TF seq.
    # The output is the enrichment and the p-value for each mutation.
"""
function get_enrichment_and_pvalues(csv_raw_in, csv_index_and_TF, α_selected_list)
    x_fold_V1_set = []; n_fold_V1_set = []; fisher_V1_set = []
    x_fold_V2_set = []; n_fold_V2_set = []; fisher_V2_set = []
    x_fold_V3_set = []; n_fold_V3_set = []; fisher_V3_set = []
    x_fold_V4_set = []; n_fold_V4_set = []; fisher_V4_set = []
    x_fold_V5_set = []; n_fold_V5_set = []; fisher_V5_set = []
    x_fold_reversion_set = []; n_fold_reversion_set = []; fisher_reversion_set = []
    x_fold_LoopD_set = []; n_fold_LoopD_set = []; fisher_LoopD_set = []
    x_fold_CD4BS_set = []; n_fold_CD4BS_set = []; fisher_CD4BS_set = []
    x_fold_N_any_set = []; n_fold_N_any_set = []; fisher_N_any_set = []
    x_fold_N_rem_set = []; n_fold_N_rem_set = []; fisher_N_rem_set = []
    x_fold_N_add_set = []; n_fold_N_add_set = []; fisher_N_add_set = []
    x_fold_N_sht_set = []; n_fold_N_sht_set = []; fisher_N_sht_set = []
    α_selected_est = []
    # ------- Make boolean vectors -------- #
    idx_type_in_all = [true for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_LD = [x ∈ idx_HXB2_LD for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_CD4BS = [x ∈ idx_HXB2_CD4BS for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_V1 = [x ∈ idx_HXB2_V1 for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_V2 = [x ∈ idx_HXB2_V2 for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_V3 = [x ∈ idx_HXB2_V3 for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_V4 = [x ∈ idx_HXB2_V4 for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_V5 = [x ∈ idx_HXB2_V5 for x in extract_integer.(csv_raw_in.HXB2_index)] 
    idx_type_in_N_any = copy(csv_raw_in.N_linked_glycan_plus_fr3 .> 0 .|| csv_raw_in.N_linked_glycan_minus_fr3 .> 0 )
    idx_type_in_N_add = copy(csv_raw_in.N_linked_glycan_plus_fr3 .> 0 )
    idx_type_in_N_rem = copy(csv_raw_in.N_linked_glycan_minus_fr3 .> 0 )
    idx_type_in_N_sht = copy(csv_raw_in.N_linked_glycan_shift_fr3 .> 0 )

    # ------- Get the nuber of mutations assuming the random prediction (= null moodel) -------- #
    idx_all_true = repeat([true], length(csv_raw_in.HXB2_index))
    N_sel_all = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_all, idx_all_true)
    #
    n_sel_all_LD = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_LD, idx_all_true) # only for Fisher's exact test!
    n_sel_all_CD4BS = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_CD4BS, idx_all_true) # only for Fisher's exact test!
    n_sel_all_V1 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V1, idx_all_true) # only for Fisher's exact test!
    n_sel_all_V2 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V2, idx_all_true) # only for Fisher's exact test!
    n_sel_all_V3 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V3, idx_all_true) # only for Fisher's exact test!
    n_sel_all_V4 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V4, idx_all_true) # only for Fisher's exact test!
    n_sel_all_V5 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V5, idx_all_true) # only for Fisher's exact test!
    n_sel_all_N_any = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_any, idx_all_true) # only for Fisher's exact test!
    n_sel_all_N_add = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_add, idx_all_true) # only for Fisher's exact test!
    n_sel_all_N_rem = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_rem, idx_all_true) # only for Fisher's exact test!
    n_sel_all_N_sht = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_sht, idx_all_true) # only for Fisher's exact test!
    (n_sel_all_rev, n_sel_all_rev_temp) = get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF, idx_all_true) # only for Fisher's exact test!

    # ------- Get the nuber of mutations assuming the random prediction (= null moodel) -------- #
    (N_null_temp_all, N_null_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_entire_gene)
    @assert N_null_temp_all == N_null_all
    (n_null_LD, N_null_LD_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_LD)
    (n_null_CD4BS, N_null_CD4BS_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_CD4BS)
    (n_null_V1, N_null_V1_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_V1)
    (n_null_V2, N_null_V2_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_V2)
    (n_null_V3, N_null_V3_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_V3)
    (n_null_V4, N_null_V4_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_V4)
    (n_null_V5, N_null_V5_all) = get_num_of_nonsyn(csv_index_and_TF, idx_HXB2_V5)
    (N_null_N_glycan_all, n_null_N_add, n_null_N_rem, n_null_N_sht) = get_glycan_plus_minus_shift_statistics(csv_index_and_TF);
    n_null_N_any = n_null_N_add + n_null_N_rem - n_null_N_sht
    (n_null_rev, N_null_rev_all) = get_num_of_nonsyn_reversion(csv_index_and_TF);

    #@printf("%d %d %d %d %d %d %d %d \n", N_null_all, N_null_LD_all, N_null_CD4BS_all, N_null_V1_all, N_null_V2_all, N_null_V3_all, N_null_V4_all, N_null_V5_all)
    #@printf("%d %d\n", N_null_N_glycan_all, N_null_rev_all)
    #@printf("%d %d %d\n\n", n_null_N_add, n_null_N_rem, n_null_N_sht);
    #@printf("%d %d %d %d %d %d %d \n", n_sel_all_LD, n_sel_all_CD4BS, n_sel_all_V1, n_sel_all_V2, n_sel_all_V3, n_sel_all_V4, n_sel_all_V5)
    #@printf("%d, %d %d %d \n", n_sel_all_rev, n_sel_all_N_add, n_sel_all_N_rem, n_sel_all_N_sht);

    # ------- Get the nuber of mutations -------- #
    # ( top x=100% && property)
    #for s_threshold_RMs in s_threshold_set
    for α_selected in α_selected_list 
    #s_threshold_RMs = 0.02
        len_tot_entries = length(csv_raw_in.HXB2_index)
	num_selected = Int( floor(len_tot_entries * α_selected) )
  		
        idx_significant = csv_raw_in.s_MPL .>= csv_raw_in.s_MPL[num_selected];
        #len_selected = count(idx_significant)
        #α_selected = len_selected / len_tot_entries;
        n_mut_tot = length(idx_significant) # For Fisher's exact test
        n_mut_sel_tot = count(idx_significant) # For Fisher's exact test

        n_sel_LD = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_LD, idx_significant)
        (n_subjected_LoopD, n_estimated_LoopD, n_tot_LoopD, x_fold_LoopD) = get_x_fold(n_sel_LD, α_selected * N_sel_all, n_null_LD, N_null_LD_all)
        # Fisher's exact test
        A = n_sel_LD; B = n_mut_sel_tot - A; C = n_null_LD - A
        D = N_null_all - (A+B+C)
        log_Fisher_LD = fishers_exact_test(A, B, C, D)

        n_sel_CD4BS = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_CD4BS, idx_significant)
        (n_subjected_CD4BS, n_estimated_CD4BS, n_tot_CD4BS, x_fold_CD4BS) = get_x_fold(n_sel_CD4BS, α_selected * N_sel_all, n_null_CD4BS, N_null_CD4BS_all)
        # Fisher's exact test
        A = n_sel_CD4BS; B = n_mut_sel_tot - A; C = n_null_CD4BS - A
        D = N_null_all - (A+B+C)
        log_Fisher_CD4BS = fishers_exact_test(A, B, C, D)

        n_sel_V1 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V1, idx_significant)
        (n_subjected_V1, n_estimated_V1, n_tot_V1, x_fold_V1) = get_x_fold(n_sel_V1,α_selected * N_sel_all, n_null_V1, N_null_V1_all)
        # Fisher's exact test
        A = n_sel_V1; B = n_mut_sel_tot - A; C = n_null_V1 - A
        D = N_null_all - (A+B+C)
        log_Fisher_V1 = fishers_exact_test(A, B, C, D)

        n_sel_V2 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V2, idx_significant)
        (n_subjected_V2, n_estimatd_V2, n_tot_V2, x_fold_V2) = get_x_fold(n_sel_V2, α_selected * N_sel_all, n_null_V2, N_null_V2_all)
        # Fisher's exact test
        A = n_sel_V2; B = n_mut_sel_tot - A; C = n_null_V2 - A
        D = N_null_all - (A+B+C)
        log_Fisher_V2 = fishers_exact_test(A, B, C, D)

        n_sel_V3 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V3, idx_significant)
        (n_subjected_V3, n_estimated_V3, n_tot_V3, x_fold_V3) = get_x_fold(n_sel_V3, α_selected * N_sel_all, n_null_V3, N_null_V3_all)
        # Fisher's exact test
        A = n_sel_V3; B = n_mut_sel_tot - A; C = n_null_V3 - A
        D = N_null_all - (A+B+C)
        log_Fisher_V3 = fishers_exact_test(A, B, C, D)


        n_sel_V4 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V4, idx_significant)
        (n_subjected_V4, n_estimated_V4, n_tot_V4, x_fold_V4) = get_x_fold(n_sel_V4, α_selected * N_sel_all, n_null_V4, N_null_V4_all)
        # Fisher's exact test
        A = n_sel_V4; B = n_mut_sel_tot - A; C = n_null_V4 - A
        D = N_null_all - (A+B+C)
        log_Fisher_V4 = fishers_exact_test(A, B, C, D)

        n_sel_V5 = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_V5, idx_significant)
        (n_subjected_V5, n_estimated_V5, n_tot_V5, x_fold_V5) = get_x_fold(n_sel_V5, α_selected * N_sel_all, n_null_V5, N_null_V5_all)    
        # Fisher's exact test
        A = n_sel_V5; B = n_mut_sel_tot - A; C = n_null_V5 - A
        D = N_null_all - (A+B+C)
        log_Fisher_V5 = fishers_exact_test(A, B, C, D)

        #  ------------ Reversion Mutations -------------  #
        (n_sel_rev, N_sel_rev) = get_n_sel_and_N_sel_reversion(csv_raw_in, csv_index_and_TF, idx_significant)
        (n_subjected_reversion, n_estimated_reversion, n_tot_reversion, x_fold_reversion) = get_x_fold(n_sel_rev, α_selected * N_sel_all, n_null_rev, N_null_rev_all)

        # Fisher's exact test
        A = n_sel_rev; B = n_mut_sel_tot - A; C = n_null_rev - A
        D = N_null_rev_all - (A+B+C)
        log_Fisher_rev = fishers_exact_test(A, B, C, D)

        #  ------------ Glycan Mutations -------------  #
        n_sel_N_add = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_add, idx_significant)
        n_sel_N_rem = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_rem, idx_significant)
        n_sel_N_sht = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_sht, idx_significant)
        n_sel_N_any = get_n_sel_and_N_sel(csv_raw_in, csv_index_and_TF, idx_type_in_N_any, idx_significant)
        #
        (n_subjected_N_add, n_estimated_N_add, n_tot_N_add, x_fold_N_add) = get_x_fold(n_sel_N_add, α_selected * N_sel_all, n_null_N_add, N_null_N_glycan_all)
        (n_subjected_N_rem, n_estimated_N_rem, n_tot_N_rem, x_fold_N_rem) = get_x_fold(n_sel_N_rem, α_selected * N_sel_all, n_null_N_rem, N_null_N_glycan_all)
        (n_subjected_N_sht, n_estimated_N_sht, n_tot_N_sht, x_fold_N_sht) = get_x_fold(n_sel_N_sht, α_selected * N_sel_all, n_null_N_sht, N_null_N_glycan_all)
        (n_subjected_N_any, n_estimated_N_any, n_tot_N_any, x_fold_N_any) = get_x_fold(n_sel_N_any, α_selected * N_sel_all, n_null_N_any, N_null_N_glycan_all)
        # Fisher's exact test
        A = n_sel_N_add; B = n_mut_sel_tot - A; C = n_null_N_add - A
        D = N_null_N_glycan_all - (A+B+C)
        log_Fisher_N_add = fishers_exact_test(A, B, C, D)
        # Fisher's exact test
        A = n_sel_N_rem; B = n_mut_sel_tot - A; C = n_null_N_rem - A
        D = N_null_N_glycan_all - (A+B+C)
        log_Fisher_N_rem = fishers_exact_test(A, B, C, D)
        # Fisher's exact test
        A = n_sel_N_sht; B = n_mut_sel_tot - A; C = n_null_N_sht - A
        D = N_null_N_glycan_all - (A+B+C)
        log_Fisher_N_sht = fishers_exact_test(A, B, C, D)
        # Fisher's exact test
        A = n_sel_N_any; B = n_mut_sel_tot - A; C = n_null_N_any - A
        D = N_null_N_glycan_all - (A+B+C)
        log_Fisher_N_any = fishers_exact_test(A, B, C, D)
        # --------------------------------------------------------------- #     
        #@printf("α:%.2f: V1:%.1f(%d,%d) V2:%.1f(%d,%d) V3:%.1f(%d,%d) V4:%.1f(%d,%d) V5:%.1f(%d,%d)\n", 
        #        100*α_selected, x_fold_V1, n_sel_V1, n_null_V1, x_fold_V2, n_sel_V2, n_null_V2, x_fold_V3, n_sel_V3, n_null_V3, x_fold_V4, n_sel_V4, n_null_V4, x_fold_V5, n_sel_V5, n_null_V5)

        push!(x_fold_LoopD_set, x_fold_LoopD); push!(n_fold_LoopD_set, n_subjected_LoopD); push!(fisher_LoopD_set, log_Fisher_LD)
        push!(x_fold_CD4BS_set, x_fold_CD4BS); push!(n_fold_CD4BS_set, n_subjected_CD4BS); push!(fisher_CD4BS_set, log_Fisher_CD4BS)
        push!(x_fold_V1_set, x_fold_V1); push!(n_fold_V1_set, n_subjected_V1); push!(fisher_V1_set, log_Fisher_V1)
        push!(x_fold_V2_set, x_fold_V2); push!(n_fold_V2_set, n_subjected_V2); push!(fisher_V2_set, log_Fisher_V2)
        push!(x_fold_V3_set, x_fold_V3); push!(n_fold_V3_set, n_subjected_V3);  push!(fisher_V3_set, log_Fisher_V3)
        push!(x_fold_V4_set, x_fold_V4); push!(n_fold_V4_set, n_subjected_V4); push!(fisher_V4_set, log_Fisher_V4)
        push!(x_fold_V5_set, x_fold_V5); push!(n_fold_V5_set, n_subjected_V5); push!(fisher_V5_set, log_Fisher_V5)

        push!(x_fold_reversion_set, x_fold_reversion); push!(n_fold_reversion_set, n_subjected_reversion); push!(fisher_reversion_set, log_Fisher_rev)
        push!(x_fold_N_rem_set, x_fold_N_rem); push!(n_fold_N_rem_set, n_subjected_N_rem); push!(fisher_N_rem_set, log_Fisher_N_rem)
        push!(x_fold_N_add_set, x_fold_N_add); push!(n_fold_N_add_set, n_subjected_N_add); push!(fisher_N_add_set, log_Fisher_N_add)
        push!(x_fold_N_sht_set, x_fold_N_sht); push!(n_fold_N_sht_set, n_subjected_N_sht); push!(fisher_N_sht_set, log_Fisher_N_sht)
        push!(x_fold_N_any_set, x_fold_N_any); push!(n_fold_N_any_set, n_subjected_N_any); push!(fisher_N_any_set, log_Fisher_N_any)
        push!(α_selected_est, α_selected)
    end


    x_fold_summary = [x_fold_LoopD_set; x_fold_CD4BS_set; x_fold_V1_set; x_fold_V2_set; x_fold_V3_set; x_fold_V4_set; 
     x_fold_V5_set; x_fold_reversion_set; x_fold_N_rem_set; x_fold_N_add_set; x_fold_N_sht_set; x_fold_N_any_set]
    #
    log_P = [fisher_LoopD_set; fisher_CD4BS_set; fisher_V1_set; fisher_V2_set; fisher_V3_set; fisher_V4_set; 
     fisher_V5_set; fisher_reversion_set; fisher_N_rem_set; fisher_N_add_set; fisher_N_sht_set; fisher_N_any_set]
    #
    α_summary = [α_selected_est; α_selected_est; α_selected_est; α_selected_est; α_selected_est; α_selected_est; 
     α_selected_est; α_selected_est; α_selected_est; α_selected_est; α_selected_est; α_selected_est]
    #
    n_fold_summary = [n_fold_LoopD_set; n_fold_CD4BS_set; n_fold_V1_set; n_fold_V2_set; n_fold_V3_set; n_fold_V4_set; 
     n_fold_V5_set; n_fold_reversion_set; n_fold_N_rem_set; n_fold_N_add_set; n_fold_N_sht_set; n_fold_N_any_set]
    #    
    types_summary = [["LD" for _ in x_fold_LoopD_set]; ["CD4BS" for _ in x_fold_CD4BS_set]; 
    ["V1" for _ in x_fold_V1_set]; ["V2" for _ in x_fold_V2_set]; ["V3" for _ in x_fold_V3_set]; ["V4" for _ in x_fold_V4_set]; 
    ["V5" for _ in x_fold_V5_set]; ["Reversion" for _ in x_fold_reversion_set]; 
    ["Hole" for _ in x_fold_N_rem_set]; ["Shield" for _ in x_fold_N_add_set]; ["Shift" for _ in x_fold_N_sht_set]; ["PNG" for _ in x_fold_N_any_set]];

    return (x_fold_summary, log_P, α_summary, n_fold_summary, types_summary)
    
end;


function normalize_coefficient_AA(coefficients_MPL, coefficients_SL, q, L_poly, seq_TF, seq_ensemble, seq_TF_num)

    s_normed_MPL = copy(coefficients_MPL)
    s_normed_SL = copy(coefficients_SL);

    for i in 1:L_poly
        a = seq_TF_num[i]
        s_normed_MPL[km.(i, 1:q, q)] .-= coefficients_MPL[km(i,a,q)]
        s_normed_SL[km.(i, 1:q, q)] .-= coefficients_SL[km(i,a,q)]
        set_aa_observed = unique(seq_ensemble[:, i])
        for b in 1:q
            if(b ∉ set_aa_observed)
                s_normed_MPL[km(i,b,q)] = 0
                s_normed_SL[km(i,b,q)] = 0
            end
        end
    end
    return (s_normed_MPL, s_normed_SL)
end;

function normalize_coefficient_AA(coefficients_MPL, coefficients_SL, q, L_poly, seq_TF, seq_TF_num)
    s_normed_MPL = copy(coefficients_MPL)
    s_normed_SL = copy(coefficients_SL);
    for i in 1:L_poly
        a = seq_TF_num[i]
        s_normed_MPL[km(i, a, q)] -= coefficients_MPL[km(i,a,q)]
        s_normed_SL[km(i, a, q)] -= coefficients_SL[km(i,a,q)]
        for b in 1:q
            if(s_normed_MPL[km(i, b, q)] != 0)
                s_normed_MPL[km(i, b, q)] -= coefficients_MPL[km(i,a,q)]
            end
            if(s_normed_SL[km(i, b, q)] != 0)
                s_normed_SL[km(i, b, q)] -= coefficients_SL[km(i,a,q)]
            end
        end
    end
    return (s_normed_MPL, s_normed_SL)
end;

function get_x1_AA(time_set, L_poly, q, seq_ensemble)
    time_unique = sort(unique(time_set))
    n_time_unique = length(time_unique)
    x1 = zeros(n_time_unique, L_poly, q);

    for i_t in 1:length(time_unique)
        t = time_unique[i_t];
        n_seq_at_t = count(time_set .== t)
        seq_at_t = seq_ensemble[time_set .== t, :];
        for n in 1:n_seq_at_t
            for i in 1:L_poly
                x1[i_t, i, seq_at_t[n, i]] += (1.0 / n_seq_at_t)
            end
        end
    end
    return x1
end

function get_first_detected(time_unique, x1, i_poly, a_poly_num)
    for i_t in 1:length(time_unique)
        if(x1[i_t, i_poly, a_poly_num]>0)
            return time_unique[i_t]
        end
    end
end;

# Input hxb2_out_set, aa_out_set, seq_TF_aa, haxb2_TF

function get_glycan_plus_minus_AA_seq(seq_TF_aa, haxb2_TF, hxb2_out_set, aa_out_set)
    glycan_plus = []; glycan_minus = []
    L_aa_temp = length(seq_TF_aa)
    for n in 1:length(hxb2_out_set)
        hxb2_temp, aa_temp = hxb2_out_set[n], aa_out_set[n]
        flag_glycan_plus, flag_glycan_minus = false, false
        idx_matched = haxb2_TF .== hxb2_temp
        if(count(idx_matched)>0)
            site_matched = collect(1:L_aa_temp)[idx_matched][1]
            if( 3<= site_matched <= L_aa_temp-2)
                seq_WT_cliped = copy(seq_TF_aa[(site_matched-2):(site_matched+2)])                
                seq_MT_cliped = copy(seq_WT_cliped)                
                seq_MT_cliped[3] = aa_temp
                (flag_glycan_plus,flag_glycan_minus) = check_glycan_shield_hole_shift(seq_MT_cliped, seq_WT_cliped)
            end
        end
        push!(glycan_plus, flag_glycan_plus); push!(glycan_minus, flag_glycan_minus)
    end
    return (glycan_plus, glycan_minus)
end;

function get_variable_site_true_false(haxb2_TF)
    v1_set_out, v2_set_out, v3_set_out = [], [], []
    v4_set_out, v5_set_out, LD_set_out, CD4BS_set_out = [], [], [], []
    for x in extract_integer.(haxb2_TF)
        push!(v1_set_out, x ∈ idx_HXB2_Pro_V1)
        push!(v2_set_out, x ∈ idx_HXB2_Pro_V2)
        push!(v3_set_out, x ∈ idx_HXB2_Pro_V3)
        push!(v4_set_out, x ∈ idx_HXB2_Pro_V4)
        push!(v5_set_out, x ∈ idx_HXB2_Pro_V5)
        push!(LD_set_out, x ∈ idx_HXB2_Pro_LD)
        push!(CD4BS_set_out, x ∈ idx_HXB2_Pro_CD4BS)
    end;
    return (v1_set_out, v2_set_out, v3_set_out, v4_set_out, v5_set_out, LD_set_out, CD4BS_set_out)
end;
