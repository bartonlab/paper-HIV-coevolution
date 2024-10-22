# --- this should be able to give HXB2 sequence automatically.
# -------- need to tell this function the begining and end HXB2 index. --------.
PRO = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"];

function index2frame(i)
    """ Return the open reading frames corresponding to a given HXB2 index. """

    frames = Int[]
    genes = String[]

    if (790 ≤ i ≤ 2292) || (5041 ≤ i ≤ 5619) || (8379 ≤ i ≤ 8469) || (8797 ≤ i ≤ 9417)
        push!(frames, 1)
        if(790 ≤ i ≤ 2292) push!(genes, "gag") end
        if(5041 ≤ i ≤ 5619) push!(genes, "vif") end
        if(8379 ≤ i ≤ 8469) push!(genes, "tat") end
        if(8797 ≤ i ≤ 9417) push!(genes, "nef") end            
    end
    if (5831 ≤ i ≤ 6045) || (6062 ≤ i ≤ 6310) || (8379 ≤ i ≤ 8653)
        push!(frames, 2)
        if(5831 ≤ i ≤ 6045) push!(genes, "tat") end
        if(6062 ≤ i ≤ 6310) push!(genes, "vpu") end
        if(8379 ≤ i ≤ 8653) push!(genes, "rev") end
    end
    if (2253 ≤ i ≤ 5096) || (5559 ≤ i ≤ 5850) || (5970 ≤ i ≤ 6045) || (6225 ≤ i ≤ 8795)
        push!(frames, 3)
        if(2253 ≤ i ≤ 5096) push!(genes, "pol") end
        if(5559 ≤ i ≤ 5850) push!(genes, "vpr") end
        if(5970 ≤ i ≤ 6045) push!(genes, "vpr") end
        if(6225 ≤ i ≤ 8795) push!(genes, "env") end
    end

    return frames, genes
end;

""" 
clip_index_start and clip_index_end are HXB2 indices that we want to ensure the output sequence belone in between them.
seq_index_start and seq_index_end are HXB2 indices that correspond to the begining and end of the current input HXB2 sequence.
"""
function get_index_set_clipping_MSA(MSA, clip_index_start, clip_index_end, seq_index_start, seq_index_end)
    INDEX_HXB2 = collect(1:size(seq_headder,1))[occursin.("HXB2", seq_headder)];
    if(length(INDEX_HXB2)==0)
        @printf "Error: This MSA does not contain HXB2 sequence!"
    end
    
    seq_HXB2 = MSA[INDEX_HXB2][1]
    size(seq_HXB2,1)
    indices_clipping = zeros(Int, size(seq_HXB2,1))
    index_HXB2 = seq_index_start - 1
    # clipping the sequences from [seq_index_start:seq_index_end] to [clip_index_start:clip_index_end]
    for n in 1:size(seq_HXB2,1)
        if(seq_HXB2[n] != "-")
            index_HXB2 +=1
        end
        if(clip_index_start<=index_HXB2 && index_HXB2 <= clip_index_end )
            indices_clipping[n] = index_HXB2
        end
    end
    
    MSA2 = []
    for n in 1:size(MSA,1)
        seq_new = copy(MSA[n][indices_clipping .!= 0] )
        MSA2 = push!(MSA2, seq_new)
    end
    indices_clipping2 = indices_clipping[indices_clipping .!= 0];
    seq_HXB2_clipped = seq_HXB2[indices_clipping .!= 0]
    return (INDEX_HXB2[1], seq_HXB2, seq_HXB2_clipped, MSA2, indices_clipping2)
end;

"""
 Extract the collecting-time from the set of headders
 Remove sequences without the collecting-time information
 Ordering MSA and headder using the extracted collecting-time. 
 Also remove the HXB2 from them. 
"""
function get_collecting_time_ordering_MSA(INDEX_HXB2, MSA, seq_headder)
    collecting_time = []
    index_containing_time = []
    ########## ---- collecting time filter start -- ##########
    # This function may need to slightly change as changing the input alignments.
    for n in 1:size(seq_headder,1)
        if(n != INDEX_HXB2)
            x = split(split(seq_headder[n], ".")[4], "_")
            t = parse(Int, split(x[2], "w")[2] )
            index_containing_time = push!(index_containing_time, n)
            collecting_time = push!(collecting_time, t)
        end
    end
    
    n_remained = size(index_containing_time,1)
    n_tot = size(seq_headder,1)-1
    n_removed = n_tot - n_remained
    # Check that the entire sequences are extracted from the 
    @printf "remove %d sequences out of %d total sequences, then left %d sequences." n_removed n_tot n_remained ;
    sort_time_ordering = sortperm(collecting_time);
    index_containing_time_sort = index_containing_time[sort_time_ordering]
    collecting_time_sort = collecting_time[sort_time_ordering];
    MSA_out = MSA[index_containing_time_sort]
    seq_headder_out = seq_headder[index_containing_time_sort]
    return (collecting_time_sort, MSA_out, seq_headder_out)
end;

""" Get the T/F sequence such that each postions are occupeid by most likely observed letters among the sequences that were collected at the first medical check. 
"""
function get_TF(MSA, index_first_seq, NUC)
    TFseq = []
    for i in 1:size(MSA[1],1)
        tot_var_on_i = []
        for n in index_first_seq
            tot_var_on_i = push!(tot_var_on_i, MSA[n][i])
        end
        a_argmax = argmax([count(tot_var_on_i .== a) for a in NUC ])
        TFseq = push!(TFseq, NUC[a_argmax])
    end
    return TFseq
end;


"""  Get the consensus sequence: each positions are occupeid by the most frequently observed letters.
"""
function get_cons_seq(MSA, seq_headder, NUC)
    #obtain sequence number for those are not HXB2
    index_nonHXB2_MSA = collect(1:size(MSA,1))[isnothing.(match.(r"HXB2", seq_headder))];
    cons_seq = []
    for i in 1:size(MSA[1],1)
        variable_at_i = []
        for n in index_nonHXB2_MSA
            variable_at_i = push!(variable_at_i, MSA[n][i]) 
        end
        a_argmax = argmax([count(variable_at_i .== a) for a in NUC ])
        cons_seq = push!(cons_seq, NUC[a_argmax])
    end
    return cons_seq
end;

"""  Get edge gap indices, where gaps are in the part of streched gap from of the edges (begining or end). 
"""
function get_edge_gap(MSA, seq_headder)
    edge_gap_list = [] 
    edge_def = 200
    index_nonHXB2_MSA = collect(1:size(MSA,1))[isnothing.(match.(r"HXB2", seq_headder))];
    MSA_nonHXB2 = copy(MSA[index_nonHXB2_MSA])

    L = length(MSA_nonHXB2[1]);
    n_seq = length(MSA_nonHXB2)

    for i in 1:L
        edge_gap = false
        gap_seqs = collect(1:n_seq)[[MSA_nonHXB2[n][i] for n in 1:n_seq] .== "-"]
        if( length(gap_seqs) > 0 && (i<=edge_def || (L-i)<=edge_def))
            edge_gap = true
            if(i<=edge_def)
                for s in gap_seqs
                    if(count(MSA_nonHXB2[s][1:i] .== "-")<i)
                        edge_gap = false
                        break
                    end
                end
            end
            if((L-i)<=edge_def)
                for s in gap_seqs
                    if(count(MSA_nonHXB2[s][i:end] .== "-")<L+1-i)
                        edge_gap = false
                        break
                    end
                end
            end
        end
        push!(edge_gap_list, edge_gap)
    end
    return edge_gap_list
end


"""
   Replace the ambigous letters based on the biogenetical knowledge as well as the most frequently observed letters at the corresponding position.
"""
function impute_ambigous_subrouti(n, i, MSA)
    x = MSA[n][i]
    ensemble_ambious_site = [MSA[m][i] for m in 1:size(MSA,1)];
    output_var = "-"
    if(x=="R")
        if(count(ensemble_ambious_site .== "A") > count(ensemble_ambious_site .== "G"))
            output_var = "A"
        else
            output_var = "G"
        end
    end
    if(x=="Y")
        if(count(ensemble_ambious_site .== "T") > count(ensemble_ambious_site .== "C"))
            output_var = "T"
        else
            output_var = "C"
        end
    end
    if(x=="K")
        if(count(ensemble_ambious_site .== "G") > count(ensemble_ambious_site .== "T"))
            output_var = "G"
        else
            output_var = "T"
        end
    end
    if(x=="M")
        if(count(ensemble_ambious_site .== "A") > count(ensemble_ambious_site .== "C"))
            output_var = "A"
        else
            output_var = "C"
        end
    end
    if(x=="S")
        if(count(ensemble_ambious_site .== "G") > count(ensemble_ambious_site .== "C"))
            output_var = "G"
        else
            output_var = "C"
        end
    end
    if(x=="W")
        if(count(ensemble_ambious_site .== "A") > count(ensemble_ambious_site .== "T"))
            output_var = "A"
        else
            output_var = "T"
        end
    end
    @printf "replace %s by %s sequence %d at site %d \n" x output_var n i
    return output_var
end;


function impute_ambigous(MSA, NUC)
    MSA_out = copy(MSA)
    for n in 1:size(MSA, 1)
        for i in 1:size(MSA[1],1)
            if(MSA[n][i] ∉ NUC)        
                MSA_out[n][i] = impute_ambigous_subrouti(n,i,MSA) 
            end
        end
    end
    MSA_out
end;


""" Filtering the sequences: removing the sequences including more than MAX_GAP_NUM. Also, removing positions where the occupancy of the gap is more than MAX_GAP_NUM. 
"""
function get_filterd_MSA(cons_seq, MSA, seq_headder, MAX_GAP_NUM, MAX_GAP_FREQ, HXB2_indices)
    #--- remove sequences those have more than MAX_GAP_NUM gaps. ------#
    index_nonHXB2_MSA = collect(1:size(MSA,1))[isnothing.(match.(r"HXB2", seq_headder))];
    index_nonHXB2_MSA_filterd = []
    n_gap_cons = count(cons_seq .== "-")
    ### ---- remove sequences that involve a lot of gaps ---- ##
    for n in index_nonHXB2_MSA
        n_gap = count(MSA[n] .== "-")
        if(n_gap < MAX_GAP_NUM+n_gap_cons)
            index_nonHXB2_MSA_filterd = push!(index_nonHXB2_MSA_filterd, n)
        end
    end

    n_eff_seq = size(index_nonHXB2_MSA_filterd,1)
    index_sites_filterd = []
    ### ---- remove positions that contains a lot of gaps, also correspond to gap in HXB2 index ---- ##
    for i in 1:size(MSA[1], 1)
        variables_at_i = []
        for n in index_nonHXB2_MSA_filterd
            variables_at_i = push!(variables_at_i, MSA[n][i])
        end
        n_gap = count(variables_at_i .== "-")

        if( n_gap/n_eff_seq<MAX_GAP_FREQ && seq_HXB2[i] != "-")
            index_sites_filterd = push!(index_sites_filterd, i)
        end
    end
    
    HXB2_indices_fliterd = HXB2_indices[index_sites_filterd]
    
    n_row = size(index_nonHXB2_MSA_filterd,1)
    n_col = size(index_sites_filterd,1)
    n_tot = n_row * n_col
    MSA_filterd = reshape(["a" for _ in 1:n_tot], (n_row, n_col));

    for m in 1:size(index_nonHXB2_MSA_filterd, 1)
        n = index_nonHXB2_MSA_filterd[m]
        MSA_filterd[m, :] = MSA[n][index_sites_filterd]
   end
    seq_headder_filterd =  seq_headder[index_nonHXB2_MSA_filterd]
    return (index_nonHXB2_MSA_filterd, index_sites_filterd, MSA_filterd, seq_headder_filterd, HXB2_indices_fliterd)
end;

# This might be able to use the set of protein letters.
function NUC2NUM(x)
    myflag = true
    if(x=="-")
        myflag = false
        return 0
    end
    if(x=="A")
        myflag = false
        return 1
    end
    if(x=="C")
        myflag = false
        return 2
    end
    if(x=="G")
        myflag = false
        return 3
    end
    if(x=="T")
        myflag = false
        return 4
    end
    if(myflag)
        @show x
    end
end;

"""
   Selecting only the polymorophic sites, where multiple letters can be observed at the specific sites. 
"""
function filter_polymorophic(MSA_filterd, TF_filterd, cons_filterd, seq_HXB2_filterd, HXB2_indices_fliterd)
    ### select sites those contain multiple variables, i.e. polymorophic sites
    num_var_each_sites = [length(unique(MSA_filterd[:,i])) for i in 1:size(MSA_filterd,2)  ];
    index_polymorophic = collect(1:size(num_var_each_sites, 1))[num_var_each_sites .> 1];

    #TF_filterd = TF_seq[index_sites_filterd]
    #cons_filterd = cons_seq[index_sites_filterd];
    #seq_HXB2_filterd = seq_HXB2_clipped[index_sites_filterd];
    @show size(MSA_filterd)
    @show size(TF_filterd), size(cons_filterd), size(seq_HXB2_filterd), size(HXB2_indices_fliterd)

    index_polymorophic = collect(1:size(num_var_each_sites, 1))[num_var_each_sites .> 1];
    n_poly = size(index_polymorophic,1)
    n_mono = size(num_var_each_sites,1) - n_poly
    @printf "Remove: %d monomorophic sites. \n Remain: %d polymorophic sites. " n_mono n_poly

    TF_polymorophic = TF_filterd[index_polymorophic]
    cons_polymorophic = cons_filterd[index_polymorophic]
    seq_HXB2_polymorophic = seq_HXB2_filterd[index_polymorophic]
    HXB2_indices_polymorophic = HXB2_indices_fliterd[index_polymorophic]
    MSA_polymorphic = MSA_filterd[:, index_polymorophic];
    return (MSA_polymorphic, HXB2_indices_polymorophic, seq_HXB2_polymorophic, TF_polymorophic, cons_polymorophic)
end;

"""
Write the final result of the MSA (after the cliping, filtering, etc.) into a fasta file that include also HXB2, CONSENSUS, and T/F sequences. 
"""
function MSA_write(fname_out_fasta, MSA_polymorphic, seq_headder, seq_headder_filterd, INDEX_HXB2, seq_HXB2_polymorophic, cons_polymorophic, TF_polymorophic)
    n_seq = size(MSA_polymorphic,1)
    L_seq = size(MSA_polymorphic,2)
    fout = open(fname_out_fasta, "w")

    # output HXB2
    println(fout, seq_headder[INDEX_HXB2])
    for i in 1:L_seq
        print(fout, seq_HXB2_polymorophic[i])
    end
    println(fout, "")

    # output CONSENSUS
    println(fout, ">TF")
    for i in 1:L_seq
        print(fout, TF_polymorophic[i])
    end
    println(fout, "")

    # output CONSENSUS
    println(fout, ">CONSENSUS")
    for i in 1:L_seq
        print(fout, cons_polymorophic[i])
    end
    println(fout, "")

        # individual alignments
    for n in 1:n_seq
        println(fout, seq_headder_filterd[n])
        for i in 1:L_seq
            print(fout, MSA_polymorphic[n, i])
        end
        println(fout, "")
    end

    close(fout);
end;

""" Write the final sequence result into a MPL analysis foramt, including the collecting time, the number of sequences and the numerically represented sequences.
"""

function MSA_MPLseq_write(fname_out_num, MSA_polymorphic, times)
    fout = open(fname_out_num, "w")
    n_seq = size(MSA_polymorphic, 1)
    L_seq = size(MSA_polymorphic, 2)
    times = collecting_time_sort[index_nonHXB2_MSA_filterd]
    for n in 1:n_seq
        print(fout, times[n], "\t", 1, "\t")
        for i in 1:L_seq
            print(fout, NUC2NUM(MSA_polymorphic[n, i]), " ")
        end
        println(fout, "")
    end
    close(fout);
    return
end;

""" This function returns a set of time steps based on the headder of seqeunces
	the key indicater is written as "something+w", this function distinguish "w+something.
"""
function get_time_precsely(seq_headder)
    seq_id_w_time = [] # this array should contain the number of sequence and the estiamted date
    seq_time_w_time = []

    for n in 1:size(seq_headder,1)
        splited_components = [split(seq_headder[n], ('_','.'))]
        length_components = length(splited_components)
        array_tells_words_match = .! isnothing.(match.(r"w", splited_components[1]))

        length_matched = count(array_tells_words_match) 

        if(length_matched>0)
            sth_including_w = splited_components[1][array_tells_words_match]
            is_Nw_or_wN = length.(split(sth_including_w[1], "w"))
            if(is_Nw_or_wN[2] != 0) # the selected one is "5w" not "w5"
                estimated_date = parse(Int, split(sth_including_w[1], "w")[2])            
                push!(seq_id_w_time, n)
                push!(seq_time_w_time, estimated_date)
            end

        end
    end
    return (seq_id_w_time, seq_time_w_time)
end;

"""
    This funciton returns possible mutations with a single mutation 
"""
function get_possible_codon_wiht_single_mut(codon_in) 
    mutants_set = []
    
    substitutions_i = filter!(x -> x != codon_in[1], copy(NUC_dna))
    mutants_set_i = copy(substitutions_i .* (codon_in[2] * codon_in[3]))
    for mut in mutants_set_i
        push!(mutants_set, mut)
    end

    substitutions_i = filter!(x -> x != codon_in[2], copy(NUC_dna))
    mutants_set_i = copy(codon_in[1] .* substitutions_i .* codon_in[3])
    for mut in mutants_set_i
        push!(mutants_set, mut)
    end

    substitutions_i = filter!(x -> x != codon_in[3], copy(NUC_dna))
    mutants_set_i = copy((codon_in[1]*codon_in[2]) .* substitutions_i)
    for mut in mutants_set_i
        push!(mutants_set, mut)
    end
    return mutants_set
end;

""" Get dictionary that gives possible codon varieties given amino acid symbol.
"""
function get_AA_to_codon()
    NUC_dna = ["A", "C", "G", "T"]
    A_set = "GC" .* NUC_dna
    C_set = "TG" .* ["C","T"]
    D_set = "GA" .* ["C","T"]
    E_set = "GA" .* ["A","G"]
    F_set = "TT" .* ["C","T"]
    G_set = "GG" .* NUC_dna
    H_set = "CA" .* ["C", "T"]
    I_set = "AT" .* ["A","C","T"]
    K_set = "AA" .* ["A","G"]
    L_set = ["TT" .* ["A","G"]; "CT" .* copy(NUC_dna)]
    M_set = ["ATG"]
    N_set = "AA" .* ["C", "T"]
    P_set = "CC" .* NUC_dna
    Q_set = "CA" .* ["A", "G"]
    R_set = ["CG" .* NUC_dna; "AG" .* ["A", "G"]]
    S_set = ["AG" .* ["C", "T"]; "TC" .* NUC_dna]
    T_set = "AC" .* NUC_dna
    V_set = "GT" .* NUC_dna
    W_set = ["TGG"]
    Y_set = "TA" .* ["C", "T"];

    AA_to_codon = Dict( "A"=>A_set, "C"=>C_set, "D"=>D_set, "E"=>E_set, "F"=>F_set,
                        "G"=>G_set, "H"=>H_set, "I"=>I_set, "K"=>K_set, "L"=>L_set,
                        "M"=>M_set, "N"=>N_set, "P"=>P_set, "Q"=>Q_set, "R"=>R_set,
                        "S"=>S_set, "T"=>T_set, "V"=>V_set, "W"=>W_set, "Y"=>Y_set);

    return AA_to_codon
end;


AA2NUC = Dict( 
    "A"=>["GCT", "GCC", "GCA", "GCG"], 
    "C"=>["TGT", "TGC"],
    "D"=>["GAT", "GAC"],
    "E"=>["GAA", "GAG"],
	"F"=>["TTC", "TTT"],
    "G"=>["GGT", "GGC", "GGA", "GGG"],
    "H"=>["CAT", "CAC"],
    "I"=>["ATA", "ATC", "ATT"],
    "K"=>["AAA", "AAG"],
    "L"=>["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
	"M"=>["ATG"],
    "N"=>["AAC", "AAT"],
    "P"=>["CCT", "CCA", "CCC", "CCG"],
    "Q"=>["CAA", "CAG"],
    "R"=>["CGA", "CGT", "CGG", "CGC", "AGG", "AGA"],
    "S"=>["TCA", "TCT", "TCG", "TCC", "AGC", "AGT"],
    "T"=>["ACA", "ACT", "ACG", "ACC"],
    "V"=>["GTA", "GTT", "GTG", "GTC"],
    "W"=>["TGG"], 
    "Y"=>["TAT", "TAC"], 
    "*"=>["TGA", "TAA", "TAG"] );

NUC2AA = Dict( 
    "GCT"=>"A", "GCC"=>"A", "GCA"=>"A", "GCG"=>"A", 
    "TGT"=>"C", "TGC"=>"C",
    "GAT"=>"D", "GAC"=>"D",
    "GAA"=>"E", "GAG"=>"E",
	"TTC"=>"F", "TTT"=>"F",
    "GGT"=>"G", "GGC"=>"G", "GGA"=>"G", "GGG"=>"G",
    "CAT"=>"H", "CAC"=>"H",
    "ATA"=>"I", "ATC"=>"I", "ATT"=>"I",
    "AAA"=>"K", "AAG"=>"K",
    "CTA"=>"L", "CTC"=>"L", "CTG"=>"L", "CTT"=>"L", "TTA"=>"L", "TTG"=>"L",
	"ATG"=>"M",
    "AAC"=>"N", "AAT"=>"N",
    "CCT"=>"P", "CCA"=>"P", "CCC"=>"P", "CCG"=>"P",
    "CAA"=>"Q", "CAG"=>"Q",
    "CGA"=>"R", "CGT"=>"R", "CGG"=>"R", "CGC"=>"R", "AGG"=>"R", "AGA"=>"R",
    "TCA"=>"S", "TCT"=>"S", "TCG"=>"S", "TCC"=>"S", "AGC"=>"S", "AGT"=>"S",
    "ACA"=>"T", "ACT"=>"T", "ACG"=>"T", "ACC"=>"T",
    "GTA"=>"V", "GTT"=>"V", "GTG"=>"V", "GTC"=>"V",
    "TGG"=>"W", 
    "TAT"=>"Y", "TAC"=>"Y", 
    "TGA"=>"*", "TAA"=>"*", "TAG"=>"*" );

# --------- Function to return Amino Acids that are accessible by a single nucleotide mutation given an Amino Acid. ------------# 
function inserting_nuc(codon_in_sp, a, i)
    if(i==1)
        return a*codon_in_sp[2]*codon_in_sp[3]
    end
    if(i==2)
        return codon_in_sp[1]*a*codon_in_sp[3]
    end
    if(i==3)
        return codon_in_sp[1]*codon_in_sp[2]*a
    end
end;

function get_possible_AA_with_single_NUC_mutation(input_A, NUC)
    coded_nuc = AA2NUC[input_A];
    single_mutated_AAs = []
    for codon_in in coded_nuc
        codon_in_sp = split(codon_in, "")
        for i in 1:3
            for a in NUC
                if(a != codon_in_sp[i])
                    new_codon = inserting_nuc(codon_in_sp, a, i)
                    AA_new = NUC2AA[new_codon]
                    if(AA_new != input_A)
                        push!(single_mutated_AAs, AA_new)
                        #@show new_codon, AA_new
                    end
                end
            end
        end
    end
    return single_mutated_AAs_unique = unique(single_mutated_AAs)
end
