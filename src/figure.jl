
# ------------------- Gene map used only for CH848 ------------------ #
gene_set_unique = ["nef", "env", "vif", "vpr", "rev", "vpu", "tat", "pol"];
label_gene_set = uppercasefirst.(gene_set_unique)

gen_to_color = Dict(
    "nef" => "orange", 
    "env" => "green", 
    "vpr" => "cyan", 
    "vpu" => "yellow", 
    "vif" => "purple",
    "rev" => "magenta",
    "tat" => "grey",
    "pol" => "brown")

gen_to_frame = Dict(
    "nef" => ["gene_fr1"], 
    "env" => ["gene_fr3"], 
    "vpr" => ["gene_fr1"], 
    "vpu" => ["gene_fr2"], 
    "vif" => ["gene_fr1"],
    "rev" => ["gene_fr2", "gene_fr3"],
    "tat" => ["gene_fr1", "gene_fr2"],
    "pol" => ["gene_fr3"],)

# ------------------------- Figure trajectotry plot ------------------------------- #
function get_trajectory_plot_CH505(csv_raw_CH505, L_fig_tot)
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    colection_time = get_collection_time(csv_raw_CH505);
    n_top = sum(csv_raw_CH505.s_MPL .> 0.01)
    n_time_max = length(colection_time)
    x1_selected = zeros(n_top, length(colection_time))

    for i_t in 1:n_time_max
        x1_selected[:, i_t] = csv_raw_CH505[1:n_top, Symbol("f_at_" * string( colection_time[i_t]))]
    end;

    color_set = []
    alpha_set = []
    for n in 1:n_top
        if(csv_raw_CH505.mutants_AA_fr3[n]=="T415K")
            push!(color_set, "red"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH505.resist_strain_specific_Abs_CH505[n] || csv_raw_CH505.common_mut_SHIV_CH505[n])
            push!(color_set, "purple"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH505.resist_mut_CH235[n])
            push!(color_set, "skyblue"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH505.resist_mut_CH103[n])
            push!(color_set, "green"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH505.N_linked_glycan_minus_fr3[n]>0)
            push!(color_set, "orange"); push!(alpha_set, 0.3)
        else
            push!(color_set, "gray"); push!(alpha_set, 1.0)
        end
    end

    # ---------- Make a subplot ---------------- # 
    my_alpha=0.5
    x1_set_for_figures_CTL = [] 
    x1_set_for_figures_CH103 = [] 
    x1_set_for_figures_CH235 = [] 
    x1_set_for_figures_strain_specific = [] 
    
    for n in sort(collect(1:n_top), rev=true)
    	if(color_set[n] == "red")
		push!(x1_set_for_figures_CTL, x1_selected[n, :])
	end
    	if(color_set[n] == "purple")
		push!(x1_set_for_figures_strain_specific, x1_selected[n, :])
	end
    	if(color_set[n] == "skyblue")
		push!(x1_set_for_figures_CH235, x1_selected[n, :])
	end
    	if(color_set[n] == "green")
		push!(x1_set_for_figures_CH103, x1_selected[n, :])
	end
    end

    p4 = Plots.plot(time_unique, 100*x1_selected[1, :], c=color_set[1], 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
    #    xlim=(0,maximum(time_unique)+30),
        grid=:false,
        margin=my_margin, 
        lw=3,
        alpha=my_alpha,
        tickfontsize=fontsize_reg,
    #    xticks=((0:300:1800)),    
        xticks=:false,
        ylabel=" \n ",
        legend=:false, 
        xlabel=" \n "
    )
    annotate!(-300, 0.3, text("Frequency (%)", :left, rotation=90, fontsize_reg))
    annotate!(900, -30, text("Days", fontsize_reg))
   

    for n in sort(collect(2:n_top), rev=true)
        if(color_set[n] != "gray")
            if(color_set[n] != "orange")
                Plots.plot!(time_unique, 100*x1_selected[n, :], c=color_set[n], 
                    alpha=my_alpha, legend=:false, lw=3,)
            #else
            #    Plots.plot!(time_unique, 100*x1_selected[n, :], c=color_set[n], alpha=alpha_set[n], legend=:false, lw=3,)
            end
        end
    end
    for x in collect(0:300:1800)
        annotate!(x, +1, text("|", 6))
        if(x%600==0)
            annotate!(x, -14, text(string(x), fontsize_reg))
        end
    end
    annotate!(-350, 120, text(L"\textbf{C}", :left, fontsize_label_reg))
    return (p4, time_unique, x1_set_for_figures_CTL, x1_set_for_figures_CH103, x1_set_for_figures_CH235, x1_set_for_figures_strain_specific)
end

function get_selection_vs_time_plot_CH505(csv_raw_CH505, L_fig_tot)
    # --------------- Plot comparing selection and collection time ----------------- #  
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    idx_T415K_mached = [!isnothing(match(r"T415K", x)) for x in csv_raw_CH505.mutants_AA_fr3 ]

    idx_N_add = (csv_raw_CH505.N_linked_glycan_plus_fr1 .> 0) .|| (csv_raw_CH505.N_linked_glycan_plus_fr2 .> 0) .|| (csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0)
    idx_N_rem = (csv_raw_CH505.N_linked_glycan_minus_fr1 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr2 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0)

    c_jitter = 10
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    my_ms = 5
    myalpha= 0.5
    idx = (idx_N_add .|| idx_N_rem)
    xaxis = csv_raw_CH505.detected_date[idx]
    Plots.scatter(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="orange",
        m=:circle,
        ms = my_ms,
        alpha=0.3,
        label="Glycan",
        subplot=1)

    idx = csv_raw_CH505.resist_mut_CH103
    xaxis = csv_raw_CH505.detected_date[idx]
    my_ms =5
    p5 = Plots.scatter(xaxis .+ c_jitter * randn(length(xaxis)), 100*csv_raw_CH505.s_MPL[idx], 
        c="green", 
        m=:circle,
        ms = my_ms,
        alpha=myalpha,
        label="CH103",
        markerstrokewidth=0,
        subplot=1
    )
    time_CH103_figure = copy(xaxis)
    selection_CH103_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = csv_raw_CH505.resist_mut_CH235
    xaxis = csv_raw_CH505.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="skyblue", 
        label="CH235",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_CH235_figure = copy(xaxis)
    selection_CH235_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = (csv_raw_CH505.resist_strain_specific_Abs_CH505 .|| csv_raw_CH505.common_mut_SHIV_CH505)
    xaxis = csv_raw_CH505.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="purple",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        label="Strain-specific Abs",
        markerstrokewidth=0,
        subplot=1)
    time_autologous_figure = copy(xaxis)
    selection_autologous_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = idx_T415K_mached
    xaxis = csv_raw_CH505.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="red",
        m=:circle,
        alpha=myalpha,
        ms = my_ms+2,
        label="CTL",
        markerstrokewidth=0,
        subplot=1, 
        #xrotation=20,
        xlim=(0, maximum(time_unique)+30),
    #    xticks=((0:900:1800)),
        xticks=:false,
        xlabel = " ",
    #    ylabel = "Selection (%)",
        ylim=(-3, 10), 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        #legend=:false, 
        legend=(0.7, 0.9), 
        margin=my_margin, 
        legendfontsize=fontsize_reg-1,
        tickfontsize=fontsize_reg,
    );
    time_CTL_figure = copy(xaxis)
    selection_CTL_figure = copy(csv_raw_CH505.s_MPL[idx])


    for x in collect(0:300:1800)
        annotate!(x, -3, text("|", 6))
    end
    annotate!(-180, 4, text("Selection (%)", fontsize_reg, rotation=90))
    annotate!(870, -4.2, text("Time mutation was first observed \n(days after infection)", fontsize_reg))
    annotate!(-300, 10, text(L"\textbf{A}", :left, fontsize_label_reg))
    return (p5, time_CTL_figure, time_CH103_figure, time_CH235_figure, time_autologous_figure, 
		selection_CTL_figure, selection_CH103_figure, selection_CH235_figure, selection_autologous_figure 
    ) 
end;

function get_violin_plot_CH505(csv_raw_CH505, L_fig_tot)
    # --------------- Box plot comparing selection distribuiton ----------------- # 
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.25*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    myalpha= 0.3
    my_ms = 6

    s_res_CH103 = csv_raw_CH505.s_MPL[csv_raw_CH505.resist_mut_CH103 .* (csv_raw_CH505.nonsynonymous .> 0)]
    s_res_CH235 = csv_raw_CH505.s_MPL[csv_raw_CH505.resist_mut_CH235 .* (csv_raw_CH505.nonsynonymous .> 0)]
    s_res_spcfc = csv_raw_CH505.s_MPL[(csv_raw_CH505.resist_strain_specific_Abs_CH505 .|| csv_raw_CH505.common_mut_SHIV_CH505) .* (csv_raw_CH505.nonsynonymous .> 0)];

    v1 = 100*s_res_spcfc; n1=length(v1)
    v2 = 100*s_res_CH235; n2=length(v2)
    v3 = 100*s_res_CH103; n3=length(v3)

    # Combine them into a single vector y
    y = vcat(v1, v2, v3)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)), fill(3, length(v3)))

    barwidth = 0.5
    width = 0.4 * barwidth
    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    # jittered x
    p6=scatter(x_jitter[1:n1], y[1:n1], c=:purple, markerstrokewidth=0, label="Strain-Specific Abs", 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        margin=my_margin, 
        ms = my_ms,
        size=(L_fig, Int(ceil(0.6*600))),
        alpha=myalpha,
        ylim=(-2,10),
        m=:circle,
    #     xticks=(1:3, ["Strain-\n specific", "CH235", "CH103"]), 
        xticks=:false,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg-2,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:skyblue, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="CH235")
    scatter!(x_jitter[(n1+n2+1):end], y[(n1+n2+1):end], c=:green, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="CH103")

    # jittering lines up with a violin plot
    boxplot!(x[1:n1], y[1:n1], alpha=0.3, c=:purple, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.3, c="skyblue", label=:false)
    boxplot!(x[(n1+n2+1):end], y[(n1+n2+1):end], alpha=0.3, c=:green, label=:false)

    annotate!(0, 10, text(L"\textbf{B}", :left, fontsize_label_reg))
    annotate!(-0.5, -4.6, text("Strain\n specific", :left, fontsize_reg, rotation=45))
    annotate!(1, -4.6, text("CH235", :left, fontsize_reg, rotation=45))
    annotate!(2.5, -4.6, text("CH103", :left, fontsize_reg, rotation=45));

    return (p6, s_res_CH103, s_res_CH235, s_res_spcfc) 
end;
# -------------------------- Figure: trajectotry plots with glycan --------------------------#
# ------------------------- Figure trajectotry plot ------------------------------- #
function get_trajectory_plot_CH505_with_glycan(csv_raw_CH505, L_fig_tot)
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    colection_time = get_collection_time(csv_raw_CH505);
    #n_top = sum(csv_raw_CH505.s_MPL .> 0.01)
    n_top = length(csv_raw_CH505.s_MPL)
    n_time_max = length(colection_time)
    x1_selected = zeros(n_top, length(colection_time))
    for i_t in 1:n_time_max
        x1_selected[:, i_t] = csv_raw_CH505[1:n_top, Symbol("f_at_" * string( colection_time[i_t]))]
    end;
    idx_glycan = (csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0 .|| csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0) .&& (csv_raw_CH505.nonsynonymous .> 0) 
    idx_CTL = [ ]; idx_CH103 = [ ]; idx_CH235 = [ ]; idx_strain_specific = [ ]
    for n in 1:n_top
        x = csv_raw_CH505.mutants_AA_fr3[n]
        if( (extract_integer(x) ∈ 409:418) && (csv_raw_CH505.nonsynonymous[n]>0))
            push!(idx_CTL, true)
        else
            push!(idx_CTL, false)
        end 
        if(csv_raw_CH505.resist_mut_CH103[n] && (csv_raw_CH505.detected_date[n] >90 ) && (csv_raw_CH505.nonsynonymous[n]>0))
            push!(idx_CH103, true)
        else
            push!(idx_CH103, false)
        end
        if(csv_raw_CH505.resist_mut_CH235[n] && (csv_raw_CH505.nonsynonymous[n]>0))
            push!(idx_CH235, true)
        else
            push!(idx_CH235, false)
        end
        if((csv_raw_CH505.resist_strain_specific_Abs_CH505[n] || csv_raw_CH505.common_mut_SHIV_CH505[n] ) && (csv_raw_CH505.nonsynonymous[n]>0))
            push!(idx_strain_specific, true)
        else
            push!(idx_strain_specific, false)
        end
    end
    
    # ---------- Make a subplot ---------------- # 
    my_alpha=0.5
    x1_set_for_figures_CTL = [] 
    x1_set_for_figures_CH103 = [] 
    x1_set_for_figures_CH235 = [] 
    x1_set_for_figures_strain_specific = [] 
    x1_set_for_figures_glycan = [] 
    mutation_labels = [] 
    for n in sort(collect(1:n_top), rev=true)
    	#if(color_set[n] == "red")
    	if(idx_CTL[n])  
		    push!(x1_set_for_figures_CTL, x1_selected[n, :])
	    end
    	#if(color_set[n] == "purple")
    	if(idx_strain_specific[n])
		    push!(x1_set_for_figures_strain_specific, x1_selected[n, :])
	    end
    	#if(color_set[n] == "skyblue")
    	if(idx_CH235[n])
		    push!(x1_set_for_figures_CH235, x1_selected[n, :])
	    end
    	#if(color_set[n] == "green")
    	if(idx_CH103[n])
		    push!(x1_set_for_figures_CH103, x1_selected[n, :])
	    end
    	#if(color_set[n] == "orange")
    	if(idx_glycan[n])
		    push!(x1_set_for_figures_glycan, x1_selected[n, :])
	    end
    end
    p4 = Plots.plot(
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
    #    xlim=(0,maximum(time_unique)+30),
        grid=:false,
        margin=my_margin, 
        lw=3,
        alpha=my_alpha,
        tickfontsize=fontsize_reg,
    #    xticks=((0:300:1800)),    
        xticks=:false,
        ylabel=" \n ",
        legend=:false, 
        xlabel=" \n "
    )
    annotate!(-300, 0.3, text("Frequency (%)", :left, rotation=90, fontsize_reg))
    annotate!(900, -30, text("Days", fontsize_reg))
   
    for n in 1:n_top 
        if(idx_CH103[n])
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:green, alpha=my_alpha, legend=:false, lw=3,)
        end
        if(idx_CH235[n])
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:skyblue, alpha=my_alpha, legend=:false, lw=3,)
        end
        if(idx_strain_specific[n])
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:purple, alpha=my_alpha, legend=:false, lw=3,)
        end
        if(idx_glycan[n])
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:orange, alpha=my_alpha, legend=:false, lw=3,)
        end
        if(idx_CTL[n])
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:red, alpha=my_alpha, legend=:false, lw=3,)
        end
    end
    for x in collect(0:300:1800)
        annotate!(x, +1, text("|", 6))
        if(x%600==0)
            annotate!(x, -14, text(string(x), fontsize_reg))
        end
    end
    annotate!(-350, 120, text(L"\textbf{C}", :left, fontsize_label_reg))
    return (p4, time_unique, x1_set_for_figures_CTL, x1_set_for_figures_CH103, x1_set_for_figures_CH235, x1_set_for_figures_strain_specific, x1_set_for_figures_glycan)
end

function get_selection_vs_time_plot_CH505_with_glycan(csv_raw_CH505, L_fig_tot)
    # --------------- Plot comparing selection and collection time ----------------- #  
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    #idx_T415K_mached = [!isnothing(match(r"T415K", x)) for x in csv_raw_CH505.mutants_AA_fr3 ]
    #idx_CTL = [extract_integer(x) ∈ 409:418 for x in csv_raw_CH505.mutants_AA_fr3 ]
    idx_CTL = [ ]
    name_mutation_CTL = []
    for n in 1:length(csv_raw_CH505.mutants_AA_fr3)
        x = csv_raw_CH505.mutants_AA_fr3[n]
        if((extract_integer(x) ∈ 409:418) * (csv_raw_CH505.nonsynonymous[n]>0) )
            push!(idx_CTL, true)
            push!(name_mutation_CTL, x)
        else
            push!(idx_CTL, false)
        end 
    end
    
    #idx_N_add = (csv_raw_CH505.N_linked_glycan_plus_fr1 .> 0) .|| (csv_raw_CH505.N_linked_glycan_plus_fr2 .> 0) .|| (csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0)
    #idx_N_rem = (csv_raw_CH505.N_linked_glycan_minus_fr1 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr2 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0)
    #idx_glycan = (csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0) 
    c_jitter = 10
    L_fig = Int(ceil(L_fig_tot*0.75*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    my_ms = 5
    myalpha= 0.5
    idx = ( (csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0) .|| (csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0) ) .&& (csv_raw_CH505.nonsynonymous .> 0  )
    xaxis = csv_raw_CH505.detected_date[idx]
    name_mutation_glycan = csv_raw_CH505.mutants_AA_fr3[idx]
    Plots.scatter(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="orange",
        m=:circle,
        ms = my_ms,
        alpha=0.3,
        markerstrokewidth=0,
        label="Glycan",
        subplot=1)
    time_glycan_figure = copy(xaxis)
    selection_glycan_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = csv_raw_CH505.resist_mut_CH103 .* (csv_raw_CH505.nonsynonymous .> 0  )
    xaxis = csv_raw_CH505.detected_date[idx]
    idx1 = xaxis .> 90 
    name_mutation_CH103 = csv_raw_CH505.mutants_AA_fr3[idx][idx1]
    my_ms =5
    p5 = Plots.scatter!(xaxis[idx1] .+ c_jitter * randn(length(xaxis[idx1])), 100*csv_raw_CH505.s_MPL[idx][idx1], 
        c="green", 
        m=:circle,
        ms = my_ms,
        alpha=myalpha,
        label="CH103",
        markerstrokewidth=0,
        subplot=1
    )
    time_CH103_figure = copy(xaxis[idx1])
    selection_CH103_figure = copy(csv_raw_CH505.s_MPL[idx][idx1])

    idx = csv_raw_CH505.resist_mut_CH235 .* (csv_raw_CH505.nonsynonymous .> 0  )
    xaxis = csv_raw_CH505.detected_date[idx]
    name_mutation_CH235 = csv_raw_CH505.mutants_AA_fr3[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="skyblue", 
        label="CH235",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_CH235_figure = copy(xaxis)
    selection_CH235_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = (csv_raw_CH505.resist_strain_specific_Abs_CH505 .|| csv_raw_CH505.common_mut_SHIV_CH505) .* (csv_raw_CH505.nonsynonymous .> 0  )
    xaxis = csv_raw_CH505.detected_date[idx]
    name_mutation_strain_specific = csv_raw_CH505.mutants_AA_fr3[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="purple",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        label="Strain-specific Abs",
        markerstrokewidth=0,
        subplot=1)
    time_autologous_figure = copy(xaxis)
    selection_autologous_figure = copy(csv_raw_CH505.s_MPL[idx])

    idx = idx_CTL .* (csv_raw_CH505.nonsynonymous .> 0  )
    xaxis = csv_raw_CH505.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH505.s_MPL[idx],
        c="red",
        m=:circle,
        alpha=myalpha,
        ms = my_ms+2,
        label="CTL",
        markerstrokewidth=0,
        subplot=1, 
        #xrotation=20,
        xlim=(0, maximum(time_unique)+30),
    #    xticks=((0:900:1800)),
        xticks=:false,
        xlabel = " ",
    #    ylabel = "Selection (%)",
        ylim=(-3, 10), 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        #legend=:false, 
        legend=(0.7, 0.9), 
        margin=my_margin, 
        legendfontsize=fontsize_reg-1,
        tickfontsize=fontsize_reg,
    );
    time_CTL_figure = copy(xaxis)
    selection_CTL_figure = copy(csv_raw_CH505.s_MPL[idx])


    for x in collect(0:300:1800)
        annotate!(x, -3, text("|", 6))
    end
    annotate!(-180, 4, text("Selection (%)", fontsize_reg, rotation=90))
    annotate!(870, -4.2, text("Time mutation was first observed \n(days after infection)", fontsize_reg))
    annotate!(-300, 10, text(L"\textbf{A}", :left, fontsize_label_reg))
    return (p5, time_CTL_figure, time_CH103_figure, time_CH235_figure, time_autologous_figure, time_glycan_figure,  
		selection_CTL_figure, selection_CH103_figure, selection_CH235_figure, selection_autologous_figure, selection_glycan_figure, 
        name_mutation_CTL, name_mutation_CH103, name_mutation_CH235, name_mutation_strain_specific, name_mutation_glycan
    ) 
end;

function get_violin_plot_CH505_with_glycan(csv_raw_CH505, L_fig_tot)
    # --------------- Box plot comparing selection distribuiton ----------------- # 
    #L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.25*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    myalpha= 0.3
    my_ms = 6
    #idx_CTL = [extract_integer(x) ∈ 409:418 for x in csv_raw_CH505.mutants_AA_fr3 ]
    idx_CTL = [ ]
    for n in 1:length(csv_raw_CH505.mutants_AA_fr3)
        x = csv_raw_CH505.mutants_AA_fr3[n]
        if((extract_integer(x) ∈ 409:418) && (csv_raw_CH505.nonsynonymous[n] > 0) )
            push!(idx_CTL, true)
        else
            push!(idx_CTL, false)
        end 
    end
    idx_CTL = Bool.(idx_CTL)
    s_res_CH103 = csv_raw_CH505.s_MPL[csv_raw_CH505.resist_mut_CH103 .&& (csv_raw_CH505.nonsynonymous .> 0) .&& (csv_raw_CH505.detected_date .> 90)]
    s_res_CH235 = csv_raw_CH505.s_MPL[csv_raw_CH505.resist_mut_CH235 .&& (csv_raw_CH505.nonsynonymous .> 0)]
    s_res_spcfc = csv_raw_CH505.s_MPL[(csv_raw_CH505.resist_strain_specific_Abs_CH505 .|| csv_raw_CH505.common_mut_SHIV_CH505) .&& (csv_raw_CH505.nonsynonymous .> 0)];
    s_res_glycan = csv_raw_CH505.s_MPL[(csv_raw_CH505.N_linked_glycan_minus_fr3 .> 0 .|| csv_raw_CH505.N_linked_glycan_plus_fr3 .> 0) .&& (csv_raw_CH505.nonsynonymous .> 0) ];
    s_res_CTL = csv_raw_CH505.s_MPL[idx_CTL];

    v1 = 100*s_res_spcfc; n1=length(v1)
    v2 = 100*s_res_CH235; n2=length(v2)
    v3 = 100*s_res_CH103; n3=length(v3)
    v4 = 100*s_res_glycan; n4=length(v4)
    v5 = 100*s_res_CTL; n5=length(v5)

    # Combine them into a single vector y
    y = vcat(v1, v2, v3, v4, v5)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)), fill(3, length(v3)), fill(4, length(v4)), fill(5, length(v5)))

    barwidth = 0.5
    width = 0.4 * barwidth
    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    # jittered x
    p6=scatter(x_jitter[1:n1], y[1:n1], c=:purple, markerstrokewidth=0, label="Strain-Specific Abs", 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        margin=my_margin, 
        ms = my_ms,
        size=(L_fig, Int(ceil(0.6*600))),
        alpha=myalpha,
        ylim=(-2,10),
        m=:circle,
    #     xticks=(1:3, ["Strain-\n specific", "CH235", "CH103"]), 
        xticks=:false,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg-2,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:skyblue, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="CH235")
    scatter!(x_jitter[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], c=:green, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="CH103")
    scatter!(x_jitter[(n1+n2+n3+1):(n1+n2+n3+n4)], y[(n1+n2+n3+1):(n1+n2+n3+n4)], c=:orange, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="Glycan")        
    scatter!(x_jitter[(n1+n2+n3+n4+1):end], y[(n1+n2+n3+n4+1):end], c=:red, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="CTL")        

    # jittering lines up with a violin plot
    boxplot!(x[1:n1], y[1:n1], alpha=0.3, c=:purple, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.3, c="skyblue", label=:false)
    boxplot!(x[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], alpha=0.3, c=:green, label=:false)
    boxplot!(x[(n1+n2+n3+1):(n1+n2+n3+n4)], y[(n1+n2+n3+1):(n1+n2+n3+n4)], alpha=0.3, c=:orange, label=:false)
    boxplot!(x[(n1+n2+n3+n4+1):end], y[(n1+n2+n3+n4+1):end], alpha=0.3, c=:red, label=:false)

    annotate!(0, 10, text(L"\textbf{B}", :left, fontsize_label_reg))
    annotate!(-0.5, -4.6, text("Specific", :left, fontsize_reg, rotation=45))
    annotate!(1, -4.6, text("CH235", :left, fontsize_reg, rotation=45))
    annotate!(2, -4.6, text("CH103", :left, fontsize_reg, rotation=45));
    annotate!(3, -4.6, text("Glycan", :left, fontsize_reg, rotation=45));
    annotate!(4, -4.6, text("CTL", :left, fontsize_reg, rotation=45));

    return (p6, s_res_CH103, s_res_CH235, s_res_spcfc, s_res_glycan, s_res_CTL) 
end;
# -----------------------------------------------------------------------------------------------------------#
# --- This is a case for the x's fitness values --- #
"""
    get_compare_fitness(fname_key_human_and_RMs, fname_dev, dir_name, 
        csv_raw_x, csv_raw_y; flag_include_HIV = false, flag_shuffle = false)
"""
function get_compare_fitness(fname_key_human_and_RMs, fname_dev, dir_name, 
        csv_raw_x, csv_raw_y; flag_include_HIV = false, flag_shuffle = false)
    # May be this cell should be written in a function because the similar process need to apply for other conditions. 
    n_H_RMs = length(fname_key_human_and_RMs) 
        
    date_unique_set = []
    errorbar_tot_x_set = []
    errorbar_tot_y_set = []
    fitness_ind_tot_x_random_set = []
    fitness_ind_tot_y_random_set = []
    mean_fitness_ind_tot_x_random_set = []
    mean_fitness_ind_tot_y_random_set = []
    bool_BNB_nonBNB = []
    for k in 1:n_H_RMs # >2
        key_RM = fname_key_human_and_RMs[k]
            if(occursin(r"^RM\d+$", key_RM) || flag_include_HIV) 
            csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
            date_unique = sort(unique(csv_raw.date))
            push!(date_unique_set, copy(date_unique))

            idx = [x ∈ csv_raw.mutation for x in csv_raw.mutation];
            n_time = length(date_unique)
            
            s_in_x = [(x,y) ∈ zip(csv_raw_x.HXB2, csv_raw_x.PRO) ? 
                    csv_raw_x.s_MPL[ (csv_raw_x.HXB2 .== x) .* (csv_raw_x.PRO .== y) ][1] : 0.0    
                    for (x,y) in zip(csv_raw.HXB2, csv_raw.PRO)];

            s_in_y = [(x,y) ∈ zip(csv_raw_y.HXB2, csv_raw_y.PRO) ? 
                    csv_raw_y.s_MPL[ (csv_raw_y.HXB2 .== x) .* (csv_raw_y.PRO .== y) ][1] : 0.0    
                    for (x,y) in zip(csv_raw.HXB2, csv_raw.PRO)];        
        
            # ----------- Start computing variance of sequencies -------------#
            seq_raw = readdlm(dir_name * key_RM * "-poly.num");
            time_temp = copy(Int.(seq_raw[:,1]))
            seq_temp = Int.(copy(seq_raw[:,3:end])) .+ 1;
            if(flag_shuffle)
                # ---- sampling while keeping the freuqncy ----# 
                seq_temp_shuffled = zeros(size(seq_temp))    
                n_seq_temp = size(seq_temp,1)
                for i in 1:size(seq_temp,2)
                    seq_temp_shuffled[:, i] = shuffle(seq_temp[:, i])
                end
                seq_temp = Int.(seq_temp_shuffled)
            end 

            # --- Set selection coefficients --- #
            L_temp = maximum(csv_raw.polymorphic)
            s_extended_y = zeros(L_temp, q_AA)
            s_extended_x = zeros(L_temp, q_AA)

            for i in 1:length(csv_raw_y.s_MPL)
                i_hxb2 = csv_raw_y.HXB2[i]
                pro_in = csv_raw_y.PRO[i]
                s_temp = csv_raw_y.s_MPL[i]    
                idx_look = (string.(csv_raw.HXB2) .== string(i_hxb2)) .* (string.(csv_raw.PRO) .== string.(pro_in))
                if(count(idx_look)>0)
                    if(count(idx_look)==1)
                        i_poly = csv_raw.polymorphic[idx_look][1]
                        a = AA2NUM[pro_in]
                        s_extended_y[i_poly, a] = s_temp
                    end
                end
            end

            for i in 1:length(csv_raw_x.s_MPL)
                i_hxb2 = csv_raw_x.HXB2[i]
                pro_in = csv_raw_x.PRO[i]
                s_temp = csv_raw_x.s_MPL[i]    
                idx_look = (string.(csv_raw.HXB2) .== string(i_hxb2)) .* (string.(csv_raw.PRO) .== string.(pro_in))
                if(count(idx_look)>0)
                    if(count(idx_look)==1)
                        i_poly = csv_raw.polymorphic[idx_look][1]
                        a = AA2NUM[pro_in]
                        s_extended_x[i_poly, a] = s_temp
                    end
                end
            end

            # --- Obtain fitness values --- #
            f_var_set_x=[]; f_var_set_y=[];
            f_mean_set_x=[]; f_mean_set_y=[];
            for i_t in 1:length(date_unique)
                idx_t = date_unique[i_t] .== time_temp
                n_t = count(idx_t)
                ensemble_t = copy(seq_temp[idx_t, :])
                f_set_x_temp = []; f_set_y_temp = []
                for i_temp in 1:n_t
                    x_t = copy(ensemble_t[i_temp, :])
                    f_y_temp = sum([s_extended_y[i, x_t[i]] for i in 1:L_temp])
                    f_x_temp = sum([s_extended_x[i, x_t[i]] for i in 1:L_temp])
                    push!(f_set_y_temp, f_y_temp)
                    push!(f_set_x_temp, f_x_temp)
                    if(key_RM ∈ fname_dev)
                        push!(bool_BNB_nonBNB, true) else
                        push!(bool_BNB_nonBNB, false)
                    end
                end
                f_var_x = std(f_set_x_temp); f_mean_x = mean(f_set_x_temp)
                f_var_y = std(f_set_y_temp); f_mean_y = mean(f_set_y_temp)
                push!(f_var_set_x, f_var_x); push!(f_mean_set_x, f_mean_x)
                push!(f_var_set_y, f_var_y); push!(f_mean_set_y, f_mean_y)            
                push!(fitness_ind_tot_y_random_set, 100*copy(float.(f_set_y_temp)))
                push!(fitness_ind_tot_x_random_set, 100*copy(float.(f_set_x_temp)))
            end;    
            push!(errorbar_tot_x_set, 100*copy(f_var_set_x))
            push!(errorbar_tot_y_set, 100*copy(f_var_set_y))
            push!(mean_fitness_ind_tot_x_random_set, 100*copy(float.(f_mean_set_x)))
            push!(mean_fitness_ind_tot_y_random_set, 100*copy(float.(f_mean_set_y)))
        end
    end
    
    # --- Make it a single long vector ---- #
    #fitness_ind_tot_x_random_set = copy(reduce(vcat, fitness_ind_tot_x_random_set))
    #fitness_ind_tot_y_random_set = copy(reduce(vcat, fitness_ind_tot_y_random_set));

    return (fitness_ind_tot_x_random_set, fitness_ind_tot_y_random_set, 
            mean_fitness_ind_tot_x_random_set, mean_fitness_ind_tot_y_random_set, 
            errorbar_tot_x_set, errorbar_tot_y_set, 
            date_unique_set, bool_BNB_nonBNB)
end;

function get_heatmap_selected_mutations_CH505(dir_name, L_fig_tot)
    set_VL_increase = ["N130D", "N279D", "K302N", "Y330H", "N334S", "H417R"];
    fname_potent = ["RM5695", "RM6070", "RM6072"]
    fname_not_potent = ["RM6697", "RM6699", "RM6701", "RM6703"]
    mutation_potent_set = []
    mutation_not_potent_set = []
    s_potent_set = []
    s_not_potent_set = []
    for k in 1:length(fname_potent)
        key_RM = fname_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_potent_set, copy(csv_raw.mutation))
    end;
    mutation_potent_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2], mutation_potent_set[3]);
    mutation_bnAbs_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_potent_union = union(mutation_potent_set[1], mutation_potent_set[2], mutation_potent_set[3]);

    for k in 1:length(fname_not_potent)
        key_RM = fname_not_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_not_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_not_potent_set, copy(csv_raw.mutation))
    end;
    mutation_not_potent_intersect = intersect(mutation_not_potent_set[1], mutation_not_potent_set[2], 
        mutation_not_potent_set[3], mutation_not_potent_set[4]);
    mutation_not_potent_union = union(mutation_not_potent_set[1], mutation_not_potent_set[2], 
        mutation_not_potent_set[3], mutation_not_potent_set[4]);
    # --- Get common mutations that observed in potent and not-potent RMs ---#\\\
    mutation_both_potent_and_not_potent = intersect(mutation_potent_union, mutation_not_potent_union);
    # --- Get mutations that are common in potent group but never seen in nont-potent group --- #
    mutation_only_in_potent = setdiff(mutation_potent_intersect, mutation_not_potent_union);

    n_max_both_potent_and_not_potent = 15
    n_max_only_in_potent = 5
    heatmap_both_potent_and_not_potent = zeros(length(fname_key_human_and_RMs), n_max_both_potent_and_not_potent)
    heatmap_only_in_potent = zeros(length(fname_key_human_and_RMs), n_max_only_in_potent);
    heatmap_both_potent_and_not_potent = -0.09 * ones(length(fname_key_human_and_RMs), n_max_both_potent_and_not_potent)
    heatmap_only_in_potent = -0.09 * ones(length(fname_key_human_and_RMs), n_max_only_in_potent);

    for k in 1:n_RMs_max
        if(k>1)
            key_RM = fname_key_human_and_RMs[k]
            csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        end
        if(k==1)
            key_RM = fname_key_human_and_RMs[k]
            csv_raw = DataFrame(CSV.File(dir_name * "RMs-poly.csv"));
        end
        for n in 1:n_max_both_potent_and_not_potent
            mut_temp = mutation_both_potent_and_not_potent[n]
            if(mut_temp ∈ csv_raw.mutation)
                idx = mut_temp .== csv_raw.mutation
                heatmap_both_potent_and_not_potent[n_RMs_max+1-k, n] = csv_raw.s_MPL[idx][1]
            end
        end
        for n in 1:n_max_only_in_potent
            mut_temp = mutation_only_in_potent[n]
            if(mut_temp ∈ csv_raw.mutation)
                idx = mut_temp .== csv_raw.mutation
                heatmap_only_in_potent[n_RMs_max+1-k, n] = csv_raw.s_MPL[idx][1]
            end
        end
    end
    # ====================== Making Figures ========================= #        


    L_fig = Int(0.6*L_fig_tot)
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    i_aa  = extract_integer.(mutation_only_in_potent[1:n_max_only_in_potent])
    idx_sort_aa = sortperm(i_aa)
    l = @layout[b{0.72w, 0.9h} a{0.24w}]
    p1 = heatmap(heatmap_only_in_potent[:, idx_sort_aa], clim=(-0.09,0.09), c=cgrad([:gray80, :white, :red]), 
        legend=:false, xmirror = true, 
        xticks=(1:n_max_only_in_potent, mutation_only_in_potent[1:n_max_only_in_potent][idx_sort_aa]) , 
        xtickfontcolor=:white,
        title=" \n \n \n",
        yaxis=:false, 
    )
    annotate!(n_max_only_in_potent/2+0.5, 0.0, text("- Only potent -", fontsize_reg))
    for i in 1:n_max_only_in_potent
        annotate!(i-0.07, 9-0.37, text(mutation_only_in_potent[idx_sort_aa[i]], :black, fontsize_reg-4, :left, rotation=90))
    end
    heatmap_only_BNB_figure = copy(heatmap_only_in_potent[:, idx_sort_aa])
    mutation_names_only_BNB_figure = mutation_only_in_potent[1:n_max_only_in_potent][idx_sort_aa]

    # ----------------- Annotate Region ------------------ #
    y_annotate = 10.9
    variable_set = ["V1", "V2", "LD", "V3", "V4", "V5"];
    idx_variable_set = [[] for _ in 1:length(variable_set)]
    for j in 1:length(variable_set)
        x = variable_set[j]
        for i in 1:length(i_aa[idx_sort_aa])
            this_var = get_variable_regions(i_aa[idx_sort_aa[i]] )
            if(this_var == x)
                push!(idx_variable_set[j], i)
            end
        end
    end
    for j in 1:length(variable_set)
        if(length(idx_variable_set[j])>0)
            v_left = minimum(idx_variable_set[j])-0.4
            v_right = maximum(idx_variable_set[j])+0.4
            x_annotate_vec = collect(v_left:0.1:v_right)
            Plots.plot!(x_annotate_vec, (y_annotate-0.4) * ones(length(x_annotate_vec)), c=:black, legend=:false, axis=:false, ticks=:false, lw=2)
            annotate!((v_right+v_left)/2-0.5, y_annotate-0.1, text(variable_set[j], fontsize_reg, :left))
        end
    end
    # ----------------------------------------------------- #


    tick_color = [x ∈ set_VL_increase ? "red" : "black"  for x in mutation_both_potent_and_not_potent[1:n_max_both_potent_and_not_potent]]
    i_aa  = extract_integer.(mutation_both_potent_and_not_potent[1:n_max_both_potent_and_not_potent])
    idx_sort_aa = sortperm(i_aa)
    idx_sort_aa = copy(idx_sort_aa[i_aa[idx_sort_aa] .!= 375])
    p2 = heatmap(heatmap_both_potent_and_not_potent[:, idx_sort_aa], clim=(-0.09,0.09), c=cgrad([:gray80, :white, :red]), 
        legend=:false, xmirror = true, 
        yticks=(1:n_H_RMs, reverse(["RMs"; fname_key_human_and_RMs[2:end]])) , 
        xticks=false,
        ylabel=" \n ", 
        ytickfontcolor=:white
    )
    for i in 1:length(idx_sort_aa)
        annotate!(i-0.07, 9-0.37, text(mutation_both_potent_and_not_potent[idx_sort_aa[i]], Symbol(tick_color[idx_sort_aa[i]]), fontsize_reg-4, :left, rotation=90))
    end

    for k in 1:n_H_RMs
        key_RM = fname_key_human_and_RMs[k]
        this_c = "black"
        if(key_RM ∈ fname_dev) this_c = "magenta" end
    #    if(key_RM == "RM6072") this_c = "pink" end
        if(k==1) key_RM = "Joint" end
        if(k!=1) key_RM = "RM"*string(k-1) end
        annotate!(-0.15, n_H_RMs+1-k, text(key_RM, Symbol(this_c), fontsize_reg))
    end
    # ----------------- Annotate Region ------------------ #
    y_annotate = 10.9
    idx_variable_set = [[] for _ in 1:length(variable_set)]
    for j in 1:length(variable_set)
        x = variable_set[j]
        for i in 1:length(i_aa[idx_sort_aa])
            this_var = get_variable_regions(i_aa[idx_sort_aa[i]] )
            if(this_var == x)
                push!(idx_variable_set[j], i)
            end
        end
    end
    for j in 1:length(variable_set)
        if(length(idx_variable_set[j])>0)
            v_left = minimum(idx_variable_set[j])-0.4
            v_right = maximum(idx_variable_set[j])+0.4
            x_annotate_vec = collect(v_left:0.1:v_right)
            Plots.plot!(x_annotate_vec, (y_annotate-0.4) * ones(length(x_annotate_vec)), c=:black, legend=:false, axis=:false, ticks=:false, lw=2)
            annotate!((v_right+v_left)/2-0.5, y_annotate-0.1, text(variable_set[j], fontsize_reg, :left))
        end
    end
    # ----------------------------------------------------- #
    heatmap_BNB_or_nonBNB_figure = copy(heatmap_both_potent_and_not_potent[:, idx_sort_aa]) 
    mutation_names_BNB_or_nonBNB_figure = mutation_both_potent_and_not_potent[idx_sort_aa]

    annotate!(n_max_both_potent_and_not_potent/2, 0.0, text("- Shared -", fontsize_reg))
    annotate!(6.6, 8.05, text("Increase in \n viral load", fontsize_reg, :red))
    annotate!(-0.5, 10.8, text(L"\textbf{B}", :left, fontsize_label_reg))

    p2tot = Plots.plot(p2, p1, layout=l, 
        size=(L_fig, Int(L_fig_tot*3.5/12)), 
        labelfontsize=fontsize_reg,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
        xrotation=25, 
        margin=my_margin,
    )
    
    return (p2tot, heatmap_BNB_or_nonBNB_figure, mutation_names_BNB_or_nonBNB_figure, heatmap_only_BNB_figure, mutation_names_only_BNB_figure) 
end;


function get_correlation_sumamry(fname_key_human_and_RMs, dir_name, title_fig, fname_key_in, L_fig, n_rand_max)
    csv_raw_CH = DataFrame(CSV.File(dir_name * fname_key_in * "-poly.csv"));
    date_unique = sort(unique(csv_raw_CH.date))
    n_time = length(date_unique)
    # ----------- Start computing variance of sequencies -------------#
    seq_raw = readdlm(dir_name * fname_key_in * "-poly.num");
    fname_key_human_and_RMs_tot = ["RMs";fname_key_human_and_RMs];
    # May be this cell should be written in a function because the similar process need to apply for other conditions. 
    fname_key_human_and_RMs_tot = [[fname_key_in, "RMs"];fname_key_human_and_RMs[2:end]]
    n_H_RMs = length(fname_key_human_and_RMs_tot) 
    cor_mat_smry = zeros(n_H_RMs, n_H_RMs)
    spearman_mat_smry = zeros(n_H_RMs, n_H_RMs)
    cor_mat_smry_std = zeros(n_H_RMs, n_H_RMs)
    spearman_mat_smry_std = zeros(n_H_RMs, n_H_RMs)

    
    for k in 1:n_H_RMs # >2
        for m in (k+1):n_H_RMs
            cor_vec_temp = []; spearman_vec_temp = []
            for n_rand in 1:10
                fitness_ind_tot_k_random_set = []
                fitness_ind_tot_m_random_set = []
                key_k = fname_key_human_and_RMs_tot[k]
                csv_raw_k = DataFrame(CSV.File(dir_name *key_k* "-poly.csv"));
                key_m = fname_key_human_and_RMs_tot[m]
                csv_raw_m = DataFrame(CSV.File(dir_name *key_m* "-poly.csv"));    
                s_in_k = [x ∈ csv_raw_k.mutation ? 
                                  csv_raw_k.s_MPL[ csv_raw_k.mutation .== x ][1] : 0.0    
                                  for x in csv_raw_CH.mutation];
                s_in_m = [x ∈ csv_raw_m.mutation ? 
                                  csv_raw_m.s_MPL[ csv_raw_m.mutation .== x ][1] : 0.0    
                                  for x in csv_raw_CH.mutation];

                time_temp = copy(Int.(seq_raw[:,1]))
                seq_temp = Int.(copy(seq_raw[:,3:end])) .+ 1;
                #"""
                # ---- sampling while keeping the freuqncy ----# 
                seq_temp_shuffled = zeros(size(seq_temp))    
                n_seq_temp = size(seq_temp,1)
                for i in 1:size(seq_temp,2)
                    seq_temp_shuffled[:, i] = shuffle(seq_temp[:, i])
                end
                seq_temp = Int.(seq_temp_shuffled)    
                #"""

                # --- Set selection coefficients --- #
                L_temp = maximum(csv_raw_CH.polymorphic)
                s_extended_k = zeros(L_temp, q_AA)
                s_extended_m = zeros(L_temp, q_AA)
                #
                for i in 1:length(csv_raw_k.s_MPL)
                    i_hxb2_in_k = csv_raw_k.HXB2[i]
                    pro_k = csv_raw_k.PRO[i]
                    s_temp = csv_raw_k.s_MPL[i]    
                    idx_look = (string.(csv_raw_CH.HXB2) .== string(i_hxb2_in_k)) .* (string.(csv_raw_CH.PRO) .== string.(pro_k))
                    if(count(idx_look)>0)
                        if(count(idx_look)==1)
                            i_poly = csv_raw_CH.polymorphic[idx_look][1]
                            a = AA2NUM[pro_k]
                            s_extended_k[i_poly, a] = s_temp
                        end
                    end
                end
                #
                for i in 1:length(csv_raw_m.s_MPL)
                    i_hxb2_in_m = csv_raw_m.HXB2[i]
                    pro_m = csv_raw_m.PRO[i]
                    s_temp = csv_raw_m.s_MPL[i]    
                    idx_look = (string.(csv_raw_CH.HXB2) .== string(i_hxb2_in_m)) .* (string.(csv_raw_CH.PRO) .== string.(pro_m))
                    if(count(idx_look)>0)
                        if(count(idx_look)==1)
                            i_poly = csv_raw_CH.polymorphic[idx_look][1]
                            a = AA2NUM[pro_m]
                            s_extended_m[i_poly, a] = s_temp
                        end
                    end
                end
                # --------------------------------- #
                # --- Obtain fitness values for macaqu sequencies --- #
                for i_t in 1:length(date_unique)
                    idx_t = date_unique[i_t] .== time_temp
                    n_t = count(idx_t)
                    ensemble_t = copy(seq_temp[idx_t, :])
                    f_set_k_temp = []; f_set_m_temp = []
                    for i_temp in 1:n_t
                        x_t = copy(ensemble_t[i_temp, :])
                        f_k_temp = sum([s_extended_k[i, x_t[i]] for i in 1:L_temp])
                        f_m_temp = sum([s_extended_m[i, x_t[i]] for i in 1:L_temp])
                        push!(f_set_k_temp, f_k_temp)
                        push!(f_set_m_temp, f_m_temp)
                    end
                    push!(fitness_ind_tot_k_random_set, 100*copy(float.(f_set_k_temp)))
                    push!(fitness_ind_tot_m_random_set, 100*copy(float.(f_set_m_temp)))
                end;
                fitness_ind_tot_k_random_set = copy(reduce(vcat, fitness_ind_tot_k_random_set))
                fitness_ind_tot_m_random_set = copy(reduce(vcat, fitness_ind_tot_m_random_set));
                cor_this = cor(fitness_ind_tot_k_random_set, fitness_ind_tot_m_random_set)
                spearman_this = corspearman(fitness_ind_tot_k_random_set, fitness_ind_tot_m_random_set)
                push!(cor_vec_temp, cor_this)
                push!(spearman_vec_temp, spearman_this)
            end
            cor_mat_smry[k,m] = mean(cor_vec_temp)
            cor_mat_smry_std[k,m] = std(cor_vec_temp)
            spearman_mat_smry[k,m] = mean(spearman_vec_temp)
            spearman_mat_smry_std[k,m] = std(spearman_vec_temp)
        end
    end

    # ----------------------- Making heatmap plot ------------------------- #
    fontsize_reg = Int(ceil((2*L_fig)/α_gen_sgl * pxl2pt * (1.0/0.8))) # because we rescale these figures
    fontsize_label_reg = Int(ceil((2*L_fig)/α_lbl_sgl * pxl2pt* (1.0/0.8)))
    my_alpha = 0.5

    fname_key_human_and_RMs_tot_label = [[ "CH", "Joint RMs"]; ["RM"*string(i) for i in 1:7]];
    p1 = heatmap(cor_mat_smry, clim=(0.0,1), c=cgrad([:white, :green]), 
        yticks=(1:n_H_RMs, fname_key_human_and_RMs_tot_label) , 
        xticks=(1:n_H_RMs, fname_key_human_and_RMs_tot_label) , 
        labelfontsize=fontsize_reg,
        markerstrokewidth=0,
        grid=:false,
        margin=my_margin,
        legend=:false,
        xrotation=30,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
    )

    for i in 1:n_H_RMs
        for j in (i+1):n_H_RMs
            if(i>1)
                str_this = string(Int(round(cor_mat_smry[i,j]*100)) ) * "±" * string(Int(round(cor_mat_smry_std[i,j]*100)) ) 
                annotate!(j, i, text( str_this, fontsize_reg-3, "Computer Modern"))
            end
            if(i==1)
                str_this = string(Int(round(cor_mat_smry[i,j]*100)) ) * "±" * string(Int(round(cor_mat_smry_std[i,j]*100)) ) 
                annotate!(j, i, text( str_this, fontsize_reg-3, "Computer Modern"))
            end
        end
    end
    annotate!(-1.3, 9.4, text(L"\textbf{A}", :left, fontsize_label_reg))
    annotate!(5, 9, text(title_fig, fontsize_reg))
    annotate!(1.5, 7, text(L"$\mathrm{Pearson's ~~} \bf{r}~\times 100$", :left, fontsize_reg));
    p_out = Plots.plot(p1, size=(L_fig, L_fig))

    return (p_out, cor_mat_smry, cor_mat_smry_std)
end;


function get_selection_vs_time_SHIV505(dir_name, L_fig_tot)
    set_VL_increase = ["N130D", "N279D", "K302N", "Y330H", "N334S", "H417R"];
    fname_potent = ["RM5695", "RM6070", "RM6072"]
    fname_not_potent = ["RM6697", "RM6699", "RM6701", "RM6703"]
    #csv_raw_RMs_merged.mutation
    mutation_potent_set = []
    mutation_not_potent_set = []
    s_potent_set = []
    s_not_potent_set = []
    date_potent_set = []
    date_not_potent_set = []
    for k in 1:length(fname_potent)
        key_RM = fname_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_potent_set, copy(csv_raw.s_MPL))
        push!(date_potent_set, copy(csv_raw.date))
        push!(mutation_potent_set, copy(csv_raw.mutation))
    end;

    mutation_potent_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2], mutation_potent_set[3]);
    mutation_bnAbs_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_potent_union = union(mutation_potent_set[1], mutation_potent_set[2], mutation_potent_set[3]);


    for k in 1:length(fname_not_potent)
        key_RM = fname_not_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(date_not_potent_set, copy(csv_raw.date))
        push!(s_not_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_not_potent_set, copy(csv_raw.mutation))
    end;


    mutation_resistant_bnAbs = ["R166K", "R169K"];
    is_it_bnAbs_RMs_set = []
    for k in 1:3
        bnAbs_resistant_temp = [x ∈ mutation_resistant_bnAbs for x in  mutation_potent_set[k]]
        push!(is_it_bnAbs_RMs_set, copy(bnAbs_resistant_temp))
    end;

    # ------------------------------------- Making figure ----------------------------------- #
    my_ms = 4
    my_alpha = 0.5
    L_fig = Int(ceil(L_fig_tot*0.7*0.5*0.9))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    c_jitter = 10

    k = 1
    p1 = scatter(date_potent_set[k] .+ c_jitter*randn(length(date_potent_set[k])), 100*s_potent_set[k], 
        c=:magenta, 
        label=:false,
        markerstrokecolor=:magenta,
        grid=:false,
        ylabel="\n Selection (%)",
        xlabel=" \n ",
        ylim=(-3, 8),
        alpha=my_alpha, 
        margin=my_margin,
        ms = my_ms,
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        xticks=(0:300:900),
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
    )

    for k in 2:3
        scatter!(date_potent_set[k] .+ c_jitter*randn(length(date_potent_set[k])), 
            100*s_potent_set[k], c=:magenta, label=:false, ms = my_ms,markerstrokecolor="magenta", alpha=my_alpha,)
    end
    for k in 1:3
        if(k==1)
                scatter!(date_potent_set[k][is_it_bnAbs_RMs_set[k]] .+ c_jitter*randn(length(date_potent_set[k][is_it_bnAbs_RMs_set[k]])), 
                100*s_potent_set[k][is_it_bnAbs_RMs_set[k]], 
            c=:green, ms=6, m=:diamond, label=" BnAbs", markerstrokecolor="green",)
        end

        scatter!(date_potent_set[k][is_it_bnAbs_RMs_set[k]] .+ c_jitter*randn(length(date_potent_set[k][is_it_bnAbs_RMs_set[k]])), 
            100*s_potent_set[k][is_it_bnAbs_RMs_set[k]], 
        c=:green, m=:diamond, label=:false, ms = 6, markerstrokecolor="green",)
    end
    annotate!(-200, 8.3, text(L"\textbf{A}", :left, fontsize_label_reg))
    annotate!(400, 8.3, text("Broad breadth", fontsize_reg))
    annotate!(400, -4.3, text("Time mutation was first observed", fontsize_reg));
    #display(p1)

    L_fig = Int(ceil(L_fig_tot*0.3*0.5*0.9)) # one plot taking 30% and another taking 70%, the combined one taking 50% for each, and saling it by 90%
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    myalpha= 0.5

    s_res_dev = [s_potent_set[1];s_potent_set[2];s_potent_set[3]]
    s_res_undev = [s_not_potent_set[1];s_not_potent_set[2];s_not_potent_set[3];s_not_potent_set[4]]

    v1 = 100*s_res_dev; n1=length(v1)
    v2 = 100*s_res_undev; n2=length(v2)

    # Combine them into a single vector y
    y = vcat(v1, v2)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)))

    barwidth = 0.5
    width = 0.4 * barwidth

    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    # jittered x
    p_violin=scatter(x_jitter[1:n1], y[1:n1], c=:magenta, markerstrokewidth=0, label="Broad breadth", 
        foreground_color_legend = nothing,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        ms = my_ms,
        size=(L_fig_tot, Int(ceil(0.5*L_fig_tot))),
        #ylabel=L"$S$ (%)",
        alpha=myalpha,
        #margin=my_margin,
        margin=1mm,
    #    yticks=(-2:2:10),
        ylim=(-3, 8),
    #    m="circle",
        xticks=(1:2, ["Broad", "Narrow"]),
        labelfontsize=fontsize_reg,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:gray, ms = my_ms,markerstrokewidth=0, 
        alpha=myalpha, label="Narrow breadth")

    # jittering lines up with a violin plot
    boxplot!(x[1:n1], y[1:n1], alpha=0.2, c=:magenta, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.2, c=:gray, label=:false);
    annotate!(0.6, 8.3, text(L"\textbf{B}", :left, fontsize_label_reg));
    #Plots.savefig("../fig/CH505/violin_plot_resistance.pdf")

    s_BNB_set = [s_potent_set[1]; s_potent_set[2]; s_potent_set[3]]
    date_BNB_set = [date_potent_set[1]; date_potent_set[2]; date_potent_set[3]]
    s_notBNB_set = [s_not_potent_set[1]; s_not_potent_set[2]; s_not_potent_set[3]; s_not_potent_set[4]]
    date_notBNB_set = [date_not_potent_set[1]; date_not_potent_set[2]; date_not_potent_set[3]; date_not_potent_set[4]]
    #
    date_bnAbs_RMs = [date_potent_set[1][is_it_bnAbs_RMs_set[1]]; 
        date_potent_set[2][is_it_bnAbs_RMs_set[2]];
        date_potent_set[3][is_it_bnAbs_RMs_set[3]]]
    s_bnAbs_RMs = [s_potent_set[1][is_it_bnAbs_RMs_set[1]]; 
        s_potent_set[2][is_it_bnAbs_RMs_set[2]];
        s_potent_set[3][is_it_bnAbs_RMs_set[3]]]

    return(p1, p_violin, s_BNB_set, date_BNB_set, s_notBNB_set, date_notBNB_set, date_bnAbs_RMs, s_bnAbs_RMs)
end;


function get_gene_idx(gen_key, csv_raw_in)
    idx = []
    for i in 1:length(csv_raw_in.gene_fr1)
        myflag=false;
        x1 = csv_raw_in.gene_fr1[i]
        x2 = csv_raw_in.gene_fr2[i]
        x3 = csv_raw_in.gene_fr3[i]
        if(!ismissing(x1)) if(x1==gen_key) myflag = true end end
        if(!ismissing(x2)) if(x2==gen_key) myflag = true end end
        if(!ismissing(x3)) if(x3==gen_key) myflag = true end end
        push!(idx, Bool(myflag))
    end
    return idx
end;

function get_trajectory_plot_CH848(csv_raw_CH848, L_fig_tot)
    L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.7))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl* pxl2pt))

    colection_time = get_collection_time(csv_raw_CH848);
    n_top = sum(abs.(csv_raw_CH848.s_MPL) .> 0.01 )
    #n_top = sum(abs.(csv_raw_CH848.s_MPL) .> 0.022 )

    n_time_max = length(colection_time)
    x1_selected = zeros(n_top, length(colection_time))

    for i_t in 1:n_time_max
        x1_selected[:, i_t] = csv_raw_CH848[1:n_top, Symbol("f_at_" * string( colection_time[i_t]))]
    end;
    
    color_set = []
    alpha_set = []
    for n in 1:n_top
        if(csv_raw_CH848.resist_strain_specific_Abs_CH848[n] || csv_raw_CH848.common_mut_SHIV_CH848[n])
            push!(color_set, "purple"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH848.resist_mut_DH475[n] && csv_raw_CH848.detected_date[n] > 300)
            push!(color_set, "cyan"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH848.resist_mut_DH272[n] && csv_raw_CH848.detected_date[n] > 300)
            push!(color_set, "skyblue"); push!(alpha_set, 1.0)
        elseif(csv_raw_CH848.resist_mut_DH270[n] && csv_raw_CH848.detected_date[n] > 900)
            push!(color_set, "green"); push!(alpha_set, 1.0)
        else
            push!(color_set, "gray"); push!(alpha_set, 1.0)
        end
    end
    
    # ---------- Make a subplot ---------------- # 
    my_alpha=0.5
    x1_set_for_figures_DH270 = [] 
    x1_set_for_figures_DH272 = [] 
    x1_set_for_figures_DH475 = [] 
    x1_set_for_figures_strain_specific = [] 

    for n in sort(collect(1:n_top), rev=true)
        if(color_set[n] == "purple")
                push!(x1_set_for_figures_strain_specific, x1_selected[n, :])
        end
        if(color_set[n] == "skyblue")
                push!(x1_set_for_figures_DH272, x1_selected[n, :])
        end
        if(color_set[n] == "green")
                push!(x1_set_for_figures_DH270, x1_selected[n, :])
        end
        if(color_set[n] == "cyan")
                push!(x1_set_for_figures_DH475, x1_selected[n, :])
        end
    end

    my_alpha=0.5
    p4 = Plots.plot(
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        margin=2mm, 
        lw=3,
        alpha=my_alpha, 
        tickfontsize=fontsize_reg, 
        xticks=:false,
        ylabel=" \n ",
        legend=:false, 
        xlabel="\n "
    )
    annotate!(-250, 0.9, text("Frequency (%)", :left, rotation=90, fontsize_reg))
    annotate!(870, -30, text("Days", :left, fontsize_reg))
    for n in sort(collect(2:n_top), rev=true)
        if(color_set[n] != "gray" && color_set[n] != "red" && color_set[n] != "orange")

            Plots.plot!(time_unique, 100*x1_selected[n, :], c=color_set[n], alpha=my_alpha, legend=:false, lw=3,)
        end
    end
    for x in collect(0:300:1800)
        annotate!(x, +1, text("|", 6))
        if(x%600==0)
            annotate!(x, -14, text(string(x), fontsize_reg))
        end
    end
    

    annotate!(-400, 110, text(L"\textbf{C}", :left, fontsize_label_reg));
    #display(p4)
    return (p4, time_unique, x1_set_for_figures_DH270, x1_set_for_figures_DH272, x1_set_for_figures_DH475, 
        x1_set_for_figures_strain_specific)
end;

function get_selection_vs_time_plot_CH848(csv_raw_CH848, L_fig_tot)

    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    my_ms = 5
    myalpha= 0.5

    idx = csv_raw_CH848.resist_mut_DH270
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 900 # matured DH270 was detectable after 3.5 years after 
    my_ms =5
    p5 = Plots.scatter(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1], 
        c="green", 
        m=:circle,
        ms = my_ms,
        alpha=myalpha,
        label="DH270",
        markerstrokewidth=0,
        subplot=1
    )
    time_DH270_figure = copy(xaxis)
    selection_DH270_figure = copy(csv_raw_CH848.s_MPL[idx])

    idx = csv_raw_CH848.resist_mut_DH272
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 300 # DH272 was first detecable after a year after
    Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
        c="skyblue", 
        label="DH272",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_DH272_figure = copy(xaxis)
    selection_DH272_figure = copy(csv_raw_CH848.s_MPL[idx])

    idx = csv_raw_CH848.resist_mut_DH475
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 300 # DH475 was first detecable after a year after
    Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1]
        ,
        c="cyan", 
        label="DH475",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_DH475_figure = copy(xaxis)
    selection_DH475_figure = copy(csv_raw_CH848.s_MPL[idx])

    idx = (csv_raw_CH848.resist_strain_specific_Abs_CH848 .|| csv_raw_CH848.common_mut_SHIV_CH848)
    xaxis = csv_raw_CH848.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH848.s_MPL[idx],
        c="purple",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        label="Strain-specific Abs",
        markerstrokewidth=0,
        subplot=1,
        #xrotation=20,
        xlim=(0, maximum(time_unique)+30),
    #    xticks=((0:900:1800)),
        xticks=:false,
        xlabel = " ",
    #    ylabel = "Selection (%)",
        ylim=(-3, 10), 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        #legend=:false, 
        legend=(0.7, 0.9), 
        margin=my_margin, 
        legendfontsize=fontsize_reg-1,
        tickfontsize=fontsize_reg,
    );
    time_autologous_figure = copy(xaxis)
    selection_autologous_figure = copy(csv_raw_CH848.s_MPL[idx])

    #annotate!(870, -5.3, text("Time mutation was first observed (days after infection)"))
    for x in collect(0:300:1800)
        annotate!(x, -3, text("|", 6))
    end
    annotate!(-180, 4, text("Selection (%)", fontsize_reg, rotation=90))
    annotate!(870, -4.2, text("Time mutation was first observed \n(days after infection)", fontsize_reg))
    annotate!(-350, 10, text(L"\textbf{A}", :left, fontsize_label_reg));

    return (p5, time_DH270_figure, time_DH272_figure, time_DH475_figure, time_autologous_figure, 
                    selection_DH270_figure, selection_DH272_figure, selection_DH475_figure, selection_autologous_figure ) 
end;

function get_violin_plot_CH848(csv_raw_CH848, L_fig_tot)
    L_fig = Int(ceil(L_fig_tot*0.25*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    myalpha= 0.3
    my_ms = 6

    s_res_DH270 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH270 .* (csv_raw_CH848.detected_date .> 900) ]
    s_res_DH272 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH272 .* (csv_raw_CH848.detected_date .> 300)]
    s_res_DH475 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH475 .* (csv_raw_CH848.detected_date .> 300)]
    s_res_spcfc = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_strain_specific_Abs_CH848 .|| csv_raw_CH848.common_mut_SHIV_CH848];

    v1 = 100*s_res_spcfc; n1=length(v1)
    v2 = 100*s_res_DH475; n2=length(v2)
    v3 = 100*s_res_DH272; n3=length(v3)
    v4 = 100*s_res_DH270; n4=length(v4)

    # Combine them into a single vector y
    y = vcat(v1, v2, v3, v4)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)), fill(3, length(v3)), fill(4, length(v4)))

    barwidth = 0.5
    width = 0.4 * barwidth
    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    # jittered x
    p6=scatter(x_jitter[1:n1], y[1:n1], c=:purple, markerstrokewidth=0, label="Strain-Specific Abs", 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        margin=my_margin, 
        ms = my_ms,
        size=(L_fig, Int(ceil(0.6*600))),
        alpha=myalpha,
        ylim=(-2,10),
        m=:circle,
    #     xticks=(1:3, ["Strain-\n specific", "DH272", "DH270"]), 
        xticks=:false,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg-2,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:cyan, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH475")
    scatter!(x_jitter[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], c=:skyblue, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH272")
    scatter!(x_jitter[(n1+n2+n3+1):end], y[(n1+n2+n3+1):end], c=:green, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH272")

    # jittering lines up with a violin plot
    boxplot!(x[1:n1], y[1:n1], alpha=0.3, c=:purple, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.3, c=:cyan, label=:false)
    boxplot!(x[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], alpha=0.3, c=:skyblue, label=:false)
    boxplot!(x[(n1+n2+n3+1):end], y[(n1+n2+n3+1):end], alpha=0.3, c=:green, label=:false)

    annotate!(0, 10, text(L"\textbf{B}", :left, fontsize_label_reg))
    annotate!(-0.5, -4.6, text("Strain\n specific", :left, fontsize_reg, rotation=45))
    annotate!(1, -4.6, text("DH272", :left, fontsize_reg, rotation=45))
    annotate!(2.5, -4.6, text("DH270", :left, fontsize_reg, rotation=45));
    #Plots.savefig("../fig/CH848/violin_plot_resistance.pdf")

    return (p6, s_res_DH270, s_res_DH272, s_res_DH475, s_res_spcfc) 
end;

# --------------------------------- With mutations that affect glycans --------------------------#
function get_trajectory_plot_CH848_with_glycan(csv_raw_CH848, L_fig_tot)
    L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.7))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl* pxl2pt))

    colection_time = get_collection_time(csv_raw_CH848);
    #n_top = sum(abs.(csv_raw_CH848.s_MPL) .> 0.01 )
    n_top = length(csv_raw_CH848.s_MPL)
    #n_top = sum(abs.(csv_raw_CH848.s_MPL) .> 0.022 )

    n_time_max = length(colection_time)
    x1_selected = zeros(n_top, length(colection_time))

    for i_t in 1:n_time_max
        x1_selected[:, i_t] = csv_raw_CH848[1:n_top, Symbol("f_at_" * string( colection_time[i_t]))]
    end;
    
    color_set = []
    alpha_set = []
    idx_glycan_alters = (csv_raw_CH848.N_linked_glycan_minus_fr3 .> 0 .|| csv_raw_CH848.N_linked_glycan_plus_fr3 .> 0) .* (csv_raw_CH848.nonsynonymous .> 0 )
   
    bool_specific_Abs = []; bool_DH475 = []; bool_DH272 = []; bool_DH270 = []; bool_glycan = []
    for n in 1:n_top
        #flag_temp = true
        #if((csv_raw_CH848.resist_strain_specific_Abs_CH848[n] || csv_raw_CH848.common_mut_SHIV_CH848[n]) * (csv_raw_CH848.nonsynonymous[n]>0)) 
        if((csv_raw_CH848.resist_strain_specific_Abs_CH848[n] || csv_raw_CH848.common_mut_SHIV_CH848[n]) && (csv_raw_CH848.nonsynonymous[n]>0)) 
            #push!(color_set, "purple"); push!(alpha_set, 1.0)
            push!(bool_specific_Abs, true)
            #flag_temp = false
        else
            push!(bool_specific_Abs, false)
        end

        #if((csv_raw_CH848.resist_mut_DH475[n] && csv_raw_CH848.detected_date[n] > 300)  * (csv_raw_CH848.nonsynonymous[n]>0))
        if((csv_raw_CH848.resist_mut_DH475[n] && csv_raw_CH848.detected_date[n] > 300) && (csv_raw_CH848.nonsynonymous[n]>0))
            #push!(color_set, "cyan"); push!(alpha_set, 1.0)
            #flag_temp = false
            push!(bool_DH475, true)
        else
            push!(bool_DH475, false)
        end

        if((csv_raw_CH848.resist_mut_DH272[n] && csv_raw_CH848.detected_date[n] > 300) && (csv_raw_CH848.nonsynonymous[n]>0) )
        #if((csv_raw_CH848.resist_mut_DH272[n]) && (csv_raw_CH848.nonsynonymous[n]>0) )
            #push!(color_set, "skyblue"); push!(alpha_set, 1.0)
            #flag_temp = false
            push!(bool_DH272, true)
        else
            push!(bool_DH272, false)
        end

        if((csv_raw_CH848.resist_mut_DH270[n] && csv_raw_CH848.detected_date[n] > 900) && (csv_raw_CH848.nonsynonymous[n]>0 ) )
        #if((csv_raw_CH848.resist_mut_DH270[n]) && (csv_raw_CH848.nonsynonymous[n]>0 ) )
            #push!(color_set, "green"); push!(alpha_set, 1.0)
            #flag_temp = false
            push!(bool_DH270, true)
        else
            push!(bool_DH270, false)
        end

        #if((idx_glycan_alters[n] && csv_raw_CH848.detected_date[n] > 900) * (csv_raw_CH848.nonsynonymous[n]>0))
        if((idx_glycan_alters[n] ) && (csv_raw_CH848.nonsynonymous[n]>0))
            #push!(color_set, "orange"); push!(alpha_set, 1.0)
            #flag_temp = false
            push!(bool_glycan, true)
        else
            push!(bool_glycan, false)
        end
    end
    
    # ---------- Make a subplot ---------------- # 
    my_alpha=0.5
    x1_set_for_figures_DH270 = [] 
    x1_set_for_figures_DH272 = [] 
    x1_set_for_figures_DH475 = [] 
    x1_set_for_figures_glycan = [] 
    x1_set_for_figures_strain_specific = [] 

    for n in sort(collect(1:n_top), rev=true)
        #if(color_set[n] == "purple")
        if(bool_specific_Abs[n])
                push!(x1_set_for_figures_strain_specific, x1_selected[n, :])
        end
        #if(color_set[n] == "skyblue")
        if(bool_DH272[n])
                push!(x1_set_for_figures_DH272, x1_selected[n, :])
        end
        #if(color_set[n] == "green")
        if(bool_DH270[n])
                push!(x1_set_for_figures_DH270, x1_selected[n, :])
        end
        #if(color_set[n] == "orange")
        if(bool_glycan[n])
            push!(x1_set_for_figures_glycan, x1_selected[n, :])
        end
        #if(color_set[n] == "cyan")
        if(bool_DH475[n])
            push!(x1_set_for_figures_DH475, x1_selected[n, :])
        end
    end

    my_alpha=0.5
    p4 = Plots.plot(
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        margin=2mm, 
        lw=3,
        alpha=my_alpha, 
        tickfontsize=fontsize_reg, 
        xticks=:false,
        ylabel=" \n ",
        legend=:false, 
        xlabel="\n "
    )
    annotate!(-250, 0.9, text("Frequency (%)", :left, rotation=90, fontsize_reg))
    annotate!(870, -30, text("Days", :left, fontsize_reg))
    for n in sort(collect(2:n_top), rev=true)
        #if(color_set[n] != "gray" && color_set[n] != "red") Plots.plot!(time_unique, 100*x1_selected[n, :], c=color_set[n], alpha=my_alpha, legend=:false, lw=3,) end
        if(bool_DH270[n]) 
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:green, alpha=my_alpha, legend=:false, lw=3,) 
        end
        if(bool_DH272[n]) 
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:skyblue, alpha=my_alpha, legend=:false, lw=3,) 
        end
        if(bool_DH475[n]) 
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:cyan, alpha=my_alpha, legend=:false, lw=3,) 
        end
        if(bool_specific_Abs[n]) 
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:purple, alpha=my_alpha, legend=:false, lw=3,) 
        end
        if(bool_glycan[n]) 
            Plots.plot!(time_unique, 100*x1_selected[n, :], c=:orange, alpha=my_alpha, legend=:false, lw=3,) 
        end
    end
    for x in collect(0:300:1800)
        annotate!(x, +1, text("|", 6))
        if(x%600==0)
            annotate!(x, -14, text(string(x), fontsize_reg))
        end
    end
    

    annotate!(-400, 110, text(L"\textbf{C}", :left, fontsize_label_reg));
    #display(p4)
    return (p4, time_unique, x1_set_for_figures_DH270, x1_set_for_figures_DH272, x1_set_for_figures_DH475, 
        x1_set_for_figures_strain_specific, x1_set_for_figures_glycan)
end;

function get_selection_vs_time_plot_CH848_with_glycan(csv_raw_CH848, L_fig_tot)

    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    my_ms = 5
    myalpha= 0.5
    #
    idx = (csv_raw_CH848.N_linked_glycan_minus_fr3 .> 0 .|| csv_raw_CH848.N_linked_glycan_plus_fr3 .> 0) .&& (csv_raw_CH848.nonsynonymous .> 0)
    mutation_name_glycan = csv_raw_CH848.mutants_AA_fr3[idx]
    xaxis = csv_raw_CH848.detected_date[idx] 
    p5 = Plots.scatter(
        xaxis, 100*csv_raw_CH848.s_MPL[idx],
        c="orange", 
        m=:circle,
        ms = my_ms,
        alpha=myalpha-0.2,
        label="Glycan",
        markerstrokewidth=0,
        subplot=1
    )
    time_glycan_figure = copy(xaxis)
    selection_glycan_figure = copy(csv_raw_CH848.s_MPL[idx])
    #
    idx = (csv_raw_CH848.resist_strain_specific_Abs_CH848 .|| csv_raw_CH848.common_mut_SHIV_CH848) .&& (csv_raw_CH848.nonsynonymous .> 0)
    mutation_name_strain_specific = csv_raw_CH848.mutants_AA_fr3[idx]
    xaxis = csv_raw_CH848.detected_date[idx]
    Plots.scatter!(xaxis, 100*csv_raw_CH848.s_MPL[idx],
        c="purple",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        label="Strain-specific Abs",
        markerstrokewidth=0,
        subplot=1)
    time_autologous_figure = copy(xaxis)
    selection_autologous_figure = copy(csv_raw_CH848.s_MPL[idx])
    #   
    idx = csv_raw_CH848.resist_mut_DH475 .&& (csv_raw_CH848.nonsynonymous .> 0 )
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 300 # DH475 was first detecable after a year after
    mutation_name_DH475 = csv_raw_CH848.mutants_AA_fr3[idx][idx1]
    Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
        c="cyan", 
        label="DH475",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_DH475_figure = copy(xaxis[idx1])
    selection_DH475_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])
    #
    idx = csv_raw_CH848.resist_mut_DH272 .&& (csv_raw_CH848.nonsynonymous .> 0 )
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 300 # DH272 was first detecable after a year after
    mutation_name_DH272 = csv_raw_CH848.mutants_AA_fr3[idx][idx1]
    Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
        c="skyblue", 
        label="DH272",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        markerstrokewidth=0,
        subplot=1)
    time_DH272_figure = copy(xaxis[idx1])
    selection_DH272_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])
    #
    idx = csv_raw_CH848.resist_mut_DH270 .&& (csv_raw_CH848.nonsynonymous .> 0 )
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = xaxis .> 900 # matured DH270 was detectable after 3.5 years after 
    mutation_name_DH270 = csv_raw_CH848.mutants_AA_fr3[idx][idx1]
    my_ms =5
    Plots.scatter!(xaxis, 
        xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1], 
        c="green",
        m=:circle,
        alpha=myalpha,
        ms = my_ms,
        label="DH270",
        markerstrokewidth=0,
        subplot=1,
        #xrotation=20,
        xlim=(0, maximum(time_unique)+30),
    #    xticks=((0:900:1800)),
        xticks=:false,
        xlabel = " ",
    #    ylabel = "Selection (%)",
        ylim=(-3, 10), 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        #legend=:false, 
        legend=(0.7, 0.9), 
        margin=my_margin, 
        legendfontsize=fontsize_reg-1,
        tickfontsize=fontsize_reg,
    );
    time_DH270_figure = copy(xaxis[idx1])
    selection_DH270_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])

    

    #annotate!(870, -5.3, text("Time mutation was first observed (days after infection)"))
    for x in collect(0:300:1800)
        annotate!(x, -3, text("|", 6))
    end
    annotate!(-180, 4, text("Selection (%)", fontsize_reg, rotation=90))
    annotate!(870, -4.2, text("Time mutation was first observed \n(days after infection)", fontsize_reg))
    annotate!(-350, 10, text(L"\textbf{A}", :left, fontsize_label_reg));

    return (p5, time_DH270_figure, time_DH272_figure, time_DH475_figure, time_autologous_figure, time_glycan_figure,  
                    selection_DH270_figure, selection_DH272_figure, selection_DH475_figure, selection_autologous_figure, selection_glycan_figure, 
                    mutation_name_DH270, mutation_name_DH272, mutation_name_DH475, mutation_name_strain_specific, mutation_name_glycan) 
end;

function get_violin_plot_CH848_with_glycan(csv_raw_CH848, L_fig_tot)
    L_fig = Int(ceil(L_fig_tot*0.25*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    myalpha= 0.3
    my_ms = 6

    s_res_DH270 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH270 .&& (csv_raw_CH848.detected_date .> 900) .&& (csv_raw_CH848.nonsynonymous .> 0 )]
    s_res_DH272 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH272 .&& (csv_raw_CH848.detected_date .> 300) .&& (csv_raw_CH848.nonsynonymous .> 0 )]
    s_res_DH475 = csv_raw_CH848.s_MPL[csv_raw_CH848.resist_mut_DH475 .&& (csv_raw_CH848.detected_date .> 300) .&& (csv_raw_CH848.nonsynonymous .> 0 )]
    s_res_spcfc = csv_raw_CH848.s_MPL[(csv_raw_CH848.resist_strain_specific_Abs_CH848 .|| csv_raw_CH848.common_mut_SHIV_CH848) .&& (csv_raw_CH848.nonsynonymous .> 0 )];
    s_res_glycan = csv_raw_CH848.s_MPL[(csv_raw_CH848.N_linked_glycan_minus_fr3 .> 0 .|| csv_raw_CH848.N_linked_glycan_plus_fr3 .> 0) .&& (csv_raw_CH848.nonsynonymous .> 0 )];

    v1 = 100*s_res_glycan; n1=length(v1)
    v2 = 100*s_res_spcfc; n2=length(v2)
    v3 = 100*s_res_DH475; n3=length(v3)
    v4 = 100*s_res_DH272; n4=length(v4)
    v5 = 100*s_res_DH270; n5=length(v5)

    # Combine them into a single vector y
    y = vcat(v1, v2, v3, v4, v5)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)), fill(3, length(v3)), fill(4, length(v4)), fill(5, length(v5)))

    barwidth = 0.5
    width = 0.4 * barwidth
    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    # jittered x
    p6=scatter(x_jitter[1:n1], y[1:n1], c=:orange, markerstrokewidth=0, label="Glycan", 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        margin=my_margin, 
        ms = my_ms,
        size=(L_fig, Int(ceil(0.6*600))),
        alpha=myalpha,
        ylim=(-3,10),
        m=:circle,
    #     xticks=(1:3, ["Strain-\n specific", "DH272", "DH270"]), 
        xticks=:false,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg-2,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:purple, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="Specific")
    scatter!(x_jitter[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], c=:cyan, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH475")
    scatter!(x_jitter[(n1+n2+n3+1):(n1+n2+n3+n4)], y[(n1+n2+n3+1):(n1+n2+n3+n4)], c=:skyblue, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH272")
    scatter!(x_jitter[(n1+n2+n3+n4+1):end], y[(n1+n2+n3+n4+1):end], c=:green, m=:circle, ms = my_ms, markerstrokewidth=0, 
        alpha=myalpha, label="DH270")

    # jittering lines up with a violin plot
    boxplot!(x[1:n1], y[1:n1], alpha=0.3, c=:orange, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.3, c=:purple, label=:false)
    boxplot!(x[(n1+n2+1):(n1+n2+n3)], y[(n1+n2+1):(n1+n2+n3)], alpha=0.3, c=:cyan, label=:false)
    boxplot!(x[(n1+n2+n3+1):(n1+n2+n3+n4)], y[(n1+n2+n3+1):(n1+n2+n3+n4)], alpha=0.3, c=:skyblue, label=:false)
    boxplot!(x[(n1+n2+n3+n4+1):end], y[(n1+n2+n3+n4+1):end], alpha=0.3, c=:green, label=:false)

    annotate!(0, 10, text(L"\textbf{B}", :left, fontsize_label_reg))
    annotate!(-0.6, -5.5, text("Glycan", :left, fontsize_reg, rotation=45))
    annotate!(0.5, -5.5, text("Specific", :left, fontsize_reg, rotation=45))
    annotate!(1.5, -5.5, text("DH475", :left, fontsize_reg, rotation=45))
    annotate!(2.5, -5.5, text("DH272", :left, fontsize_reg, rotation=45))
    annotate!(3.5, -5.5, text("DH270", :left, fontsize_reg, rotation=45));
    #Plots.savefig("../fig/CH848/violin_plot_resistance.pdf")

    return (p6, s_res_DH270, s_res_DH272, s_res_DH475, s_res_spcfc, s_res_glycan) 
end;
# ------------------------------------------------------------------------------------------#

function get_trajectory_plot_CH848_reduced(csv_raw_CH848, L_fig_tot)    
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    coeff_th = 0.015

    colection_time = get_collection_time(csv_raw_CH848);
    n_top = sum(csv_raw_CH848.s_MPL .> coeff_th)
    n_time_max = length(colection_time)
    x1_selected = zeros(n_top, length(colection_time))

    for i_t in 1:n_time_max
        x1_selected[:, i_t] = csv_raw_CH848[1:n_top, Symbol("f_at_" * string( colection_time[i_t]))]
    end;

    gene_set_temp = [gen_to_frame[x] for x in gene_set_unique];
    color_set = [gen_to_color[x] for x in gene_set_unique];
    color_all_set = []
    alpha_set = []
    n_gene_set = length(gene_set_unique)
    idx_to_take = []
    for n in 1:n_top
        for k in 1:n_gene_set
            if(length(gene_set_temp[k])==1)
                if(string(csv_raw_CH848[n, Symbol(gene_set_temp[k][1])])==gene_set_unique[k])
                    if(string(csv_raw_CH848[n, Symbol(gene_set_temp[k][1])]) ∈ ["nef", "env"])
                        push!(color_all_set, color_set[k]); 
                    else
                        push!(color_all_set, "gray"); 
                    end
                    push!(alpha_set, 0.5)
                    push!(idx_to_take, n)
                end
            else
                if(string(csv_raw_CH848[n, Symbol(gene_set_temp[k][1])])==gene_set_unique[k]||
                   string(csv_raw_CH848[n, Symbol(gene_set_temp[k][2])])==gene_set_unique[k] )
                    push!(color_all_set, "gray"); 
                    push!(alpha_set, 0.5)
                    push!(idx_to_take, n)
                end
            end
        end
    end
    # Update because some doesn't belowng to any gene
    x1_selected = copy(x1_selected[idx_to_take, :] );
    n_top = length(color_all_set);

    # ---------- Make a subplot ---------------- # 
    my_alpha=0.5
    x1_set_for_figures_nef = [] 
    x1_set_for_figures_env = [] 
    x1_set_for_figures_else = [] 

    for n in sort(collect(1:n_top), rev=true)
        if(color_all_set[n] == "orange")
                push!(x1_set_for_figures_nef, x1_selected[n, :])
        end
        if(color_all_set[n] == "green")
                push!(x1_set_for_figures_env, x1_selected[n, :])
        end
        if(color_all_set[n] == "gray")
                push!(x1_set_for_figures_else, x1_selected[n, :])
        end
    end

    my_alpha = 0.5
    p4 = Plots.plot(time_unique, 100*x1_selected[1, :], c=color_all_set[1], 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        #xlim=(0,maximum(time_unique)+30),
        grid=:false,
        margin=my_margin, 
        lw=3,
        alpha=my_alpha,
        tickfontsize=fontsize_reg,
        xticks=:false,
    #    xticks=((0:300:1800)),    
        ylabel=" \n ",
        legend=:false, 
        xlabel=" \n "
    )
    annotate!(-300, 0.3, text("Frequency (%)", :left, rotation=90, fontsize_reg))
    annotate!(900, -30, text("Days", fontsize_reg))
    for n in sort(collect(2:n_top), rev=true)
        Plots.plot!(time_unique, 100*x1_selected[n, :], c=color_all_set[n], alpha=alpha_set[n], legend=:false, lw=3,)
    end
    for x in collect(0:300:1700)    
        annotate!(x, +1, text("|", 6))
        if(x%600==0)
            annotate!(x, -14, text(string(x), fontsize_reg))
        end
    end
    annotate!(-400, 120, text(L"\textbf{F}", :left, fontsize_label_reg));
    # return (p4, time_unique, x1_set_for_figures_CTL, x1_set_for_figures_CH103, x1_set_for_figures_CH235, x1_set_for_figures_strain_specific)
    #display(p4)

    return (p4, time_unique, x1_set_for_figures_nef, x1_set_for_figures_env, x1_set_for_figures_else)
end;

function get_selection_vs_time_plot_CH848_reduced(csv_raw_CH848, L_fig_tot)

    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    color_set = [gen_to_color[x] for x in gene_set_unique];
    gene_set_temp = [gen_to_frame[x] for x in gene_set_unique];
    markershape_set = ["circle" for _ in 1:8];

    coeff_th = 0.015
    my_alpha= 0.5
    my_ms = 5
    k = 1
    idx = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][1])]]
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = abs.(100*csv_raw_CH848.s_MPL[idx]) .> coeff_th*100
    p5 = Plots.scatter(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
        c=color_set[k], 
        label=label_gene_set[k],
        m=Symbol(markershape_set[k]),
        ms = my_ms, 
        subplot=1,
        alpha=my_alpha,
        xlim=(0, maximum(time_unique)+30),
        #xticks=((0:900:1800)),
        xticks=:false,
        xlabel = " ",
        ylabel = " ",
        ylim=(-3, 13), 
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        markerstrokewidth=0,
        grid=:false,
        #legend=:false, 
        legend=(0.8, 0.9), 
        margin=my_margin, 
        legendfontsize=fontsize_reg-1,
        tickfontsize=fontsize_reg,
    )
    time_nef_figure = copy(xaxis[idx1])
    selection_nef_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])

    # --- Env ---- #
    k = 2
    idx = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][1])]]
    if(size(gene_set_temp[k],1) > 1)
        idx_temp1 = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][1])]] 
        idx_temp2 = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][2])]] 
        idx = idx_temp1 .|| idx_temp2
    end
    xaxis = csv_raw_CH848.detected_date[idx]
    idx1 = abs.(100*csv_raw_CH848.s_MPL[idx]) .> coeff_th*100
    Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
        c=color_set[k], 
        label=label_gene_set[k],
        alpha=my_alpha,
        m=Symbol(markershape_set[k]),
        ms = my_ms,
        margin=my_margin, 
        markerstrokewidth=0,
        subplot=1,
    )
    time_env_figure = copy(xaxis[idx1])
    selection_env_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])

    time_els_figure = []
    selection_els_figure = []
    for k in 3:length(gene_set_unique)-1
        if(length(gene_set_temp[k])==1)
            idx = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][1])]]
        end
        if(size(gene_set_temp[k],1) > 1)
            idx_temp1 = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][1])]] 
            idx_temp2 = [!ismissing(x) ? x == gene_set_unique[k] : false for x in csv_raw_CH848[:, Symbol(gene_set_temp[k][2])]] 
            idx = idx_temp1 .|| idx_temp2
        end
        xaxis = csv_raw_CH848.detected_date[idx]
        idx1 = abs.(100*csv_raw_CH848.s_MPL[idx]) .> coeff_th*100

        if(k==3)
            Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
                c=:gray, 
                label="Else", alpha=my_alpha, m=Symbol(markershape_set[k]),
                ms = my_ms, markerstrokewidth=0, subplot=1,)
            time_els_figure = copy(xaxis[idx1])
            selection_els_figure = copy(csv_raw_CH848.s_MPL[idx][idx1])
        else
            Plots.scatter!(xaxis[idx1], 100*csv_raw_CH848.s_MPL[idx][idx1],
                c=:gray, 
                label=:false, alpha=my_alpha, m=Symbol(markershape_set[k]),
                ms = my_ms, markerstrokewidth=0, subplot=1,)
        end        
    end

    for x in collect(0:300:1700)
        annotate!(x, -3, text("|", 6))
    end
    annotate!(-230, 4, text("Selection (%)", fontsize_reg, rotation=90))
    annotate!(870, -4.6, text("Time mutation was first observed \n(days after infection)", fontsize_reg))
    annotate!(-400, 13, text(L"\textbf{D}", :left, fontsize_label_reg));


    return (p5, time_nef_figure, time_env_figure, time_els_figure, 
                selection_nef_figure, selection_env_figure, selection_els_figure)
end;

function get_violin_plot_CH848_reduced(csv_raw_CH848, L_fig_tot)
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    s_set = []
    for i in 1:length(gene_set_unique)
        gen_key = gene_set_unique[i]
        idx = Bool.(get_gene_idx(gen_key, csv_raw_CH848))
        s_in = 100*csv_raw_CH848.s_MPL[idx]
        idx_cut = abs.(s_in) .> 0.1
        push!(s_set, copy(s_in)[idx_cut])
    end
    s_set = copy([s_set[1], s_set[2], 
        vcat(s_set[3], s_set[4], s_set[5], 
             s_set[6], s_set[7], s_set[8]) ])

    n_set = [length(s_temp) for s_temp in s_set]
    n_set = [n_set[1], n_set[2], sum(n_set[3:end])]
    n_cum = cumsum(n_set)
    # Combine them into a single vector y
    y = vcat(s_set[1], s_set[2], s_set[3])
    x = vcat(fill(1, length(s_set[1])), fill(2, length(s_set[2])), fill(3, length(s_set[3])) );

    L_fig = Int(ceil(L_fig_tot*0.25*0.5))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    color_set = [gen_to_color[x] for x in gene_set_unique];
    markershape_set = ["circle" for _ in 1:8];

    myalpha= 0.3
    my_ms = 6
    barwidth = 0.5
    width = 0.4 * barwidth
    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))

    s_res_nef = copy(x[1:n_cum[1]])
    s_res_env = copy(x[(n_cum[1]+1):n_cum[2]])
    s_res_els = copy(x[(n_cum[2]+1):end])

    # jittered x
    p6 = scatter(x_jitter[1:n_cum[1]], y[1:n_cum[1]], c=color_set[1], markerstrokewidth=0.0,markerstrokealpha=0.0, 
        label=label_gene_set[1], 
        foreground_color_legend = nothing,
        markerstrokecolor=color_set[1], 
        labelfontsize=fontsize_reg,
        grid=:false,
        alpha=0.3, 
        ms=5,
        legend=:false, 
        yticks=:false, 
        yaxis=:false, 
        m=:circle,
        margin=my_margin, 
        #xticks=(1:3, ["Nef", "Env", "Else"]), 
        xticks=:false,
        legendfontsize=fontsize_reg,
        tickfontsize=fontsize_reg,
    )
    k = 1
    scatter!(x_jitter[(n_cum[k]+1):n_cum[k+1]], y[(n_cum[k]+1):n_cum[k+1]], 
        c=color_set[k+1], m=Symbol(markershape_set[k+1]), ms=5,
        markerstrokewidth=0.0, markerstrokealpha=0.0, markerstrokecolor=color_set[k+1], 
        alpha=0.3, label=label_gene_set[k+1], )

    k = 2
    scatter!(x_jitter[(n_cum[k]+1):end], y[(n_cum[k]+1):end], 
        c=:gray, m=:circle, ms=5,
        markerstrokewidth=0.0, markerstrokealpha=0.0, markerstrokecolor=:gray, 
        alpha=0.3, label="Else" )

    # jittering lines up with a violin plot
    boxplot!(x[1:n_cum[1]], y[1:n_cum[1]], alpha=0.3, c=color_set[1], label=:false, subplot=1, markerstrokewidth=0.0, )
    k = 1
    boxplot!(x[(n_cum[k]+1):n_cum[k+1]], y[(n_cum[k]+1):n_cum[k+1]], alpha=0.3, c=color_set[k+1], label=:false, markerstrokewidth=0.0, )
    k = 2
    boxplot!(x[(n_cum[k]+1):end], y[(n_cum[k]+1):end], alpha=0.3, c=:gray, label=:false, markerstrokewidth=0.0, )

    annotate!(0, 13, text(L"\textbf{E}", :left, fontsize_label_reg))
    annotate!(1, -5, text("Nef", fontsize_reg))
    annotate!(2, -5, text("Env", fontsize_reg))
    annotate!(3, -5, text("Else", fontsize_reg));
    #Plots.savefig("../fig/CH848/violin_plot_resistance.pdf")
    return (p6, s_res_nef, s_res_env, s_res_els)
end;

function get_heatmap_selected_mutations_CH848(dir_name, L_fig_tot)
    fname_potent = ["RM6163", "RM6167"]
    fname_not_potent = ["RM6700", "RM6713", "RM6714", "RM6720"]
    mutation_potent_set = []
    mutation_not_potent_set = []
    s_potent_set = []
    s_not_potent_set = []
    for k in 1:length(fname_potent)
        key_RM = fname_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_potent_set, copy(csv_raw.mutation))
    end;

    mutation_potent_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_bnAbs_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_potent_union = union(mutation_potent_set[1], mutation_potent_set[2]);


    for k in 1:length(fname_not_potent)
        key_RM = fname_not_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_not_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_not_potent_set, copy(csv_raw.mutation))
    end;
    mutation_not_potent_intersect = intersect(mutation_not_potent_set[1], mutation_not_potent_set[2], 
        mutation_not_potent_set[3], mutation_not_potent_set[4]);

    mutation_not_potent_union = union(mutation_not_potent_set[1], mutation_not_potent_set[2], 
        mutation_not_potent_set[3], mutation_not_potent_set[4]);
    
    # --- Get common mutations that observed in potent and not-potent RMs ---#\\\
    # This condition is too strong not intersect and 
    mutation_both_potent_and_not_potent = intersect(mutation_potent_union, mutation_not_potent_union);
    # --- Get mutations that are common in potent group but never seen in nont-potent group --- #
    mutation_only_in_potent = setdiff(mutation_potent_intersect, mutation_not_potent_union);
    mutation_only_in_bnAbs = setdiff(mutation_bnAbs_intersect, mutation_not_potent_union);
    #mutation_only_in_not_potent = setdiff(mutation_not_potent_intersect, mutation_potent_intersect);
    mutation_both_potent_and_not_potent;
    mutation_only_in_potent
    n_max_both_potent_and_not_potent = 15
    n_max_only_in_potent = 5

    heatmap_both_potent_and_not_potent = -0.09 * ones(length(fname_key_human_and_RMs), n_max_both_potent_and_not_potent)
    heatmap_only_in_potent = -0.09 * ones(length(fname_key_human_and_RMs), n_max_only_in_potent);

    for k in 1:n_RMs_max
        if(k>1)
            key_RM = fname_key_human_and_RMs[k]
            csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        end
        if(k==1)
            key_RM = fname_key_human_and_RMs[k]
            csv_raw = DataFrame(CSV.File(dir_name * "RMs-poly.csv"));
        end
        for n in 1:n_max_both_potent_and_not_potent
            mut_temp = mutation_both_potent_and_not_potent[n]
            if(mut_temp ∈ csv_raw.mutation)
                idx = mut_temp .== csv_raw.mutation
                heatmap_both_potent_and_not_potent[n_RMs_max+1-k, n] = csv_raw.s_MPL[idx][1]
            end
        end
        for n in 1:n_max_only_in_potent
            mut_temp = mutation_only_in_potent[n]
            if(mut_temp ∈ csv_raw.mutation)
                idx = mut_temp .== csv_raw.mutation
                heatmap_only_in_potent[n_RMs_max+1-k, n] = csv_raw.s_MPL[idx][1]
            end
        end
    end

    # ====================== Making Figures ========================= #  
    L_fig = Int(0.6*L_fig_tot)
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))

    i_aa  = extract_integer.(mutation_only_in_potent[1:n_max_only_in_potent])
    idx_sort_aa = sortperm(i_aa)
    l = @layout[b{0.72w, 0.9h} a{0.24w}]
    p1 = heatmap(heatmap_only_in_potent[:, idx_sort_aa], clim=(-0.09,0.09), c=cgrad([:gray80, :white, :red]), 
        legend=:false, xmirror = true, 
        xticks=(1:n_max_only_in_potent, mutation_only_in_potent[1:n_max_only_in_potent][idx_sort_aa]) , 
        xtickfontcolor=:white,
        title=" \n \n \n",
        yaxis=:false, 
    )
    annotate!(n_max_only_in_potent/2+0.5, 0.0, text("- Only potent -", fontsize_reg))
    for i in 1:n_max_only_in_potent
        annotate!(i-0.07, 8-0.37, text(mutation_only_in_potent[idx_sort_aa[i]], :black, fontsize_reg-4, :left, rotation=90))
    end
    heatmap_only_BNB_figure = copy(heatmap_only_in_potent[:, idx_sort_aa])
    mutation_names_only_BNB_figure = mutation_only_in_potent[1:n_max_only_in_potent][idx_sort_aa]

    # ----------------- Annotate Region ------------------ #
    y_annotate = 10.0
    variable_set = ["V1", "V2", "LD", "V3", "V4", "V5"];
    idx_variable_set = [[] for _ in 1:length(variable_set)]
    for j in 1:length(variable_set)
        x = variable_set[j]
        for i in 1:length(i_aa[idx_sort_aa])
            this_var = get_variable_regions(i_aa[idx_sort_aa[i]] )
            if(this_var == x)
                push!(idx_variable_set[j], i)
            end
        end
    end
    for j in 1:length(variable_set)
        if(length(idx_variable_set[j])>0)
            v_left = minimum(idx_variable_set[j])-0.4
            v_right = maximum(idx_variable_set[j])+0.4
            x_annotate_vec = collect(v_left:0.1:v_right)
            Plots.plot!(x_annotate_vec, (y_annotate-0.5) * ones(length(x_annotate_vec)), c=:black, legend=:false, axis=:false, ticks=:false, lw=2)
            annotate!((v_right+v_left)/2-0.5, y_annotate-0.2, text(variable_set[j], fontsize_reg, :left))
        end
    end
    # ----------------------------------------------------- #

    i_aa  = extract_integer.(mutation_both_potent_and_not_potent[1:n_max_both_potent_and_not_potent])
    idx_sort_aa = sortperm(i_aa)
    idx_sort_aa = copy(idx_sort_aa[i_aa[idx_sort_aa] .!= 375])
    p2 = heatmap(heatmap_both_potent_and_not_potent[:, idx_sort_aa], clim=(-0.09,0.09), c=cgrad([:gray80, :white, :red]), 
        legend=:false, xmirror = true, 
        yticks=(1:n_H_RMs, reverse(["RMs"; fname_key_human_and_RMs[2:end]])) , 
        xticks=false,
        ylabel=" \n ", 
        ytickfontcolor=:white
    )
    for i in 1:length(idx_sort_aa)
        annotate!(i-0.07, 8-0.37, text(mutation_both_potent_and_not_potent[idx_sort_aa[i]], fontsize_reg-4, :left, rotation=90))
    end

    for k in 1:n_H_RMs
        key_RM = fname_key_human_and_RMs[k]
        this_c = "black"
        if(key_RM ∈ fname_dev) this_c = "magenta" end
    #    if(key_RM == "RM6072") this_c = "pink" end
        if(k==1) key_RM = "Joint" end
        if(k!=1) key_RM = "RM"*string(k-1) end
        annotate!(-0.05, n_H_RMs+1-k, text(key_RM, Symbol(this_c), fontsize_reg))
    end
    # ----------------- Annotate Region ------------------ #
    y_annotate = 10.0
    idx_variable_set = [[] for _ in 1:length(variable_set)]
    for j in 1:length(variable_set)
        x = variable_set[j]
        for i in 1:length(i_aa[idx_sort_aa])
            this_var = get_variable_regions(i_aa[idx_sort_aa[i]] )
            if(this_var == x)
                push!(idx_variable_set[j], i)
            end
        end
    end
    for j in 1:length(variable_set)
        if(length(idx_variable_set[j])>0)
            v_left = minimum(idx_variable_set[j])-0.4
            v_right = maximum(idx_variable_set[j])+0.4
            x_annotate_vec = collect(v_left:0.1:v_right)
            Plots.plot!(x_annotate_vec, (y_annotate-0.5) * ones(length(x_annotate_vec)), c=:black, legend=:false, axis=:false, ticks=:false, lw=2)
            annotate!((v_right+v_left)/2-0.5, y_annotate-0.2, text(variable_set[j], fontsize_reg, :left))
        end
    end
    # ----------------------------------------------------- #
    heatmap_BNB_or_nonBNB_figure = copy(heatmap_both_potent_and_not_potent[:, idx_sort_aa]) 
    mutation_names_BNB_or_nonBNB_figure = mutation_both_potent_and_not_potent[idx_sort_aa]

    annotate!(n_max_both_potent_and_not_potent/2, 0.0, text("- Shared -", fontsize_reg))
    annotate!(-0.5, 9.5, text(L"\textbf{B}", :left, fontsize_label_reg))
    p2tot = Plots.plot(p2, p1, layout=l, 
        size=(L_fig, Int(L_fig_tot*3.5/12)), 
        labelfontsize=fontsize_reg,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
        xrotation=25, 
        margin=1mm,
    )
    return (p2tot, heatmap_BNB_or_nonBNB_figure, mutation_names_BNB_or_nonBNB_figure, heatmap_only_BNB_figure, mutation_names_only_BNB_figure) 
end;

function get_selection_vs_time_SHIV848(dir_name, L_fig_tot)
    fname_potent = ["RM6163", "RM6167"]
    fname_not_potent = ["RM6700", "RM6713", "RM6714", "RM6720"]
    #csv_raw_RMs_merged.mutation
    mutation_potent_set = []
    mutation_not_potent_set = []
    s_potent_set = []
    s_not_potent_set = []
    date_potent_set = []
    date_not_potent_set = []
    for k in 1:length(fname_potent)
        key_RM = fname_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(s_potent_set, copy(csv_raw.s_MPL))
        push!(date_potent_set, copy(csv_raw.date))
        push!(mutation_potent_set, copy(csv_raw.mutation))
    end;

    mutation_potent_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_bnAbs_intersect = intersect(mutation_potent_set[1], mutation_potent_set[2]);
    mutation_potent_union = union(mutation_potent_set[1], mutation_potent_set[2]);


    for k in 1:length(fname_not_potent)
        key_RM = fname_not_potent[k]
        csv_raw = DataFrame(CSV.File(dir_name * key_RM * "-poly.csv"));
        push!(date_not_potent_set, copy(csv_raw.date))
        push!(s_not_potent_set, copy(csv_raw.s_MPL))
        push!(mutation_not_potent_set, copy(csv_raw.mutation))
    end;


    #, "I326T" appeared to be resistant but this mutaiton wasn't observed before w80, when bnAbs were detectable
    mutation_resistant_bnAbs = ["N301D", "D325N"];

    is_it_bnAbs_RMs_set = []
    for k in 1:2
        bnAbs_resistant_temp = [x ∈ mutation_resistant_bnAbs for x in  mutation_potent_set[k]]
        push!(is_it_bnAbs_RMs_set, copy(bnAbs_resistant_temp))
    end;

    my_ms = 4
    my_alpha = 0.5
    L_fig_tot = 1200
    L_fig = Int(ceil(L_fig_tot*0.7*0.5*0.9))
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    fontsize_label_reg = Int(ceil(L_fig_tot/α_lbl_sgl * pxl2pt))
    c_jitter = 10

    k = 1
    p1 = scatter(date_potent_set[k] .+ c_jitter*randn(length(date_potent_set[k])), 100*s_potent_set[k], 
        c=:magenta, 
        label=:false,
        #title="Broad breadth",
        markerstrokecolor=:magenta,
        grid=:false,
        ylabel="\n Selection (%)",
        xlabel=" \n ",
        ylim=(-2, 6),
        #titlefontsize=fontsize_reg,
        ms = my_ms,
        alpha=my_alpha,
        margin=my_margin,
        foreground_color_legend = nothing,
        labelfontsize=fontsize_reg,
        xticks=(0:300:900),
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
    )
    for k in 2:2
        scatter!(date_potent_set[k] .+ c_jitter*randn(length(date_potent_set[k])), 
            100*s_potent_set[k], c=:magenta, label=:false, ms = my_ms,markerstrokecolor="magenta",)
    end
    for k in 1:2
        if(k==1)
                scatter!(date_potent_set[k][is_it_bnAbs_RMs_set[k]] .+ c_jitter*randn(length(date_potent_set[k])), 
                100*s_potent_set[k][is_it_bnAbs_RMs_set[k]], 
            c=:green, ms=6, m=:diamond, label=" BnAbs", markerstrokecolor="green",)
        end

        scatter!(date_potent_set[k][is_it_bnAbs_RMs_set[k]] .+ c_jitter*randn(length(date_potent_set[k])), 
            100*s_potent_set[k][is_it_bnAbs_RMs_set[k]], 
        c=:green, m=:diamond, label=:false, ms = 6, markerstrokecolor="green",)
    end
    annotate!(-300, 6.2, text(L"\textbf{C}", :left, fontsize_label_reg))
    annotate!(400, 6.2, text("Broad breadth", fontsize_reg))
    annotate!(450, -3.1, text("Time mutation was first observed", fontsize_reg));

    # =============================== Making Figure ====================================#
    L_fig = Int(ceil(L_fig_tot*0.3*0.5*0.9)) # one plot taking 30% and another taking 70%, the combined one taking 50% for each, and saling it by 90%
    fontsize_reg = Int(ceil(L_fig_tot/α_gen_sgl * pxl2pt))
    myalpha= 0.5

    s_res_dev = [s_potent_set[1];s_potent_set[2]]
    s_res_undev = [s_not_potent_set[1];s_not_potent_set[2];s_not_potent_set[3];s_not_potent_set[4]]

    v1 = 100*s_res_dev; n1=length(v1)
    v2 = 100*s_res_undev; n2=length(v2)

    # Combine them into a single vector y
    y = vcat(v1, v2)
    x = vcat(fill(1, length(v1)), fill(2, length(v2)))

    barwidth = 0.5
    width = 0.4 * barwidth

    ngroups = length(unique(x))

    k = Array{UnivariateKDE}(undef, ngroups)
    for i in 1:ngroups
        k[i] = kde(y[x .== i])
    end
    max_dens = map(x -> maximum(x.density), k)

    x_jitter = x .+ width./max_dens[x] .* rand(length(x)) .* pdf.(k[x], y) .* rand([-1, 1], length(x))
    p_violin=scatter(x_jitter[1:n1], y[1:n1], c=:magenta, markerstrokewidth=0, label="Broad breadth", 
        foreground_color_legend = nothing,
        grid=:false,
        yticks=:false,
        yaxis=:false,
        legend=:false,
        ms = my_ms,
        size=(L_fig_tot, Int(ceil(0.5*L_fig_tot))),
        alpha=myalpha,
        ylim=(-2, 6),
        margin=my_margin,
        xticks=(1:2, ["Broad", "Narrow"]),
        labelfontsize=fontsize_reg,
        legendfontsize=fontsize_reg, 
        tickfontsize=fontsize_reg,
    )
    scatter!(x_jitter[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], c=:gray, ms = my_ms,markerstrokewidth=0, 
        alpha=myalpha, label="Narrow breadth")

    boxplot!(x[1:n1], y[1:n1], alpha=0.2, c=:magenta, label=:false)
    boxplot!(x[(n1+1):(n1+n2)], y[(n1+1):(n1+n2)], alpha=0.2, c=:gray, label=:false);
    annotate!(0.5, 6.2, text(L"\textbf{D}", :left, fontsize_label_reg));

    s_BNB_set = [s_potent_set[1]; s_potent_set[2]]
    date_BNB_set = [date_potent_set[1]; date_potent_set[2]]
    s_notBNB_set = [s_not_potent_set[1]; s_not_potent_set[2]; s_not_potent_set[3]; s_not_potent_set[4]]
    date_notBNB_set = [date_not_potent_set[1]; date_not_potent_set[2]; date_not_potent_set[3]; date_not_potent_set[4]]
    #
    date_bnAbs_RMs = [date_potent_set[1][is_it_bnAbs_RMs_set[1]]; 
        date_potent_set[2][is_it_bnAbs_RMs_set[2]]]
    s_bnAbs_RMs = [s_potent_set[1][is_it_bnAbs_RMs_set[1]]; 
        s_potent_set[2][is_it_bnAbs_RMs_set[2]]]

    return(p1, p_violin, s_BNB_set, date_BNB_set, s_notBNB_set, date_notBNB_set, date_bnAbs_RMs, s_bnAbs_RMs)
end;
