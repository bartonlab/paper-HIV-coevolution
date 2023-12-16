function split_string(s::String)
    this_match = match(r"\d+", s)  # Finds the number in the string
    if this_match == nothing
        error("No number found in string.")
    end

    number_position = this_match.offset  # Position where the number starts
    number = parse(Int, this_match.match)  # The number itself
    left_part = split(s[1:number_position-1], "/")
    right_part = split(s[number_position + length(this_match.match):end], "/")

    # Returns three parts: left, number, right
    return (unique(left_part), number, unique(right_part))
end

#using Plots, StatsPlots, KernelDensity
x = randn(1000)
y = randn(1000)
#title_scatter, title_x, title_y = "", "", ""
function get_scatter_density_plot(x, y, title_scatter, title_x, title_y)
    # Kernel Density Estimations
    x_kde = kde(x)
    y_kde = kde(y)
    #x_kde.density = x_kde.density[x_kde.density .> 0]
    #y_kde.density = y_kde.density[y_kde.density .> 0]
    # Creating the plot layout
    l = @layout [a{0.8w} _; c{0.8h} d]
    # Plotting
    p = plot(x, y, st=:scatter, legend=false, title=title_scatter, subplot=2, layout=l, markerstrokewidth=0, markeralpha=0.4)
    plot!(p[1], x_kde.x, x_kde.density, legend=false, title=title_x, subplot=1, layout=l) #, yscale=:log10)
    plot!(p[3], y_kde.density, y_kde.x, legend=false, title=title_y, subplot=3, layout=l) #, xscale=:log10)
    return p
end;

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

function get_collection_time(csv_raw_in)
    name_time_collection = names(csv_raw_in)
    colection_time = []
    for x in split.(name_time_collection, "f_at_")
        if(length(x)>1)
            push!(colection_time, parse(Int, x[2]) )
        end;
    end
    return colection_time 
end;

# Note csv_raw_in is summarized yet containing trajecory, but csv_raw_full_in is 
function get_poly_site(q, L, csv_raw_in, csv_raw_full_in, th_below=0.05, th_above=0.95)
    colection_time = get_collection_time(csv_raw_in)
    qL = q*L
    
    n_time_max = length(colection_time)
    x1_in = zeros(qL, length(colection_time))
    for i_t in 1:n_time_max
        x1_in[:, i_t] = csv_raw_full_in[1:qL, Symbol("f_at_" * string( colection_time[i_t]))]
    end;

    poly_in = []
    for i in 1:qL
        n_below = count(x1_in[i, :] .<= th_below)
        n_above = count(x1_in[i, :] .>= th_above)
        if(n_below!=n_time_max && n_above!=n_time_max)
            push!(poly_in, true)
        end
        if(n_below==n_time_max || n_above==n_time_max)
            push!(poly_in, false)
        end
    end;
    return (poly_in, colection_time)
end;
linreg2(x,y) = hcat(fill!(similar(x),1),x) \ y;

get_selection(num, cov, γ) = (cov + γ*I) \ num;

function get_cov_num(q, L, fname)
    csv_raw = readdlm(fname);
    cov = zeros(q*L, q*L)
    num = zeros(q*L);
    dx = zeros(q*L);
    for n in 1:size(csv_raw,1)
        x = csv_raw[n, :]
        if(x[1]=="C")
            i,j,v = x[2], x[3], x[4]
            cov[i,j] = v
            cov[j,i] = v
        end
        if(x[1]=="num")
            i,v = x[2],x[3]
            num[i]=v
        end
        if(x[1]=="dx")
            i,v = x[2],x[3]
            dx[i]=v
        end
    end;
    csv_raw = [];
    return (cov, num, dx)
end;
function get_color(qL, poly_common, index_likely_positive, index_likely_negative)
    plot_color = []
    for i in 1:qL
        if(poly_common[i])
            if(index_likely_positive[i])
                push!(plot_color, "red")
            elseif(index_likely_negative[i])
                push!(plot_color, "blue")
            else
                push!(plot_color, "gray")
            end
        end
    end
    return plot_color
end;

n = 5
function get_traject_plot(csv_raw_in, colection_time, x1_selected, n)
    #this_title = @sprintf("AA: %s\nDNA: %s", csv_raw_in.mutants_AA[n] , csv_raw_in.mutants_nuc[n] )
    this_title = @sprintf("%s", csv_raw_in.mutants_AA_filtered_fr3[n] )
    this_label= @sprintf("s = %.2f (%%)", 100*csv_raw_in.s_MPL[n])
    p1 = Plots.plot(colection_time, 100*x1_selected[n, :], title = this_title, m=:o, 
        label=this_label, foreground_color_legend = nothing, 
        markerstrokewidth=0, ms=3, ylim=(-5, 105), 
#        xlabel="Days", ylabel="Frequency (%)", margin=2mm
    )
    return p1
end;
# Create a 3x4 grid of subplots
function get_traject_plot(csv_raw_in, colection_time, x1_selected, nrows, ncols, size_x, size_y)
    # Initialize the plot
    p = plot(layout = (nrows, ncols), size = (size_x, size_y), margin=5mm)

    # Data (example data - replace with your own)
    for i in 1:nrows
        for j in 1:ncols
            plot_idx = (i - 1) * ncols + j
            this_title = @sprintf("%s", csv_raw_in.mutants_AA_filtered_fr3[plot_idx] )
            this_label= @sprintf("s = %.2f (%%)", 100*csv_raw_in.s_MPL[plot_idx])
            # Plot the data
            plot!(p[plot_idx], 
                colection_time, 100*x1_selected[plot_idx, :], title = this_title, m=:o, 
                label=this_label, foreground_color_legend = nothing, 
                markerstrokewidth=0, ms=3, ylim=(-5, 105))

            # Add x-axis label for plots on the bottom row
            if i == nrows
                xlabel!(p[plot_idx], "Days")
            end

            # Add y-axis label for plots in the leftmost column
            if j == 1
                ylabel!(p[plot_idx], "Frequency")
            end
        end
    end

    # Show the plot
    display(p)
end;

function get_freq(csv_raw_full, row_eff)
    headers = names(csv_raw_full)
    headers_collet_time = []
    for x in headers
        if( split(x, "at_")[1]=="f_")
           push!(headers_collet_time, x) 
        end
    end
    freq_mut = []
    for x in headers_collet_time
        push!(freq_mut, csv_raw_full[:, x][row_eff])
    end
    return freq_mut
end;

function get_lin_integration(x,y,k_max)
    Δx = [i>1 ? x[i]-x[i-1] : x[i] for i in 1:length(x)]    
    return lin_int = sum(y[1:k_max] .* Δx[1:k_max])
end


# Function to extract integer from a string
function extract_integer(str)
    # Use a regular expression to find all sequences of digits
    matches = match(r"\d+", str)

    # If there are matches, convert the first match to an integer
    if matches !== nothing
        return parse(Int, matches.match)
    else
        return nothing # or some default value/error handling
    end
end;

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h]);

# Custom transformation function
function transform_scale(x; a = 1) return log10(a + x) end
#[transform_scale(x, a=1) for x in rand(10)]

function get_variable_regions(i)
    if(i ∈ 131:157) return "V1" end
    if(i ∈ 158:196) return "V2" end
    if(i ∈ 275:283) return "LD" end
    if(i ∈ 296:331) return "V3" end
    if(i ∈ 385:418) return "V4" end
    if(i ∈ 460:470) return "V5" end
end;






