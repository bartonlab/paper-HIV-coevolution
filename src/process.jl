get_time_list(data) = Int.(unique(data[:,1]));

km(i,a,q) = (i-1)*q+a

function get_data_time(file_key1, file_key2, data_dir, id_ensemble, time_upto=300)
    fname_in = data_dir * file_key1 *"_id-"*string(id_ensemble) * file_key2
    data = readdlm(fname_in);
    read_upto = count(data[:,1] .<= time_upto)
    data = copy(data[1:read_upto, :])
    time_list = get_time_list(data);
    return (data, time_list)
end;

function get_sample_at_t(data, t_get)
    n_t = []
    sample_t = []
    for n in 1:size(data,1)
        if(Int(data[n, 1]) == t_get)
            #n_t = push!(n_t, data[n,2]) 
            #sample_t = push!(sample_t, data[n,4:end]) 
            if(length(sample_t)>0)
                sample_t = hcat(sample_t, data[n,4:end]) 
                n_t = vcat(n_t, data[n,2])
            end
            if(length(sample_t)==0)
                sample_t = copy(data[n,4:end])
                n_t = data[n,2]
            end
        end

        if(data[n, 1] > t_get)
            break
        end
    end
    return (n_t, sample_t)
end;


# mappinf function that converts (i,j,a,b) indices of mattrix to a scalar index 
# note need to be j>i
stepfunc(x) = x>0 ? 1 : 0
G(i,j,a,b,q,L) = stepfunc(j-i)*Int(q^2*(i-1)*(L-i/2.0) + (a-1)*q*(L-i) + (j-i-1)*q + b);

# Note that suppose seq_in is not binary seqeune
function get_x_1w2_cat(seq_in, q, L, qL, qqLLhalf)
    x1 = zeros(qL)
    x2_vec = zeros(qqLLhalf);
    for i in 1:L
        a = seq_in[i]
        u = km(i,a,q)
        x1[u] = 1
        for j in (i+1):L
            b = seq_in[j]
            ξ = G(i,j,a,b,q,L)
            x2_vec[ξ] = 1
        end
    end
    return ([x1;x2_vec])
end;

function get_d_mu(x_1w2, q, L, qL, qqLLhalf)
    d_mu_temp = zeros(qL+qqLLhalf)
    for i in 1:L
        for a in 1:q
	    u = km(i,a,q)
	    for c in 1:q 
	     	if(c != a)	
			u_c = km(i,c,q)
			d_mu_temp[u] += x_1w2[u_c] - x_1w2[u]
		end
            end 
	    
	    for j in (i+1):L
                for b in 1:q
                    ξ = qL + G(i,j,a,b,q,L)

		    for c in 1:q
                    	ξ_ac = qL + G(i,j,a,c,q,L)
                    	ξ_cb = qL + G(i,j,c,b,q,L)
			if(c != a)
				d_mu_temp[ξ] += x_1w2[ξ_ac] - x_1w2[ξ]
			end
			if(c != b)
				d_mu_temp[ξ] += x_1w2[ξ_cb] - x_1w2[ξ]
			end
		    end 
                end
            end
        end
    end
    return d_mu_temp
end;


function get_d_mu(x_1w2, q, L, qL, qqLLhalf, muMat)
    d_mu_temp = zeros(qL+qqLLhalf)
    for i in 1:L
        for a in 1:q
	    u = km(i,a,q)
	    for c in 1:q 
	     	if(c != a)	
			u_c = km(i,c,q)
			d_mu_temp[u] += muMat[c,a] * x_1w2[u_c] - muMat[a,c] * x_1w2[u]
		end
            end 
	    
	    for j in (i+1):L
                for b in 1:q
                    ξ = qL + G(i,j,a,b,q,L)

		    for c in 1:q
                    	ξ_ac = qL + G(i,j,a,c,q,L)
                    	ξ_cb = qL + G(i,j,c,b,q,L)
			if(c != a)
				d_mu_temp[ξ] += muMat[c,a] * x_1w2[ξ_ac] - muMat[a,c] * x_1w2[ξ]
			end
			if(c != b)
				d_mu_temp[ξ] += muMat[c,b] * x_1w2[ξ_cb] - muMat[b,c] * x_1w2[ξ]
			end
		    end 
                end
            end
        end
    end
    return d_mu_temp
end;

# The indices is tricky, aligned site indices (or polygenetic sites) are not exactly the sequence indices.
function get_d_rec(x_1w2, q, L, qL, qqLLhalf)
    d_rec_temp = zeros(qL+qqLLhalf)
    for i in 1:L
        for a in 1:q
            u = km(i,a,q)
            for j in (i+1):L
                for b in 1:q
            	    v = km(j,b,q)
                    ξ = qL + G(i,j,a,b,q,L)
                    #this should be modified for categorical variables case. 
                    d_rec_temp[ξ] = (j-i) * (x_1w2[ξ] - x_1w2[u] * x_1w2[v]) 
                end
            end
        end
    end
    return d_rec_temp
end;

# The indices is tricky, aligned site indices (or polygenetic sites) are not exactly the sequence indices.
function get_d_rec(x_1w2, q, L, qL, qqLLhalf, indices_mapped)
    d_rec_temp = zeros(qL+qqLLhalf)
    for i in 1:L
        for a in 1:q
            u = km(i,a,q)
            for j in (i+1):L
                for b in 1:q
            	    v = km(j,b,q)
                    ξ = qL + G(i,j,a,b,q,L)
                    #this should be modified for categorical variables case. 
                    i_eff, j_eff = indices_mapped[i], indices_mapped[j]
		    #@printf("%d, %d, %d, %.3f, %.3f\n", (j_eff - i_eff), j_eff, i_eff, x_1w2[ξ] - x_1w2[u] * x_1w2[v], x_1w2[ξ]) 
		    d_rec_temp[ξ] = (j_eff - i_eff) * (x_1w2[ξ] - x_1w2[u] * x_1w2[v]) 
                end
            end
        end
    end
    return d_rec_temp
end

# Note that suppose seq_cat_t is not binary seqeune
function get_4th_Cov_mut_random_compression(seq_cat_t, n_t, q, L, qL, qqLLhalf, rank_x, rank_S, S)
    Neff = Int(sum(n_t))
    #LLhalf = Int(q^2*L*(L-1)/2.0)
    n_species = size(seq_cat_t,1);
    
    x_1w2 = zeros(qL+qqLLhalf) # it is not necssary to have a huge matrix x_3w4.
    for m in 1:n_species
        n_t_scale = n_t[m] 
        
        x_1w2_m = get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf)
        x_1w2 += n_t_scale * x_1w2_m
    end
    x_1w2 /= Neff
    d_mu = get_d_mu(x_1w2, q, L, qL, qqLLhalf);

    # check also identical S case. 
    CS_tot = zeros(rank_S, rank_S)
    #C_tot = zeros(rank_x, rank_x)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff

        if(m>n_species)
            # centerize and normalize
            x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf) - x_1w2 )
            # random matrix compression 
            S_x_1w2_m = S' * x_1w2_m
        end

        # centerize and normalize
        x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf) - x_1w2 )
        # random matrix compression 
        S_x_1w2_m = S' * x_1w2_m
        
        # Want to avoid inv function    
        C_S = S_x_1w2_m * S_x_1w2_m'
        #C = x_1w2_m * x_1w2_m'
        CS_tot += C_S
        #C_tot += C
    end
    #matrix can be recoverd by P * CS_tot * P' ~ C_tot;
    return (CS_tot, x_1w2, d_mu)
end;


function get_x_1w2_set_and_d_mu_d_rec_random_compression(seq_cat_t, n_t, q, L, qL, qqLLhalf, rank_x, muMat)
    Neff = Int(sum(n_t))
    n_species = size(seq_cat_t,1);
    x_1w2 = zeros(qL+qqLLhalf) # it is not necssary to have a huge matrix x_3w4.
    for m in 1:n_species
        n_t_scale = n_t[m]     
        x_1w2_m = get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf)
        x_1w2 += n_t_scale * x_1w2_m
    end
    
    x_1w2 /= Neff
    
    d_mu = get_d_mu(x_1w2, q, L, qL, qqLLhalf, muMat);
    d_rec = get_d_rec(x_1w2, q, L, qL, qqLLhalf);
    
    x_1w2_set = zeros(n_species, rank_x)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        # centerize and normalize
        x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf) - x_1w2 )
        x_1w2_set[m, :] = copy(x_1w2_m)
    end
    return (x_1w2_set, x_1w2, d_mu, d_rec)
end;

# The indices_mapped maps compressed sequence to original allignment indices. 
function get_x_1w2_set_and_d_mu_d_rec_random_compression(seq_cat_t, indices_mapped, n_t, q, L, qL, qqLLhalf, rank_x, muMat)
    Neff = Int(sum(n_t))
    n_species = size(seq_cat_t,1);
    x_1w2 = zeros(qL+qqLLhalf) # it is not necssary to have a huge matrix x_3w4.
    for m in 1:n_species
        n_t_scale = n_t[m]     
        x_1w2_m = get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf)
        x_1w2 += n_t_scale * x_1w2_m
    end
    
    x_1w2 /= Neff
    
    d_mu = get_d_mu(x_1w2, q, L, qL, qqLLhalf, muMat);
    d_rec = get_d_rec(x_1w2, q, L, qL, qqLLhalf, indices_mapped);
    
    x_1w2_set = zeros(n_species, rank_x)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        # centerize and normalize
        x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf) - x_1w2 )
        x_1w2_set[m, :] = copy(x_1w2_m)
    end
    return (x_1w2_set, x_1w2, d_mu, d_rec)
end;

function get_x_1w2_set_and_d_mu_random_compression(seq_cat_t, n_t, q, L, qL, qqLLhalf, rank_x)
    Neff = Int(sum(n_t))
    n_species = size(seq_cat_t,1);
    x_1w2 = zeros(qL+qqLLhalf) # it is not necssary to have a huge matrix x_3w4.
    for m in 1:n_species
        n_t_scale = n_t[m]     
        x_1w2_m = get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf)
        x_1w2 += n_t_scale * x_1w2_m
    end
    
    x_1w2 /= Neff
    d_mu = get_d_mu(x_1w2, q, L, qL, qqLLhalf);
    x_1w2_set = zeros(n_species, rank_x)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        # centerize and normalize
        x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2_cat(seq_cat_t[m, 1:L], q, L, qL, qqLLhalf) - x_1w2 )
        x_1w2_set[m, :] = copy(x_1w2_m)
    end

    return (x_1w2_set, x_1w2, d_mu)
end;



function get_4th_Cov_random_compression(S, x_1w2_set)
    S_x1w2_set = S' * x_1w2_set' # (m x d) x (d x n_species) --> (m x n_species) mat
    CS = S_x1w2_set * S_x1w2_set' # \sum_i^n_species (m x m)_i
    return CS
end;


function check_removing_pairs_condition(q,L, x2_mask, x_1w2, cond1=1.0/3, cond2=1.0/4^2)
	for i in 1:L
	for a in 1:q
		x_i = x_1w2[km(i,a,q)]
		for j in (i+1):L
		for b in 1:q
			x_j = x_1w2[km(j,b,q)]
			ξ = G(i,j,a,b,q,L)
			x_ij = x_1w2[qL+ξ]
			#if(x_ij > cond2 || x_i > cond1 || x_j > cond1)
			if(x_ij > cond2 || (x_i > cond1 && x_j > cond1))
				x2_mask[km(i,a,q), km(j,b,q)] = 1
				x2_mask[km(j,b,q), km(i,a,q)] = 1
			end
		end
		end
	end
	end
	return x2_mask
end;
