# simulate both gene counts and peak counts based on a joint bivariate expression-measurement model
sim_exper <- function(rho, n, s_gene, s_peak, mu_vec, var_vec, n_reps, seed=2023, copula = T){
    set.seed(seed)
    x_peak_mat <- matrix(0,n_reps,n)
    x_gene_mat <- matrix(0,n_reps,n)
    theta_vec <- mu_vec^2 / var_vec
    sigma_vec <- var_vec/mu_vec

    for(i in 1:n_reps){
        if(copula){
            # gaussian copula
            ## Note: z should be simulated for each iteration
            gp_mat <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1), nrow=2))
            z_gene <- qgamma(pnorm(gp_mat[,1]), shape=theta_vec[1], scale=sigma_vec[1])
            z_peak <- qgamma(pnorm(gp_mat[,2]), shape=theta_vec[2], scale=sigma_vec[2])   
        }else{
            # linear combination of Gamma's
            # 5.1 (ii) from https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1442340
            Sigma <- matrix(c(1,rho,rho,1),2,2)
            mat_de <- svd(Sigma) 
            A <- (sqrt(mat_de$d) * mat_de$u) %>% t # A %*% t(A) == Sigma
            z_01 <- rgamma(n, shape=10, 1)
            z_02 <- rgamma(n, shape=10, 1)
            # mean(z) = mu, var(z) = Sigma
            # however, simulated z's are not guaranteed to be positive
            z <- A %*% matrix(c(z_01, z_02), nrow = 2, byrow = T) + matrix(rep(mu_vec - rowSums(A),n), nrow=2, byrow=F) 
            z_gene <- z[1,]
            z_peak <- z[2,]
        }
        x_peak <- rpois(n, z_peak * s_peak)
        x_gene <- rpois(n, z_gene * s_gene)
        x_peak_mat[i,] <- x_peak
        x_gene_mat[i,] <- x_gene
    }

    # save a gene matrix of n_reps * n, f n_reps * n
    return(list(x_gene_mat=x_gene_mat,x_peak_mat=x_peak_mat)) # , copula_cor=cor(z_peak, z_gene)
}

sim_exper_with_batch_effects <- function(rho, n, s_gene, s_peak, mu_mat, var_vec, n_reps, batch_inds_list, seed=2023,
                                        copula = T){
    set.seed(seed)
    x_peak_mat <- matrix(0,n_reps,n)
    x_gene_mat <- matrix(0,n_reps,n)
    batches <- names(batch_inds_list)

    z_gene_full <- z_peak_full <- matrix(0,n_reps,n)
    for(batch in batches){
        theta_vec <- mu_mat[,batch]^2 / var_vec
        sigma_vec <- var_vec / mu_mat[,batch]
        n_batch <- length(batch_inds_list[[batch]])
        binds <- batch_inds_list[[batch]]

        if(!copula){
            mu_vec <- mu_mat[,batch]
                covar <- rho * sqrt(prod(var_vec))
                a_2 <- (mu_vec[1] - covar)^2 / (var_vec[1] - covar)
                a_3 <- (mu_vec[2] - covar)^2 / (var_vec[2] - covar)
                b_1 <- (var_vec[1] - covar) / (mu_vec[1] - covar)
                c_2 <- (var_vec[2] - covar) / (mu_vec[2] - covar)
                a_1 <- covar
            
            print(sprintf('%.4f, %.4f, %.4f', a_1, a_2, a_3))
            print(sprintf('%.4f, %.4f', b_1, c_2))
        }
        for(i in 1:n_reps){
            # gaussian copula
            if(copula){
                ## Note: z should be simulated for each iteration
                gp_mat <- mvrnorm(n_batch, c(0,0), matrix(c(1,rho,rho,1), nrow=2))
                z_gene_full[i,binds] <- qgamma(pnorm(gp_mat[,1]), shape=theta_vec[1], scale=sigma_vec[1])
                z_peak_full[i,binds] <- qgamma(pnorm(gp_mat[,2]), shape=theta_vec[2], scale=sigma_vec[2])   
            }else{
                
                alpha_0 <- rgamma(n_batch, shape = a_1, scale = 1)
                alpha_1 <- rgamma(n_batch,shape = a_2, scale = 1)
                alpha_2 <- rgamma(n_batch,shape = a_3, scale = 1)
                z_gene_full[i,binds] <- alpha_0 + b_1 * alpha_1
                z_peak_full[i,binds] <- alpha_0 + c_2 * alpha_2
            }
        }
    }
    for(i in 1:n_reps){
        x_gene_mat[i, ] <- rpois(n, z_gene_full[i, ] * s_gene)
        x_peak_mat[i, ] <- rpois(n, z_peak_full[i, ] * s_peak)
    }

    # save a gene matrix of n_reps * n, f n_reps * n
    return(list(x_gene_mat=x_gene_mat,x_peak_mat=x_peak_mat, z_gene_mat=z_gene_full, z_peak_mat=z_peak_full)) # , copula_cor=cor(z_peak, z_gene)
}

sim_exper_with_batch_effects_new <- function(rho, n, s_gene, s_peak, mu_mat, var_vec, n_reps, batch_inds_list, seed=2023,
                                        copula = T){
    set.seed(seed)
    x_peak_mat <- matrix(0,n_reps,n)
    x_gene_mat <- matrix(0,n_reps,n)
    batches <- names(batch_inds_list)

    z_gene_full <- z_peak_full <- matrix(0,n_reps,n)
    i_min_batch_gene <- which.min(mu_mat[1,])
    i_min_batch_peak <- which.min(mu_mat[2,])
    mu_vec <- c(mu_mat[1,i_min_batch_gene],mu_mat[2,i_min_batch_peak])
    theta_vec <- mu_vec^2 / var_vec
    sigma_vec <- var_vec / mu_vec
    
    gamma_cor <- numeric(n_reps)
    for(i in 1:n_reps){
        gp_mat <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1), nrow=2))
        z_gene_full[i,] <- qgamma(pnorm(gp_mat[,1]), shape=theta_vec[1], scale=sigma_vec[1])
        z_peak_full[i,] <- qgamma(pnorm(gp_mat[,2]), shape=theta_vec[2], scale=sigma_vec[2])   
        gamma_cor[i] <- cor(z_gene_full[i,], z_peak_full[i,])
    }
    z_gene_full_raw <- z_gene_full
    z_peak_full_raw <- z_peak_full
    print(summary(gamma_cor))

    for(batch in batches){
        binds <- batch_inds_list[[batch]]
        z_gene_full[,binds] <- z_gene_full[,binds] + mu_mat[1,batch] - mu_vec[1]
        z_peak_full[,binds] <- z_peak_full[,binds] + mu_mat[2,batch] - mu_vec[2]
    }
    print(mu_mat)
    print(mu_mat[1,] - mu_mat[1,i_min_batch_gene])
    print(mu_mat[2,] - mu_mat[2,i_min_batch_peak])
    
    for(i in 1:n_reps){
        x_gene_mat[i, ] <- rpois(n, z_gene_full[i, ] * s_gene)
        x_peak_mat[i, ] <- rpois(n, z_peak_full[i, ] * s_peak)
    }

    # save a gene matrix of n_reps * n, f n_reps * n
    return(list(x_gene_mat=x_gene_mat,x_peak_mat=x_peak_mat, z_gene_mat=z_gene_full, z_peak_mat=z_peak_full,
               z_gene_full_raw=z_gene_full_raw, z_peak_full_raw=z_peak_full_raw)) # , copula_cor=cor(z_peak, z_gene)
}


load_simulated_data <- function(mode, method, i_permu, permu_peak=T, sparse_factor=1){
    if(mode %in% c('permutation', 'permutation_batch_effect')){
        # load permuted rna data
        rna_counts <- sprintf('%s_RNA_%i_pairs_%s_sampled_i_permu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table %>% t
        peak_counts <- ct_obj[[pn]]@counts[sampled_pairs$peak,]
        if(permu_peak){
            rna_counts <- ct_obj[['RNA']]@counts[sampled_pairs$gene,]
            peak_counts <- sprintf('%s_peak_%i_pairs_%s_sampled_i_permu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table %>% t
            print('load permuted peak data')
            print(sprintf('%s_peak_%i_pairs_%s_sampled_i_permu_%i.txt', output_dir, n_pairs, sampling, i_permu) )
        }
    }else if(mode %in% c('simulation', 'simulation_null', 'simulation_batch_effect', 'simulation_null_batch_effect')){
        # load simulated data
        rna_counts <- sprintf('%s_RNA_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table %>% t
        sparse_factor_suffix <- ifelse(sparse_factor == 1, '', sprintf('s*%.1f_', sparse_factor))
	if(sparse_factor != 1) print('load simulated data with modified seq depth')
        peak_counts <- sprintf('%s_peak_%i_pairs_%s_sampled_%si_simu_%i.txt', output_dir, n_pairs, sampling, sparse_factor_suffix, i_permu) %>% read.table %>% t
    }else if(mode %in% c('simulation_null', 'simulation_null_batch_effect')){
        # load simulated data
        rna_counts <- sprintf('%s_RNA_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table %>% t
        # peak_counts <- sprintf('%s_peak_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table %>% t
        peak_counts <- ct_obj[[pn]]@counts[sampled_pairs$peak,]
    }else if(mode %in% c('permutation_peak_counts', 'permutation_peak', 'permutation_counts_batch_effect')){
        # extract counts
        rna_counts <- ct_obj[['RNA']]@counts[sampled_pairs$gene,]
        if(mode %in% c('permutation_peak_counts', 'permutation_counts_batch_effect')){
            fn_str <- 'peak_counts'
        }else{
            fn_str <- 'peak'
        }
        fn <- sprintf('%s_%s_%i_pairs_%s_sampled_i_permu_%i.txt', 
                               permu_output_dir, 
                               fn_str,
                               n_pairs, sampling, i_permu)
        if(i_permu == 1) print(fn)
        peak_counts <- fn %>% read.table %>% t
    }
    if(method %in% c('SCENT', 'Signac')){
        if(!class(rna_counts)[1] %in% c('dtCMatrix', 'dgCMatrix', 'CsparseMatrix')){
            rna_counts <- Matrix(rna_counts, sparse = TRUE) 
        }
        if(!class(peak_counts)[1] %in% c('dtCMatrix', 'dgCMatrix', 'CsparseMatrix')){
            peak_counts <- Matrix(peak_counts, sparse = TRUE)
        }
        rownames(rna_counts) <- sampled_pairs$gene
        rownames(peak_counts) <- sampled_pairs$peak
    }else if(method == 'proposed'){
        rna_counts <- rna_counts %>% as.matrix %>% t
        peak_counts <- peak_counts %>% as.matrix %>% t
    }
    return(list(gene = rna_counts, peak = peak_counts))
}

#    else if(mode == 'simulation_null' | mode == 'simulation_null_permutation'){
        # use gene counts independently simulated from peak counts
        # rna_counts <- sprintf('%s_RNA_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table
        # rna_counts <- Matrix(t(rna_counts), sparse = TRUE) 
        # rownames(rna_counts) <- sampled_pairs$gene
        # peak_counts <- ct_obj[[pn]]@counts

        # use peak counts independently simulated from gene counts
        # peak_counts <- sprintf('%s_peak_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% read.table
        # peak_counts <- Matrix(t(peak_counts), sparse = TRUE)
#    	print(sprintf('%s: use peak counts independently simulated from gene counts', mode))
#    	peak_counts <- sprintf('%s_peak_%i_pairs_%s_sampled_i_simu_%i.txt', output_dir, n_pairs, sampling, i_permu) %>% readMM
#        peak_counts <- as(peak_counts, 'CsparseMatrix')
#        rna_counts <- ct_obj[['RNA']]@counts[sampled_pairs$gene,]
#    }
