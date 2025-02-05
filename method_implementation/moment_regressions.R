mom_marginal <- function(X, seq_depth, irls=T, verbose=T){
    if (is.null(seq_depth)) {
        seq_depth = apply(X, 1, sum, na.rm = T)
    }
    if(nrow(X) != length(seq_depth)){
        stop('The length of the sequencing depth must match the number of cells.')
    }
    n_cell = nrow(X)
    n_gene = ncol(X)
    seq_depth_sq = seq_depth^2
    seq_2 = sum(seq_depth_sq)
    seq_4 = sum(seq_depth^4)
    mu = colSums(X * seq_depth)/seq_2
    M = outer(seq_depth, mu)
    X_centered = X - M
    sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4

    if(irls){
        theta = mu^2/sigma2
    j = 0
    delta = Inf
    
    while( delta > 0.05 & j <= 10 ){
        theta_previous = theta
        theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
        theta[theta < 0] = Inf
        w = M + outer(seq_depth_sq, mu^2/theta_median)
        w[is.na(w)|w <= 0] = 1
        mu = colSums((X/w) * seq_depth)/colSums(seq_depth_sq/w)
        M = outer(seq_depth, mu)
        X_centered = X - M
        h = (M^2/theta_median + M)^2
        h[h <= 0] = 1
        sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/colSums(seq_depth_sq^2/h)
        theta = mu^2/sigma2
        j = j+1
        # print(paste0("Iteration: ", j, ", Median: ", theta_median))
        delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
        # print(delta)
    }
    if(j == 10 & delta > 0.05 & verbose){
        print('IRLS failed to converge after 10 iterations. Please check your data.')
    }else if(verbose){
        print(sprintf('IRLS converged after %i iterations.', j))
    }
    theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    theta[theta < 0] = Inf
    w = M + outer(seq_depth_sq, mu^2/theta_median)
    w[is.na(w)|w <= 0] = 1
    }else{
        w <- NULL
    }
    
    return(list(mu = mu, sigma_sq = sigma2, w = w, X_centered = X - outer(seq_depth, mu),
               X_var = outer(seq_depth, mu) + outer(seq_depth^2, sigma2)))
}

# as there are only one gene/peak in the regression
# theta_median will be NA if sigma_sq is estimated to be negative
# in this case, maybe we can just use s\mu+s^2\sigma^2 as the weight
mom_marginal_one_pair <- function (X, seq_depth, irls = T, verbose = T) 
{
    if (is.null(seq_depth)) {
        seq_depth = apply(X, 1, sum, na.rm = T)
    }
    if (nrow(X) != length(seq_depth)) {
        stop("The length of the sequencing depth must match the number of cells.")
    }
    n_cell = nrow(X)
    n_gene = ncol(X)
    seq_depth_sq = seq_depth^2
    seq_2 = sum(seq_depth_sq)
    seq_4 = sum(seq_depth^4)
    mu = colSums(X * seq_depth)/seq_2
    M = outer(seq_depth, mu)
    X_centered = X - M
    sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4
    if (irls) {
        j = 0
        delta = Inf
        while (delta > 0.05 & j <= 10) {
            mu_prev <- mu
            w <- M + outer(seq_depth_sq, sigma2)
            mu = colSums((X/w) * seq_depth)/colSums(seq_depth_sq/w)
            M = outer(seq_depth, mu)
            X_centered = X - M
            h = (M + outer(seq_depth_sq, sigma2))^2
            sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/colSums(seq_depth_sq^2/h)
            sigma2[sigma2 < 0] <- 0
            # theta = mu^2/sigma2
            j = j + 1
            delta = max(abs(log(mu/mu_prev)))
        }
        if (j == 10 & delta > 0.05 & verbose) {
            print("IRLS failed to converge after 10 iterations. Please check your data.")
        }
        else if (verbose) {
            print(sprintf("IRLS converged after %i iterations.", 
                j))
        }
        w = M + outer(seq_depth_sq, sigma2)
    }
    else {
        w <- NULL
    }
    return(list(mu = mu, sigma_sq = sigma2, w = w, X_centered = X - outer(seq_depth, mu),
               X_var = outer(seq_depth, mu) + outer(seq_depth_sq, sigma2)))
}

mom_covariance <- function(x_gene, X_peak, 
                           w_gene, W_peak, 
                           gene_seq_depth, peak_seq_depth,
                           var_gene, var_peak){
    s_gp <- gene_seq_depth * peak_seq_depth
    g <- w_gene * W_peak
    nume <- colSums(x_gene * X_peak / g * s_gp)
    deno <- colSums(var_gene * var_peak / g^2 * s_gp^2)
    est <- nume / colSums(s_gp^2 / g)
    ts <- nume / sqrt(deno)
    return(list(covar = est, test_stat = ts))
}

mom_marginal_adj_for_batch <- function(X, seq_depth, batch, irls=T, verbose=T){
    if (is.null(seq_depth)) {
        seq_depth = apply(X, 1, sum, na.rm = T)
    }
    if(nrow(X) != length(seq_depth)){
        stop('The length of the sequencing depth must match the number of cells.')
    }
    n_cell = nrow(X)
    n_gene = ncol(X)
    seq_depth_sq = seq_depth^2
    seq_2 = sum(seq_depth_sq)
    seq_4 = sum(seq_depth^4)

    # batch variable
    batch_logical <- matrix(as.logical(batch), nrow = nrow(batch), ncol = ncol(batch))
    n_batches <- ncol(batch)
    # mu = colSums(X * seq_depth)/seq_2
    # M = outer(seq_depth, mu)
    M <- matrix(nrow = n_cell, ncol = n_gene)
    mu_vec <- matrix(NA, nrow = n_gene, ncol = n_batches)
    tmp <- X * seq_depth
    for(i_b in 1:n_batches){
        b_inds <- batch_logical[,i_b]
        mu_vec[, i_b] <- colSums(tmp[b_inds,,drop=F])/sum(seq_depth_sq[b_inds])
        M[b_inds, ] <- outer(seq_depth[b_inds], mu_vec[, i_b])
    }
    X_centered = X - M
    sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4

    if(irls){
        # theta = mu^2/sigma2
        theta = mu_vec^2/sigma2
    j = 0
    delta = Inf
    
    while( delta > 0.05 & j <= 10 ){
        theta_previous = theta
        #theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
        ## new
        theta_median <- numeric(n_cell)
        for(i_b in 1:n_batches){
            # for all samples in the same batch
            # obtain the median of mu^2/sigma^2 as a common median
            theta_median[batch_logical[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
        }
        theta[theta < 0] = Inf
        # no batch effect:
        # # w = M + outer(seq_depth_sq, mu^2/theta_median)
        # with adjustment for batch effect:
        w = M + M^2/theta_median
        w[is.na(w)|w <= 0] = 1
        # no batch effect:
        # mu = colSums((X/w) * seq_depth)/colSums(seq_depth_sq/w)
        # M = outer(seq_depth, mu)
        # with adjustment for batch effect:
        mu_vec <- matrix(NA, nrow = n_gene, ncol = n_batches)
        num_tmp <- (X/w) * seq_depth
        deno_tmp <- seq_depth_sq/w
        for(i_b in 1:n_batches){
            b_inds <- batch_logical[,i_b]
            mu_vec[, i_b] <- colSums(num_tmp[b_inds,,drop=F])/colSums(deno_tmp[b_inds,,drop=F])
            M[b_inds, ] <- outer(seq_depth[b_inds], mu_vec[, i_b])
        }
        X_centered = X - M
        h = (M^2/theta_median + M)^2
        h[h <= 0] = 1
        sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/colSums(seq_depth_sq^2/h)
        # theta = mu^2/sigma2
        theta = mu_vec^2/sigma2
        j = j+1
        # print(paste0("Iteration: ", j, ", Median: ", theta_median))
        delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
        # print(delta)
    }
    if(j == 10 & delta > 0.05 & verbose){
        print('IRLS failed to converge after 10 iterations. Please check your data.')
    }else if(verbose){
        print(sprintf('IRLS converged after %i iterations.', j))
    }
    #theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    ## new
    theta_median <- numeric(n_cell)
    for(i_b in 1:n_batches){
        # for all samples in the same batch
        # obtain the median of mu^2/sigma^2 as a common median
        theta_median[batch_logical[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
    }
        
    theta[theta < 0] = Inf
    #w = M + outer(seq_depth_sq, mu^2/theta_median)
    w = M + M^2/theta_median
    w[is.na(w)|w <= 0] = 1
    }else{
        w <- NULL
    }
    
    return(list(mu = mu_vec, sigma_sq = sigma2, w = w, X_centered = X_centered,
               X_var = M + outer(seq_depth^2, sigma2), M = M))
}



mom_marginal_adj_for_batch_constant_od <- function(X, seq_depth, batch, irls=T, verbose=T){
    if (is.null(seq_depth)) {
        seq_depth = apply(X, 1, sum, na.rm = T)
    }
    if(nrow(X) != length(seq_depth)){
        stop('The length of the sequencing depth must match the number of cells.')
    }
    n_cell = nrow(X)
    n_gene = ncol(X)
    seq_depth_sq = seq_depth^2
    seq_2 = sum(seq_depth_sq)
    seq_4 = sum(seq_depth^4)

    # batch variable
    batch_logical <- matrix(as.logical(batch), nrow = nrow(batch), ncol = ncol(batch))
    n_batches <- ncol(batch)
    # mu = colSums(X * seq_depth)/seq_2
    # M = outer(seq_depth, mu)
    M <- matrix(nrow = n_cell, ncol = n_gene)
    mu_vec <- matrix(NA, nrow = n_gene, ncol = n_batches)
    tmp <- X * seq_depth
    for(i_b in 1:n_batches){
        b_inds <- batch_logical[,i_b]
        mu_vec[, i_b] <- colSums(tmp[b_inds,,drop=F])/sum(seq_depth_sq[b_inds])
        M[b_inds, ] <- outer(seq_depth[b_inds], mu_vec[, i_b])
    }
    mu_marginal <- colSums(tmp) / sum(seq_depth_sq)
    X_centered = X - M
    sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4

    if(irls){
        # theta = mu^2/sigma2
        # replace 'new' with constant od
        # theta = mu_vec^2/sigma2
        theta <- mu_marginal^2 / sigma2
    j = 0
    delta = Inf
    
    while( delta > 0.05 & j <= 10 ){
        theta_previous = theta
        #theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
        ## new
        #theta_median <- numeric(n_cell)
        #for(i_b in 1:n_batches){
        #    # for all samples in the same batch
        #    # obtain the median of mu^2/sigma^2 as a common median
        #    theta_median[batch_logical[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
        #}
        ## replace 'new' with constant od parameter
        theta_median <- stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
        theta[theta < 0] = Inf
        # no batch effect:
        # # w = M + outer(seq_depth_sq, mu^2/theta_median)
        # with adjustment for batch effect:
        w = M + M^2/theta_median
        w[is.na(w)|w <= 0] = 1
        # no batch effect:
        # mu = colSums((X/w) * seq_depth)/colSums(seq_depth_sq/w)
        # M = outer(seq_depth, mu)
        # with adjustment for batch effect:
        mu_vec <- matrix(NA, nrow = n_gene, ncol = n_batches)
        num_tmp <- (X/w) * seq_depth
        deno_tmp <- seq_depth_sq/w
        for(i_b in 1:n_batches){
            b_inds <- batch_logical[,i_b]
            mu_vec[, i_b] <- colSums(num_tmp[b_inds,,drop=F])/colSums(deno_tmp[b_inds,,drop=F])
            M[b_inds, ] <- outer(seq_depth[b_inds], mu_vec[, i_b])
        }
        mu_marginal <- colSums(num_tmp) / colSums(deno_tmp)
        X_centered = X - M
        h = (M^2/theta_median + M)^2
        h[h <= 0] = 1
        sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/colSums(seq_depth_sq^2/h)
        # theta = mu^2/sigma2
        ## new
        # theta = mu_vec^2/sigma2
        ## replace 'new' with constant od
        theta <- mu_marginal^2 / sigma2
        j = j+1
        # print(paste0("Iteration: ", j, ", Median: ", theta_median))
        delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
        # print(delta)
    }
    if(j == 10 & delta > 0.05 & verbose){
        print('IRLS failed to converge after 10 iterations. Please check your data.')
    }else if(verbose){
        print(sprintf('IRLS converged after %i iterations.', j))
    }
    #theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    ## new
    # theta_median <- numeric(n_cell)
    # for(i_b in 1:n_batches){
    #     # for all samples in the same batch
    #     # obtain the median of mu^2/sigma^2 as a common median
    #     theta_median[batch_logical[,i_b]] <- stats::quantile(theta[, i_b][theta[, i_b] > 0], na.rm = T, probs = 0.5)
    # }
    # replace 'new' with constant od
    theta[theta < 0] = Inf
    theta_median <- stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    #w = M + outer(seq_depth_sq, mu^2/theta_median)
    w = M + M^2/theta_median
    w[is.na(w)|w <= 0] = 1
    }else{
        w <- NULL
    }
    
    return(list(mu = mu_vec, sigma_sq = sigma2, w = w, X_centered = X_centered,
               X_var = M + outer(seq_depth^2, sigma2), M = M))
}

mom_covariance_adj_for_batch <- function(x_gene, X_peak, 
                                         w_gene, W_peak, 
                                   gene_seq_depth, peak_seq_depth,
                                   var_gene, var_peak,
                                        batch_mat){
    s_gp <- gene_seq_depth * peak_seq_depth
    U <- cbind(s_gp, batch_mat)
    g <- w_gene * W_peak
    y <- x_gene * X_peak
    nume <- solve(t(U/g) %*% U) %*% t(U/g) %*% y
    #colSums(x_gene * X_peak / g * s_gp)
    #deno <- colSums(var_gene * var_peak / g^2 * s_gp^2)
    deno <- solve(t(U/g) %*% U) %*% t(U/g)
    est <- nume / colSums(1 / g * s_gp^2)
    ts <- nume / sqrt(deno)
    return(list(covar = est, test_stat = ts))
}