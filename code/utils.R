get_index_mat <- function(nChar, nStates, nRateClass){
  # create the appropriate sized rate matrix and dataset
  all_states <- list()
  for(i in 1:nChar){
    all_states[[i]] <- 1:nStates
  }
  min_df <- data.frame(sp="NA", expand.grid(all_states))
  tmp <- getStateMat4Dat(min_df)
  rate_mat <- tmp$rate.mat
  legend <- tmp$legend
  if(nRateClass>1){
    rate_cat_mat <- getRateCatMat(nRateClass)
    rate_mat_list <- list()
    for(i in 1:nRateClass){
      rate_mat_list[[i]] <- rate_mat
    }
    full_rate_mat <- getFullMat(rate_mat_list, rate_cat_mat)
  }else{
    full_rate_mat <- rate_mat
  }
  out <- list(legend=legend, 
              nChar=nChar, 
              nStates=nStates,
              nRateClass=nRateClass,
              min_df=min_df, 
              rate_mat=rate_mat, 
              full_rate_mat=full_rate_mat)
  return(out)
}

get_par_table <- function(index_mat, nSim, mean=0, sd=1){
  nPar <- max(index_mat$full_rate_mat)
  par_table <- do.call(rbind, 
                       lapply(1:nSim, function(x) 
                         10^rnorm(nPar, mean=mean, sd=sd)))
  par_names <- c()
  for(i in 1:nPar){
    index <- which(i == index_mat$full_rate_mat, arr.ind = TRUE)
    from <- rownames(index_mat$full_rate_mat)[index[1]]
    to <- colnames(index_mat$full_rate_mat)[index[2]]
    par_names <- c(par_names, paste0(from, "_", to))
  }
  colnames(par_table) <- par_names
  return(par_table)
}

get_rate_mats <- function(index_mat, par_table){
  rate_mats <- list()
  rate_index_mat <- index_mat$full_rate_mat
  for(i in 1:nrow(par_table)){
    pars <- par_table[i,]
    rate_mat <- index_mat$full_rate_mat
    for(j in 1:length(pars)){
      par_index <- which(rate_index_mat == j, arr.ind = TRUE)
      for(k in 1:nrow(par_index)){
        rate_mat[par_index[k,1], par_index[k,2]] <- unlist(pars[j])
      }
    }
    diag(rate_mat) <- -rowSums(rate_mat)
    rate_mats[[i]] <- rate_mat
  }
  return(rate_mats)
}

get_sim_data <- function(phy, rate_mat, index_mat, root.p=NULL){
  if(is.null(root.p)){
    root.p <- c(Null(rate_mat))
  }
  dat <- corHMM:::simMarkov(phy, rate_mat, root.p)
  return(dat)
}

get_formatted_data <- function(dat, index_mat){
  states <- sapply(dat, function(x) colnames(index_mat)[x])
  tmp_df <- do.call(rbind, strsplit(states, "\\|"))
  tmp_state_names <- do.call(rbind, strsplit(colnames(index_mat), "\\|"))
  unique_states <- apply(tmp_state_names, 2, unique)
  cor_dat <- data.frame(sp=names(dat))
  for(i in 1:ncol(unique_states)){
    cor_dat[,i+1] <- factor(tmp_df[,i], unique_states[,i])
  }
  rownames(cor_dat) <- names(dat)
  return(cor_dat)
}

get_solution_from_res <- function(res, index.mat = NULL){
  if(class(res) == "mcmc"){
    tmp <- summary(res)$quantiles
    p <- tmp[,3]
  }
  if(class(res) == "try-error"){
    p <- NA
  }
  if(class(res) == "corhmm"){
    if(!is.null(index.mat)){
      if(dim(index.mat)[1] == dim(res$index.mat)[1]){
        p <- sapply(1:max(res$index.mat, na.rm = TRUE), function(x) 
          na.omit(c(res$solution))[na.omit(c(res$index.mat) == x)][1])
      }else{
        p <- NA
      }
    }else{
      p <- sapply(1:max(res$index.mat, na.rm = TRUE), function(x) 
        na.omit(c(res$solution))[na.omit(c(res$index.mat) == x)][1])
    }
  }
  return(p)
}

get_par_from_rate_mat <- function(dat, index_mat){
  p <- sapply(1:max(index_mat$full_rate_mat, na.rm = TRUE), function(x) 
    na.omit(c(dat$par))[na.omit(c(index_mat$full_rate_mat) == x)][1])
  return(p)
}

get_better_df <- function(df, col_nm, type, ntips){
  df <- as.data.frame(df)
  colnames(df) <- col_nm
  df <- cbind(ntips = ntips, type = type, df)
  df_longer <- df %>%
    pivot_longer(cols = col_nm, names_to = "par", values_to = "value")
  return(df_longer)
}

get_rate_mat <- function(index_mat, pars){
  rate_index_mat <- index_mat$full_rate_mat
  rate_mat <- index_mat$full_rate_mat
  for(j in 1:length(pars)){
    par_index <- which(rate_index_mat == j, arr.ind = TRUE)
    for(k in 1:nrow(par_index)){
      rate_mat[par_index[k,1], par_index[k,2]] <- unlist(pars[j])
    }
  }
  diag(rate_mat) <- -rowSums(rate_mat)
  return(rate_mat)
}

model_test_dep <- function(cor_obj, index_mat){
  Q <- cor_obj$solution
  Q[is.na(Q)] <- 0
  test_1 <- Q[index_mat == 8] > Q[index_mat == 3] # ovi->vivi|cold > |warm (q1)
  test_2 <- Q[index_mat == 1] > Q[index_mat == 6] # vivi->ovi|warm > |cold (q2)
  test_3 <- all(index_mat[order(Q, na.last = TRUE, decreasing = TRUE)[1:2]] == c(8,1)) #q1>q2>qij
  return(c(test_1=test_1, test_2=test_2, test_3=test_3))
}


model_test_ord <- function(cor_obj, index_mat){
  Q <- cor_obj$solution
  test_1 <- is.na(Q[1,3]) # must transition through intermediate?
  test_2 <- !is.na(Q[2,1]) & !is.na(Q[2,3]) # intermediate can go back and forward?
  test_3 <- all(is.na(Q[3,])) # final state cannot transition back 
  test_4 <- max(cor_obj$index.mat, na.rm = TRUE) == 1 # equal rates for all
  return(c(test_1=test_1, test_2=test_2, test_3=test_3, test_4=test_4))
}

model_test_hmm <- function(cor_obj, index_mat){
  Q <- cor_obj$solution
  if(dim(Q)[1] > 2){
    test_1 <- TRUE # evidence of hidden rates
    R1 <- mean(Q[c(1,2), c(1,2)], na.rm = TRUE) 
    R2 <- mean(Q[c(3,4), c(3,4)], na.rm = TRUE)
    test_2 <- sort(c(R1, R2))[2]/sort(c(R1, R2))[1] >= 100 # fast rate class 100x faster
    fastest_rate <- sort(MatrixToPars(cor_obj), decreasing = TRUE)[1]
    rate_index <- which(Q == fastest_rate, arr.ind = TRUE)
    test_3 <- !any(abs(rate_index[,1] - rate_index[,2]) > 1) # fastest rate is observed trans
    R1_sd <- sd(cor_obj$index.mat[c(1,2), c(1,2)], na.rm = TRUE)
    R2_sd <- sd(cor_obj$index.mat[c(3,4), c(3,4)], na.rm = TRUE)
    test_4 <- (R1_sd == 0) & (R2_sd == 0) # equal rates for both rate classes
  }else{
    test_1=test_2=test_3=test_4 <- FALSE
  }
  return(c(test_1=test_1, test_2=test_2, test_3=test_3, test_4=test_4))
}


fit_models <- function(data_list, pen_type, lambda, grad, file_name, overwrite, mccores) {
  if (!file_exists(file_name) || overwrite) {
    res <- mclapply(data_list, function(x) 
      corHMM:::corHMMDredgeBase(x$phy, x$cor_dat, 1, pen.type = pen_type, 
        lambda = lambda, root.p = "maddfitz", grad = grad), 
      mc.cores = mccores)
    saveRDS(res, file = file_name)
    return(res)
  } else {
    return(readRDS(file_name))
  }
}

file_exists <- function(file_path) {
  file_path %in% dir(dirname(file_path))
}

get_full_path <- function(folder_name) {
  return(normalizePath(file.path(getwd(), folder_name), winslash = "/", mustWork = FALSE))
}
