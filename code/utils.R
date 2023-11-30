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

get_sim_data <- function(phy, rate_mat, index_mat){
  eq_freq <- c(Null(rate_mat))
  dat <- corHMM:::simMarkov(phy, rate_mat, eq_freq)
  return(dat)
}

get_formatted_data <- function(dat, index_mat){
  tip_data <- dat$TipStates
  tip_states <- sapply(tip_data, function(x) colnames(index_mat$full_rate_mat)[x])
  tip_states <- gsub(",.*", "", tip_states)
  tip_states <- gsub("\\(", "", tip_states)
  tip_states <- as.numeric(gsub("\\)", "", tip_states))
  obs_states <- index_mat$legend[tip_states]
  cor_dat <- data.frame(sp=names(tip_data), do.call(rbind, strsplit(obs_states, "_")))
  rownames(cor_dat) <- NULL
  colnames(cor_dat) <- c("sp", 1:index_mat$nChar)
  return(cor_dat)
}

get_solution_from_res <- function(res){
  p <- sapply(1:max(res$index.mat, na.rm = TRUE), function(x) 
    na.omit(c(res$solution))[na.omit(c(res$index.mat) == x)][1])
  return(p)
}

# Define the posterior distribution function
log_posterior <- function(params, tree, data, rate.cat){
  q_prior <- list(rate=0.1)
  lp_q <- dexp(params, q_prior$rate, log=TRUE)
  lp_like <- corHMM(tree, data, rate.cat = rate.cat, model = "ARD", p = params, node.states = "none")$loglik
  print(params)
  print(lp_like)
  lp_posterior <- sum(lp_q) + lp_like
  if(is.nan(lp_posterior)){
    lp_posterior <- -1e10
  } 
  return(lp_posterior)
}
