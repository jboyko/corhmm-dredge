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
  p <- sapply(1:max(index_mat$rate_mat, na.rm = TRUE), function(x) 
    na.omit(c(dat$par))[na.omit(c(index_mat$rate_mat) == x)][1])
  return(p)
}

# Define the posterior distribution function
log_posterior <- function(params, tree, data, rate.cat, rate.mat){
  q_prior <- list(rate=0.1)
  lp_q <- dexp(params, q_prior$rate, log=TRUE)
  lp_like <- corHMM(tree, data, rate.cat = rate.cat, rate.mat=rate.mat, p = params, node.states = "none")$loglik
  print(params)
  print(lp_like)
  print(sum(lp_q))
  lp_posterior <- sum(lp_q) + lp_like
  if(is.nan(lp_posterior) | lp_posterior < -1e6){
    lp_posterior <- -1e10
  } 
  return(lp_posterior)
}

get_better_df <- function(df, col_nm, type, ntips){
  df <- as.data.frame(df)
  colnames(df) <- col_nm
  df <- cbind(ntips = ntips, type = type, df)
  df_longer <- df %>%
    pivot_longer(cols = col_nm, names_to = "par", values_to = "value")
  return(df_longer)
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x - width / 2,
                     xmax = x)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, 
                              xmaxv = x,
                              xminv = x + violinwidth * (xmin - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

