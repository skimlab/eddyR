#' Calculate DDN score for glasso_ddn
#'
#' @param glasso_ddn object returned from \code{\link{to_glasso_DDN}}
#' @return output of \code{\link{calc_DDN_score}}
calc_glasso_DDN_score <- function(glasso_ddn) {
  calc_DDN_score(glasso_ddn$ddn_tbl)
}

#' Calculate DDN score for ddn_tbl
#'
#' @param ddn_tbl object returned from \code{\link{to_edge_list}}
#' @return list of prob, rewiring, and other parameters computed from ddn_tbl
calc_DDN_score <- function(ddn_tbl) {
  if (is.null(ddn_tbl)) {
    return(
      list(
        prob = NA,
        rewiring = NA,
        n_nodes = 0,
        n_edges_C1 = 0,
        n_edges_C2 = 0,
        n_edges_common = 0,
        n_edges = 0
      )
    )
  }

  ddn_tbl_c <-
    ddn_tbl %>%
    filter(tolower(condition) != "common")

  ddn_tbl_s <-
    ddn_tbl %>%
    filter(tolower(condition) == "common")

  conditions <- unique(ddn_tbl_c$condition)

  # number of edges in G1 and G2
  n_edges_C1 <-
    ddn_tbl %>% filter(condition == conditions[1] |
                         tolower(condition) == "common") %>% nrow()
  n_edges_C2 <-
    ddn_tbl %>% filter(condition == conditions[2] |
                         tolower(condition) == "common") %>% nrow()
  n_edges_common <- nrow(ddn_tbl_s)
  n_edges <- n_edges_C1 + n_edges_C2 - n_edges_common


  # n1 <- sum(abs(sign(p1_x)))/2
  # n2 <- sum(abs(sign(p2_x)))/2
  # q <- sum(abs(sign(p1_x*p2_x)))/2
  # k <- length(a_set)
  # n <- k*(k-1)/2

  n_nodes = length(union(ddn_tbl$node_src, ddn_tbl$node_dst))

  list(
    prob = phyper(
      q = n_edges_common,
      m = n_edges_C1,
      n = n_edges_C2,
      k = n_edges_C2
    ),
    rewiring = 1 - n_edges_common / n_edges,
    n_nodes = n_nodes,
    n_edges_c1 = n_edges_C1,
    n_edges_C2 = n_edges_C2,
    n_edges_common = n_edges_common,
    n_edges = n_edges
  )
}

#' Train GLASSO over cross-validation
#'
#' @param X data
#' @param name name for GLASSO
#' @param K k for k-fold cross validation
#' @param iter_cv iteration of CV
#' @param rholist list of rho values.  If NULL, they will be automatically created.
#' @param thr threshold for convergence in \code{\link{glasso}}
#' @param maxit passed to \code{\link{glasso}}
#' @param approx passed to \code{\link{glasso}}
#' @param penalized_diagonal passed to \code{\link{glasso}}
#' @param rho_min_ratio used to set minimum rho value
#' @param trace passed to \code{\link{glasso}}
#' @param cores number of cores for parallel processes. Use all available cores if 0.
cv_glasso <-
  function(X,
           name = "cv_glasso",
           K = 5,
           iter_cv = 10,
           rholist = NULL,
           thr = 1.0e-4,
           maxit = 1e4,
           approx = FALSE,
           penalize_diagonal = TRUE,
           rho_min_ratio = 0.01,
           trace = 0,
           progress_bar = TRUE,
           cores = 1) {


    # start time
    s_t <- Sys.time()

    c1 <- ncol(X)
    X <- prep_data_for_Kfolds(X, K)
    if (c1 > ncol(X)) {
      warning(
        sprintf(
          "%d (out of %d) features were removed to ensure %d-fold cross-validation work well...",
          c1 - ncol(X),
          c1,
          K
        )
      )
    }


    n <- nrow(X)

    s <- var(X)
    # set up rholist so that it will be consistent for all CV iteration
    if (is.null(rholist)) {
      # calculate lam.max and lam.min
      nrho = 10
      rho_max = max(abs(s))
      rho_min = rho_min_ratio * rho_max

      # calculate grid of lambda values
      rholist = 10 ^ seq(log10(rho_min), log10(rho_max), length = nrho)

      # from glassopath itself, replaced by the above code taken from CV.glasso
      # rholist <- seq(max(abs(s)) / 100, max(abs(s)), length = 10)
    }
    rholist <- sort(rholist)

    train.size <- round(n * (K - 1) / K)

    folds_idx_list <-
      lapply(1:iter_cv, function(i, n, K) {
        create_folds_idx(n, K)
      }, n = n, K = K) %>%
      unlist(recursive = FALSE)

    train_idx_list <-
      lapply(folds_idx_list,
             function(a_fold) {
               setdiff(1:n, a_fold)
             })

    # use all available cores
    if (cores == 0) {
      cores <- availableCores()
    }

    if (progress_bar)
      pb <- progressr::progressor(along = seq_along(train_idx_list))

    #
    # internal function to train glasso over CV-fold
    # BEGIN
    #
    cv_glasso_iter <- function(k) {
      X_train <- X[train_idx_list[[k]], ]
      X_test <- X[-train_idx_list[[k]], ]

      S_train <- var(X_train)
      S_test <- var(X_test)

      if (!matrixcalc::is.positive.definite(S_train)) {
        warning(
          "Covariance matrix, S_train, is not positive-definite which could render 'glasso' hung...  So, S.train is amended to be positive-definite."
        )
        w_init = to_PD_matrix(S_train, diagonal = penalize_diagonal)
        wi_init = diag(ncol(S_train))
      } else {
        w_init = S_train
        wi_init = diag(ncol(S_train))
      }

      glist <- lapply(rholist,
                      function(rho) {
                        glasso::glasso(
                          s = S_train,
                          rho = rho,
                          thr = thr,
                          maxit = maxit,
                          approx = approx,
                          penalize.diagonal = penalize_diagonal,
                          start = "warm",
                          w.init = w_init,
                          wi.init = wi_init,
                          trace = trace
                        )
                      })

      glp <- list(
        rholist = rholist,
        approx = approx,
        wi = lapply(glist, `[[`, "wi"),
        w = lapply(glist, `[[`, "w"),
        del = lapply(glist, `[[`, "del"),
        niter = lapply(glist, `[[`, "niter")
      )

      glp$errs <-
        sapply(seq_along(glp$rholist),
               function(i) {
                 # (nrow(X_train) / 2) * (sum(wi * S_test) - determinant(wi, logarithm = TRUE)$modulus[1])
                 calc_loglik(
                   n = nrow(X_train),
                   wi = glp$wi[[i]],
                   S = S_test,
                   rho = glp$rholist[i]
                 )
               })

      glp
    }
    #
    # internal function to train glasso over CV-fold
    # END
    #

    if (cores > 1) {
      max_cores <- availableCores()

      if (max_cores < cores) {
        message(
          sprintf(
            "The number of available cores (%d) is less than what is requested (%d)...",
            max_cores,
            cores
          )
        )
        cores <- max_cores
      }

      plan(multisession, workers = cores)

      glasso_trained <-
        future_lapply(seq_along(train_idx_list),
                      function(k) {
                        if (progress_bar) {
                          # report progress, before returning
                          pb()
                        }

                        cv_glasso_iter(k)
                      })
    } else {
      glasso_trained <-
        lapply(seq_along(train_idx_list),
               function(k) {
                 if (progress_bar) {
                   # report progress, before returning
                   pb()
                 }

                 cv_glasso_iter(k)
               })
    }

    names(glasso_trained) <- 1:length(train_idx_list)

    w_init <- s
    wi_init <- diag(ncol(s))

    # summary of glasso trained over CV
    cv_gls <-
      list(
        name = name,
        X = X,
        S = s,
        rho = glasso_trained[[1]]$rholist,
        thr = thr,
        maxit = maxit,
        approx = approx,
        penalize.diagonal = penalize_diagonal,
        w_init = w_init,
        wi_init = wi_init,
        log10rho = log10(glasso_trained[[1]]$rholist),
        cv_errors = lapply(glasso_trained, `[[`, "errs") %>% bind_cols()
      )

    cv_gls[["avg_errors"]] <- rowMeans(cv_gls$cv_errors)
    cv_gls[["sd_errors"]] <- apply(cv_gls$cv_errors, MARGIN = 1, sd)


    #optimal_rho <- glasso_trained[[1]]$rholist[which.min(glp[["avg_errors"]]) + 2]
    optimal_rho <-
      find_optimal_rho(
        rholist = glasso_trained[[1]]$rholist,
        avg_errors = cv_gls[["avg_errors"]],
        sd_errors = cv_gls[["sd_errors"]]
      )

    if (optimal_rho$rho_loc_min == 1 ||
        optimal_rho$rho_loc_min == length(glasso_trained[[1]]$rholist)) {
      optimal_rho$found <- FALSE
      warning(
        sprintf(
          "Mininum 'rho' is found at the boundary, loc_min = %d.  Consider extending the range of 'rho'...",
          optimal_rho$rho.loc_min
        ),
        optimal_rho$rho_loc_min
      )
    } else {
      optimal_rho$found <- TRUE
    }

    cv_gls[["optimal"]] <- optimal_rho

    # final estimate
    gls_final <-
      glasso::glasso(
        s = s,
        rho = optimal_rho$rho,
        thr = thr,
        maxit = maxit,
        approx = approx,
        penalize.diagonal = penalize_diagonal,
        start = "warm",
        w.init = w_init,
        wi.init = wi_init
      )

    cv_gls$wi <- gls_final$wi
    rownames(cv_gls$wi) <- colnames(cv_gls$wi) <- colnames(X)

    cv_gls$w <- gls_final$w
    rownames(cv_gls$w) <- colnames(cv_gls$w) <- colnames(X)

    cv_gls$loglik <-
      calc_loglik(
        n = nrow(cv_gls$X),
        wi = cv_gls$wi,
        S = cv_gls$S,
        cv_gls$optimal$rho
      )

    cv_gls$niter <- gls_final$niter

    # end time
    e_t <- Sys.time()

    cv_gls$compute_time <-
      list(start_time = s_t,
           end_time = e_t,
           duration_secs = as.integer(difftime(e_t, s_t, units = "sec")))

    class(cv_gls) <- "cv_glasso"
    return(cv_gls)
  }


#' Get glasso with rho
#'
#' @param cv_gls output of cv_glasso
#' @param rho rho value (penalty).  If NULL, use the optimal rho found from cv_glasso.
#' @return glasso
snapshot_cv_glasso <- function(cv_gls, rho = NULL) {
  if (is.null(rho)) {
    rho <- cv_gls$optimal[["rho"]]
  }

  glasso::glasso(
    s = cv_gls$S,
    rho = rho,
    thr = cv_gls$thr,
    maxit = cv_gls$maxit,
    approx = cv_gls$approx,
    penalize.diagonal = cv_gls$penalize.diagonal,
    start = "warm",
    w.init = cv_gls$w_init,
    wi.init = cv_gls$wi_init
  ) -> gls

  rownames(gls$wi) <- colnames(gls$wi) <- colnames(cv_gls$X)
  rownames(gls$w) <- colnames(gls$w) <- colnames(cv_gls$X)

  #
  # glasso is supposed to calculate loglik but somehow, it didn't.
  #
  gls$loglik <- calc_loglik(n = nrow(cv_gls$X), wi = gls$wi, S = cv_gls$S, rho = rho)

  gls
}

#' Convert adjacent matrix into edge list
#'
#' @param adj_mat adjacent matrix
#' @param diag passed to \code{\link{shrink_adj_mat}}
#' @param symmetric force adj_mat to symmetric, also passed to \code{\link{shrink_adj_mat}}
#' @param signed if TRUE, edge will be marked either "inhibitory" or "synergistic".  If FALSE, all edges will be "connected".
#' @param simplify if TRUE, all disconnected edges will be removed
#' @return edge list (data.frame)
to_edge_list <- function(adj_mat, diag = "off", symmetric = TRUE, signed = FALSE, simplify = TRUE) {
  adj_mat <- shrink_adj_mat(adj_mat, diag_elem = diag, symmetric = symmetric, signed = signed)

  if (prod(dim(adj_mat)) == 0)  # in case an empty network
    return(NULL)

  if (symmetric)
    adj_mat[upper.tri(adj_mat)] <- 0

  data.frame(node_src = rownames(adj_mat),
             adj_mat, check.names = FALSE) %>%
    pivot_longer(
      cols = -starts_with("node_src"),
      names_to = "node_dst",
      values_to = "dependency"
    ) ->
    el_df

  if (signed) {
    el_df %>%
      mutate(dir = factor(c("inhibitory", "-", "synergistic")[2+sign(dependency)])) ->
      el_df
  } else {
    el_df %>%
      mutate(dir = "connected") ->
      el_df
  }

  el_df %>%
    mutate(edge_id = paste(node_src, node_dst, sep = "_")) ->
    el_df

  if (simplify) {
    el_df %>%
      filter(dependency != 0) ->
      el_df
  }
  el_df
}

#' Convert glasso or cv_glasso into glassonet
#'
#' @param gls glassor or cv_glasso
#' @return glassonet (glasso, precision and partial correlation matrices, edge list)
to_glassonet <- function(gls) {
  # gls: output of glasso::glasso() or cv_glasso()
  l <- list(
    name = gls$name,
    glasso = gls,
    precision = gls$wi,
    pcorr = gls$wi %>% qgraph::wi2net() %>% as.matrix(),
    el_df = gls$wi %>% to_edge_list()
  )
  class(l) <- "glassonet"

  # return a list object
  l
}

#' convert list of glassonets to DDN
#'
#' @param glnet_list a list of glassonets
#' @return glasso_DDN list (a list of glassonets and ddn_tbl)
to_glasso_DDN <- function(glnet_list) {
  a_DDN <- glnet_list

  el_list <- lapply(glnet_list, `[[`, "el_df")

  a_DDN[["ddn_tbl"]] <- merge_edges_list_DDN(el_list)
  if (is.null(a_DDN[["ddn_tbl"]])) {
    warning(sprintf("DDN: %s is empty...", a_DDN[[1]]$name))
  }

  class(a_DDN) <- "DDN"

  a_DDN
}

#' Generate summary table
#'
#' @param glass_DDN_list list of glasso_DDNs
#' @return summary table
to_glasso_DDN_summary <- function(glasso_DDN_list, selected_names = NULL) {

  if (!is.null(selected_names)) {
    glasso_DDN_list <- glasso_DDN_list[selected_names]
  }

  DDN_summary <-
    lapply(glasso_DDN_list,
           function(ddn) {
             calc_glasso_DDN_score(ddn)
           }) %>% bind_rows(.id = "pathway")

  if (!('known_dependency' %in% colnames(DDN_summary))) {
    DDN_summary[["known_dependency"]] <- 0     # all statistical dependencies
  }

  DDN_summary %>%
    dplyr::rename(prob.raw = prob) %>%
    mutate(prob.adj = p.adjust(prob.raw)) %>%
    mutate(prob = prob.adj) %>%    # set prob to prob.adj for glassnoet
    arrange(prob, -rewiring)
}

#' Combine all DDNs into a single table
#'
#' @param glasso_DDN_list list of glass DDNs
#' @return data frame, DDN table
to_glasso_DDN_tbl <- function(glasso_DDN_list) {
  # remove empty DDNs
  to_keep <-
    glasso_DDN_list %>% lapply(function(ddn) {
      !is.null(ddn$ddn_tbl) && nrow(ddn$ddn_tbl) > 0
    }) %>% unlist()

  lapply(glasso_DDN_list[which(to_keep)],
         `[[`, "ddn_tbl") %>%
    bind_rows(.id = "pathway")
}


#
# internal -----
#
merge_edges_list_DDN <- function(edges_list) {
  ddn_df <-
    bind_rows(edges_list, .id = "condition")

  if ("condition" %in% colnames(ddn_df)) {
    ddn_df <- relocate(ddn_df, condition, .after = last_col())
  } else {
    return(NULL)
  }

  which_duplicated <- which(duplicated(ddn_df$edge_id))
  duplicated_el <- ddn_df$edge_id[which_duplicated]

  ddn_df %>%
    mutate(condition = ifelse(edge_id %in% duplicated_el, "common", condition)) -> ddn_df

  if (length(which_duplicated) > 0)
    ddn_df <- ddn_df[-which_duplicated, ]
  ddn_df
}

merge_edges_list_DDN_v1 <- function(edges_list) {
  el1 <- edges_list[[1]]
  el2 <- edges_list[[2]]

  c1 <- names(edges_list)[1]
  c2 <- names(edges_list)[2]

  merge_edges_DDN(el1, el2, c1, c2)
}



merge_edges_DDN <- function(el1, el2, c1, c2) {
  if (is.null(el1) & is.null(el2)) {
    warning("DDN is empty...")
    return(NULL)
  }

  el <-
    bind_rows(bind_cols(el1, condition = c1),
              bind_cols(el2, condition = c2))

  which_duplicated <- which(duplicated(el$edge_id))
  duplicated_el <- el$edge_id[which_duplicated]

  el %>%
    mutate(condition = ifelse(edge_id %in% duplicated_el, "common", condition)) ->
    el

  if (length(which_duplicated) > 0)
    el <- el[-which_duplicated, ]
  el
}


prep_data_for_Kfolds <- function(x, K) {
  # scale data
  x_scaled <- scale(x)

  nr <- nrow(x_scaled)
  n_nonzeros <- colSums(x_scaled > 0, na.rm = 0)  # count the number of non-zeros in each column
  zero_columns <- which(n_nonzeros < nr*(1/K))  # identify the column with the count < 1/K samples/rows
  na_columns <- which(colSums(is.na(x_scaled)) > 0)  # NA columns

  col2remove <- union(zero_columns, na_columns)

  if (length(col2remove) > 0)
    x <- x[, -col2remove]

  return(x)
}

create_folds <- function(samples, k) {
  n <- length(samples)

  blksz <- floor(n/k)
  leftover <- n - blksz*k

  bagsz <- rep(blksz, k)

  if (leftover > 0) {
    bagsz[1:leftover] <- bagsz[1:leftover] + 1
  }
  sx <- c(1, bagsz[-length(bagsz)]) %>% cumsum()
  ex <- bagsz %>% cumsum()

  lapply(1:length(sx),
         function(i) {
           samples[seq(sx[i], ex[i])]
         })
}

create_folds_idx <- function(n, k) {
  create_folds(permute::shuffle(1:n), k)
}


to_PD_matrix <- function(s, diagonal = TRUE, rho = max(abs(s))) {
  # specify initial estimate for Sigma
  if (diagonal) {

    # simply force init to be positive definite final diagonal
    # elements will be increased by lam
    init = s + diag(rho, nrow(s))

  } else {

    s_minus <- s
    diag(s_minus) <- 0

    # provide estimate that is pd and dual feasible
    alpha = min(c(rho/max(abs(s_minus)), 1))
    init = (1 - alpha) * s
    diag(init) = diag(s)
  }
  init
}


find_optimal_rho <-
  function(rholist,
           avg_errors,
           sd_errors,
           method = c("sd", "min_error")) {

    min_loc <- which.min(avg_errors)

    min_rho <- rholist[min_loc]
    min_sd <- sd_errors[min_loc]

    x <- avg_errors - avg_errors[min_loc]
    x[1:min_loc] <- 0

    sd_loc <- min(which(x > min_sd), na.rm = T)

    choice <- pmatch(method[1], c("sd", "min_error"))

    switch(
      choice,
      "1" = list(rholist[sd_loc], sd_loc, c("sd", "min_error")[choice]),
      "2" = list(rholist[min_loc], min_loc, c("sd", "min_error")[choice])
    ) ->
      res

    names(res) <- c("rho", "rho_loc", "method")

    c(
      res,
      list(
        rho_min = rholist[min_loc],
        rho_loc_min = min_loc,
        rho_sd = rholist[sd_loc],
        rho_loc_sd = sd_loc
      )
    )
  }



calc_loglik <- function(n, wi, S, rho, penalize_diagonal = TRUE) {
  # option to penalize diagonal
  if (penalize_diagonal) {
    C = 1
  } else {
    C = 1 - diag(ncol(S))
  }

  (n / 2) * (sum(diag(as.matrix(wi) %*% as.matrix(S))) + rho*sum(abs(C*wi)) - determinant(wi, logarithm = TRUE)$modulus[1])
}


shrink_adj_mat <- function(adj_mat, diag_elem = "off", symmetric = FALSE, signed = TRUE) {
  if (diag_elem == "off")
    diag(adj_mat) <- 0

  adj_conn <- sign(abs(adj_mat))

  which_conn <- which(rowSums(adj_conn) > 1)
  adj_mat <- adj_mat[which_conn, which_conn, drop = FALSE]

  if (symmetric) {
    if (!isSymmetric(adj_mat)) {
      warning("Adjacency matrix is not symmetric... now being forced to symmetric...")
      adj_mat <- (adj_mat + t(adj_mat))/2
    }

    # adj_mat[upper.tri(adj_mat)] <- 0
  }
  adj_mat
}

# END -----

