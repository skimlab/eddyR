#' Plot CV error vs rho from cv_glasso
#'
#' @param
ggplot.cv_glasso <- function(glp) {
  plot.cv_glasso(glp)
}

plot.cv_glasso <- function(glp) {
  df <- data.frame(rho = glp$rho,
                   avg_error = glp$avg_errors,
                   sd = glp$sd_errors)
  ggplot(df,
         aes(x = rho, y = avg_error)) +
    geom_point() +
    geom_errorbar(aes(ymin=avg_error-sd, ymax=avg_error+sd),
                  position=position_dodge(0.05)) +
    geom_path() +
    geom_vline(xintercept = glp$optimal$rho,
               linetype="dotted",
               color = "blue") ->
    gp

  if (glp$optimal$method == "sd") {
    gp <- gp +
      geom_hline(yintercept = glp$sd_errors[glp$optimal$rho_loc_min] + glp$avg_errors[glp$optimal$rho_loc_min],
                 linetype="dotted",
                 color = "blue")
  }

  gp + scale_x_log10()
}


# wi is a precision matrix or
#   partial correlation matrix which is obtained by applying wi2net to a precision matrix
heatmap.glassonet <-
  function(glassonet,
           signed = TRUE,
           symmetric = TRUE,
           diag = "off",
           heatmap_size = unit(3, "inch")) {

    pcorr <- glassonet$pcorr

    pcorr %>%
      shrink_adj_mat(diag = diag,
                     symmetric = symmetric) %>%
      sign -> pcorr

    if (!signed)
      pcorr <- abs(pcorr)

    Heatmap(
      pcorr,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      show_row_names = F,
      show_column_names = F,
      width = heatmap_size,
      height = heatmap_size
    )
  }




plot.netDegree <- function(net, include_zero = FALSE) {
  net <- abs(sign(net))
  netDegree <-
    data.frame(name = rownames(net),
               degree = rowSums(net))

  if (!include_zero)
    netDegree <- filter(netDegree, degree > 0)

  ggplot(netDegree,
         aes(x = degree)) +
    geom_histogram(binwidth = 1) +
    geom_vline(aes(color = "red", xintercept = median(degree))) +
    geom_vline(aes(color = "blue", xintercept = mean(degree))) +
    geom_label(aes(x = median(degree) + 0.02*max(degree),
                   y = max(table(degree))/2,
                   hjust = 0,
                   label = sprintf("median degree = %.1f", median(degree)))) +
    geom_label(aes(x = mean(degree) + 0.02*max(degree),
                   y = max(table(degree))/3,
                   hjust = 0,
                   label = sprintf("mean degree = %.1f", mean(degree))))

}



plot_graph_glassonet <- function(glassonet, ...) {
  plot_DDN_tbl_graph(glassonet$el_df, ...)
}


plot_glasso_DDN <- function(glasso_ddn, ...) {
  plot_DDN_tbl_graph(glasso_ddn$ddn_tbl, ...)
}


pathway_annotate_group <- function(g, group_column = "group", genesets, categories = c("REACTOME")) {
  g %>% activate(nodes) %>%
    data.frame() -> g_nodes_df

  x <- g_nodes_df[[group_column]] %>% unique()
  x <- sort(x[!is.na(x)])

  sapply(x,
         function(i) {
           nn <- g_nodes_df[["name"]][which(g_nodes_df[[group_column]] == i)]
           do_enricher(gene = nn,
                       geneset_df = genesets,
                       db_categories = categories) %>%
             data.frame() -> x.df
           x.df$ID[1]
         }) ->
    grp_to_pathway

  grp_to_pathway <- data.frame(group = x, pathway = grp_to_pathway)

  g %>% activate(nodes) %>%
    left_join(grp_to_pathway, by = c(group.brief = "group"))
}

