#' Post process EDDY-GPU outputs
#'
#' @param eddy_run_dir ...
#' @param create_ddn_graph ...
#' @param mapping_samples_to_condition ...
#' @param db_name_prefix ...
#' @return ...
post_proc_EDDY <-
  function(eddy_run_dir,
           create_ddn_graph = FALSE,
           mapping_samples_to_condition = NULL,
           db_name_prefix = NA) {

    results_tbl <- read_results_txt(eddy_run_dir) %>% arrange(P)
    # (results_tbl <- mutate(results_tbl, Pathway = reactome_cleanup(Pathway)))

    DDNs <- read_DDNs(eddy_run_dir) %>%
      filter(condition != "Neither")

    if (!is.null(mapping_samples_to_condition)) {
      DDNs <-
        mutate(DDNs, condition = mapping_samples_to_condition[condition])
    }

    DDNs %>%
      group_by(pathway) %>%
      summarise(
        n_actual = length(union(node_src, node_dst)),
        known_dependency = sum(prior != "NONE") / n(),
        rewiring = sum(toupper(condition) != "BOTH") / n()
      ) %>%
      left_join(results_tbl, by = c("pathway" = "Pathway")) %>%
      rename(n_nodeset = n) %>%
      arrange(P) ->
      results_tbl

    DDNs %>%
      compute_DDN_mediators_by_pathway() -> mediator_tbl

    mediator_tbl %>%
      flatten_DDN_mediators() %>%
      left_join(results_tbl, ., by = "pathway") ->
      results_tbl

    eddy_post_proc_list <- list(summary = results_tbl,
                                DDNs = DDNs)

    if (!is.na(db_name_prefix)) {
      eddy_post_proc_list[["summary"]] %>%
        dplyr::mutate(pathway_orig = pathway) %>%
        dplyr::mutate(pathway = db_name_cleanup(pathway_orig, db_name = db_name_prefix)) %>%
        dplyr::relocate(pathway, .before = n_actual) %>%
        dplyr::relocate(pathway_orig, .after = last_col()) ->
        eddy_post_proc_list[["summary"]]

      eddy_post_proc_list[["DDNs"]] %>%
        dplyr::mutate(pathway_orig = pathway) %>%
        dplyr::mutate(pathway = db_name_cleanup(pathway_orig, db_name = db_name_prefix)) ->
        eddy_post_proc_list[["DDNs"]]
    }

    if (create_ddn_graph) {
      eddy_post_proc_list[["list_DDN_graph"]] <-
        split(eddy_post_proc_list[["DDNs"]],
              eddy_post_proc_list[["DDNs"]]$pathway) %>%
        lapply(
          FUN = function(ddn) {
            plot_DDN_tbl_graph(
              ddn,
              show_mediators = "both",
              show_node_label = FALSE,
              show_pathway_group = FALSE
            ) +
              ggtitle(ddn$pathway[1])
          }
        )

      eddy_post_proc_list[["DDN_graph_aggregated"]] <-
        plot_DDN_tbl_graph(
          eddy_post_proc_list[["DDNs"]],
          show_mediators = mediator_tbl,
          show_node_label = FALSE,
          show_pathway_group = TRUE
        )

      eddy_post_proc_list[["DDN_graph_aggregated_1"]] <-
        plot_DDN_tbl_graph(
          eddy_post_proc_list[["DDNs"]],
          show_mediators = "both",
          show_node_label = FALSE,
          show_pathway_group = TRUE
        )
    }

    eddy_post_proc_list
  }


#' Write EDDY post-processed results
#'
#' @param eddy_post_processed ...
#' @param out_dir ...
#' @param fig.width ...
#' @param fig.height ...
#' @return ...s
write_eddy_post_processed <- function(eddy_post_processed, out_dir = ".", fig.width = 12, fig.height = 9) {
  if ("DDNs" %in% names(eddy_post_processed)) {
    readr::write_csv(eddy_post_processed$DDNs, file = file.path(out_dir, "DDNs.csv"))
  }

  if ("summary" %in% names(eddy_post_processed)) {
    eddy_post_processed$summary %>%
      select(-ends_with(".l")) %>%
      readr::write_csv(file = file.path(out_dir, "summary.csv"))
  }

  if ("list_DDN_graph" %in% names(eddy_post_processed)) {
    for (nn in names(eddy_post_processed$list_DDN_graph)) {
      ggsave(filename = file.path(out_dir, sprintf("DDN_%s.pdf", nn)),
             plot = eddy_post_processed$list_DDN_graph[[nn]],
             width = fig.width, height = fig.height)
    }
  }

  if ("DDN_graph.aggregated" %in% names(eddy_post_processed)) {
    ggsave(filename = file.path(out_dir, "aggregated_DDNs_separate_mediators.pdf"),
           plot = eddy_post_processed$DDN_graph_aggregated,
           width = 2*fig.width, height = 2*fig.height)
  }

  if ("DDN_graph_aggregated_1" %in% names(eddy_post_processed)) {
    ggsave(filename = file.path(out_dir, "aggregated_DDNs_collected_mediators.pdf"),
           plot = eddy_post_processed$DDN_graph_aggregated_1,
           width = 2*fig.width, height = 2*fig.height)
  }
}

