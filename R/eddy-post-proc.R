# EDDY Post-processing -----


#' Post process EDDY DDNs and summary table
#'
#' @param eddy_run_dir string (path)
#' @param glasso_DDN_list list of glasso DDNs
#' @param summary_tbl data.frame
#' @param DDNs data.frame (DDN table)
#' @param create_ddn_graph logical (default is FALSE)
#' @param mapping_samples_to_condition named vector
#' @param db_name_prefix string, like REACTOME
#' @return list of 1) summary tabla, 2) DDNs (data frame),
#'         3) a list of DDN graphs,
#'         4) aggregated DDN graph (full), and
#'         5) aggregated DDN graph (p < 0.05)
post_proc_EDDY_DDNs <-
  function(DDNs,
           summary_tbl,
           create_ddn_graph = FALSE,
           mapping_conditions = NA,
           db_name_prefix = NA) {

    conditions <- setdiff(unique(DDNs$condition),
                          c("common", "Common", "COMMON", "both", "Both", "BOTH"))

    if (!isTRUE(is.na(mapping_conditions))) {
      DDNs <-
        DDNs %>%
        mutate(condition = mapping_conditions[condition]) %>%
        mutate(condition = ifelse(is.na(condition), "common", condition))

      conditions <- mapping_conditions[conditions]
    }
    # DDNs (table done)

    mediator_tbl <-
      compute_DDN_mediators_by_pathway(DDNs)

    mediator_tbl %>%
      flatten_DDN_mediators() %>%
      left_join(summary_tbl, ., by = "pathway") ->
      summary_tbl

    eddy_post_proc_list <- list(
      summary = summary_tbl,
      conditions = conditions,
      DDNs = DDNs
    )

    if (!isTRUE(is.na(db_name_prefix))) {
      eddy_post_proc_list[["summary"]] %>%
        dplyr::mutate(pathway_orig = pathway) %>%
        dplyr::mutate(pathway = db_name_cleanup(pathway_orig, db_name = db_name_prefix)) %>%
        dplyr::relocate(pathway, .before = n_nodes) %>%
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

      # Full DDNs
      eddy_post_proc_list[["DDN_graph_aggregated"]] <-
        plot_DDN_tbl_graph(
          eddy_post_proc_list[["DDNs"]],
          show_mediators = "both",
          show_node_label = FALSE,
          show_pathway_group = TRUE
        )

      # DDNs, prob < 0.05
      eddy_post_proc_list[["summary"]] %>%
        filter(prob < 0.05) %>%
        select(pathway) %>%
        unlist() %>% unname() ->
        pathway_list

      eddy_post_proc_list[["DDNs"]] %>%
        filter(pathway %in% pathway_list) %>%
        plot_DDN_tbl_graph(
          show_mediators = "both",
          show_node_label = FALSE,
          show_pathway_group = TRUE
        ) ->
        eddy_post_proc_list[["DDN_graph_aggregated_p0_05"]]

    }

    eddy_post_proc_list
  }


#' Post process EDDY-GPU outputs
#'
#' @rdname post_proc_EDDY_DDNs
post_proc_EDDY_folder <-
  function(eddy_run_dir,
           create_ddn_graph = FALSE,
           mapping_conditions = NA,
           db_name_prefix = NA) {
    results_tbl <-
      read_results_txt(eddy_run_dir) %>%
      dplyr::rename(prob.raw = P) %>%
      mutate(prob.adj = p.adjust(prob.raw)) %>%
      mutate(prob = prob.raw) %>%    # for EDDY-GPU run, we will use prob, instead of prob.adj
      arrange(prob)

    # (results_tbl <- mutate(results_tbl, Pathway = reactome_cleanup(Pathway)))

    DDNs <- read_DDNs(eddy_run_dir) %>%
      filter(condition != "Neither")

    conditions <- setdiff(unique(DDNs$condition),
                          c("common", "Common", "COMMON", "both", "Both", "BOTH"))

    DDNs %>%
      group_by(pathway) %>%
      summarise(
        n_nodes = length(union(node_src, node_dst)),
        n_edges_C1 = sum(tolower(condition) != conditions[1]),
        n_edges_C2 = sum(tolower(condition) != conditions[2]),
        n_edges_common = sum(tolower(condition) == "common"),
        rewiring = sum(tolower(condition) != "common") / n(),
        known_dependency = sum(prior != "NONE") / n()
      ) %>%
      left_join(results_tbl, by = c("pathway" = "Pathway")) %>%
      rename(n_nodeset = n) %>%
      filter(!is.na(prob)) %>%
      arrange(prob) ->
      summary_tbl

    post_proc_EDDY_DDNs(
      DDNs = DDNs,
      summary_tbl = summary_tbl,
      create_ddn_graph = create_ddn_graph,
      mapping_conditions = mapping_conditions,
      db_name_prefix = db_name_prefix
    )
  }

#' Post process EDDY-glasso-DDN outputs
#'
#' @rdname post_proc_EDDY_DDNs
post_proc_EDDY_glasso <-
  function(glasso_DDN_list,
           create_ddn_graph = FALSE,
           mapping_conditions = NA,
           db_name_prefix = NA) {

    # p-value thresholds
    p_val_threshold_loose <- 0.1
    p_val_threshold_strict <- 0.05

    summary_tbl <-
      to_glasso_DDN_summary(glasso_DDN_list) %>%
      filter(prob < p_val_threshold_loose)

    # DDNs, prob < 0.10
    pathway_list <-
      summary_tbl %>%
      filter(prob < p_val_threshold_loose) %>%
      select(pathway) %>%
      unlist() %>% unname()

    DDNs <- to_glasso_DDN_tbl(glasso_DDN_list)

    DDNs %>%
      filter(pathway %in% pathway_list) -> DDNs

    post_proc_EDDY_DDNs(
      DDNs = DDNs,
      summary_tbl = summary_tbl,
      create_ddn_graph = create_ddn_graph,
      mapping_conditions = mapping_conditions,
      db_name_prefix = db_name_prefix
    )
  }


#' Write EDDY post-processed results into tables
#'
#' @param eddy_postproc ...
#' @param output_dir ...
#' @return ...s
write_eddy_postproc_csv <- function(eddy_postproc, output_dir = ".") {
  # create output_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

  if ("DDNs" %in% names(eddy_postproc)) {
    readr::write_csv(eddy_postproc$DDNs, file = file.path(output_dir, "DDNs.csv"))
  }

  if ("summary" %in% names(eddy_postproc)) {
    eddy_postproc$summary %>%
      select(-ends_with(".l")) %>%
      readr::write_csv(file = file.path(output_dir, "summary.csv"))
  }
}


#' Write EDDY DDNs into PDF
#'
#' @param eddy_postproc ...
#' @param output_dir ...
#' @param fig.width ...
#' @param fig.height ...
#' @return ...s
write_eddy_postproc_DDN_to_pdf <- function(eddy_postproc, output_dir = ".", fig.width = 12, fig.height = 9) {
  # create output_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

  if ("list_DDN_graph" %in% names(eddy_postproc)) {
    for (nn in names(eddy_postproc$list_DDN_graph)) {
      # some pathway names has "/" and need to be fixed.
      ddn_file_name <- sprintf("DDN_%s.pdf", nn)
      ddn_file_name <- gsub("/+", "_", ddn_file_name)

      ggsave(filename = file.path(output_dir, ddn_file_name),
             plot = eddy_postproc$list_DDN_graph[[nn]],
             width = fig.width, height = fig.height)
    }
  }

  if ("DDN_graph_aggregated_p0_05" %in% names(eddy_postproc)) {
    ggsave(filename = file.path(output_dir, "aggregated_DDNs_p0_05.pdf"),
           plot = eddy_postproc$DDN_graph_aggregated_p0_05,
           width = 2*fig.width, height = 2*fig.height)
  }

  if ("DDN_graph_aggregated" %in% names(eddy_postproc)) {
    ggsave(filename = file.path(output_dir, "aggregated_DDNs.pdf"),
           plot = eddy_postproc$DDN_graph_aggregated,
           width = 2*fig.width, height = 2*fig.height)
  }
}


#' Write EDDY DDNs into JSON
#'
#' @param eddy_postproc list: ...
#' @param json_dir string: ...
#' @return ...s
write_eddy_postproc_DDN_to_json <- function(eddy_postproc, json_dir = ".") {
  # create json_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(json_dir, showWarnings = TRUE, recursive = TRUE)

  #
  pathway_list <- unique(eddy_postproc$DDNs$pathway)

  for (a_pathway in pathway_list) {
    # some pathway names has "/" and need to be fixed.
    json_file_name <- sprintf("DDN_%s.json", a_pathway)
    json_file_name <- gsub("/+", "_", json_file_name)

    eddy_postproc$DDNs %>%
      filter(pathway == a_pathway) %>%
      ddn_to_json(filepath = file.path(json_dir, json_file_name))
  }

  #
  p_val_max <- max(eddy_postproc$summary$prob, na.rm = TRUE)

  eddy_postproc$summary %>%
    filter(prob < 0.05) %>%
    select(pathway) %>%
    unlist %>% unname ->
    pathway_list

  eddy_postproc$DDNs %>%
    filter(pathway %in% pathway_list) %>%
    ddn_to_json(filepath = file.path(json_dir, "aggregated_p0_05.json"))

  aggregated_markdown <-
    paste(
      sprintf("# %s", paste(eddy_postproc$conditions, collapse = " vs ")),
      "### Aggregated DDNs (could be slow)\n",
      "<a href=\"ddngraph.html?DDN=aggregated_p0_05\" target=\"_blank\">DDNs (P.val < 0.05)</a>",
      sep = "\n")

  if (p_val_max > 0.05) {
    eddy_postproc$DDNs %>%
      ddn_to_json(filepath = file.path(json_dir, "aggregated.json"))

    aggregated_markdown <-
      paste(
        aggregated_markdown,
        "| <a href=\"ddngraph.html?DDN=aggregated\" target=\"_blank\">DDNs (Full)</a>\n",
        sep = "\n")
  }

  eddy_postproc[["summary"]] %>%
    mutate(specificity_mediators_html.l = add_GeneCard_link(specificity_mediators.l, style = "html"),
           essentiality_mediators_html.l = add_GeneCard_link(essentiality_mediators.l, style = "html")) %>%
    mutate(specificity_mediators_html = paste(specificity_mediators_html.l, collapse = ", "),
           essentiality_mediators_html = paste(essentiality_mediators_html.l, collapse = ", ")) ->
    eddy_postproc[["summary"]]

  eddy_postproc[["aggregated_markdown"]] <- aggregated_markdown

  eddy_postproc
}

#' @rdname write_eddy_postproc_json
#' Deprecated
obsolete_write_eddy_glasso_postproc_DDN_to_json <- function(eddy_postproc, json_dir = ".") {
  # create json_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(json_dir, showWarnings = TRUE, recursive = TRUE)

  #
  pathway_list <- unique(eddy_postproc$DDNs$pathway)

  for (a_pathway in pathway_list) {
    # some pathway names has "/" and need to be fixed.
    json_file_name <- sprintf("DDN_%s.json", a_pathway)
    json_file_name <- gsub("/+", "_", json_file_name)

    eddy_postproc$DDNs %>%
      filter(pathway == a_pathway) %>%
      ddn_to_json(filepath = file.path(json_dir, json_file_name))
  }

  #
  p_val_max <- max(eddy_postproc$summary$prob.adj)

  eddy_postproc$summary %>%
    filter(prob.adj < 0.05) %>%
    select(pathway) %>%
    unlist %>% unname ->
    pathway_list

  eddy_postproc$DDNs %>%
    filter(pathway %in% pathway_list) %>%
    ddn_to_json(filepath = file.path(json_dir, "aggregated_p0_05.json"))

  aggregated_markdown <-
    paste(
      "# C1 vs C2\n",
      "### Aggregated DDNs (could be slow)\n",
      "<a href=\"ddngraph.html?DDN=aggregated_p0_05\" target=\"_blank\">DDNs (P.val < 0.05)</a>",
      sep = "\n")

  if (p_val_max > 0.05) {
    eddy_postproc$DDNs %>%
      ddn_to_json(filepath = file.path(json_dir, "aggregated.json"))

    aggregated_markdown <-
      paste(
        aggregated_markdown,
        "| <a href=\"ddngraph.html?DDN=aggregated\" target=\"_blank\">DDNs (Full)</a>\n",
        sep = "\n")
  }

  eddy_postproc[["summary"]] %>%
    mutate(specificity_mediators_html.l = add_GeneCard_link(specificity_mediators.l, style = "html"),
           essentiality_mediators_html.l = add_GeneCard_link(essentiality_mediators.l, style = "html")) %>%
    mutate(specificity_mediators_html = paste(specificity_mediators_html.l, collapse = ", "),
           essentiality_mediators_html = paste(essentiality_mediators_html.l, collapse = ", ")) ->
    eddy_postproc[["summary"]]

  eddy_postproc[["aggregated_markdown"]] <- aggregated_markdown

  eddy_postproc
}


#' Write summary table in markdown format
#'
#' @param eddy_postproc list: ...
#' @param output_dir string: ...
#' @return ...s
write_eddy_summary_table_markdown <- function(eddy_postproc, output_dir) {
  # create json_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

  eddy_postproc$summary %>%
    select(
      pathway,
      n_nodes,
      known_dependency,
      rewiring,
      prob,
      essentiality_mediators_html.l,
      specificity_mediators_html.l
    ) %>%
    dplyr::rename(
      n = n_nodes,
      P.value = prob,
      'known dependency (%)' = known_dependency,
      'rewiring (%)' = rewiring,
      'essentiality mediators' = essentiality_mediators_html.l,
      'specificity mediators' = specificity_mediators_html.l
    ) %>%
    mutate(DDN = add_DDN_link(pathway, style = "html")) %>%
    mutate(pathway = add_MSigDB_link(pathway, style = "html")) %>%
    arrange(P.value) %>%
    knitr::kable(format = "markdown") %>%
    gsub("[' ']+,", ",", .) ->      # to get rid of extra spaces between words
    summary_table_markdown

  eddy_postproc[["summary_table_markdown"]] <-
    c(eddy_postproc[["aggregated_markdown"]],
      summary_table_markdown)

  eddy_postproc[["summary_table_markdown"]] %>%
    write(file = file.path(output_dir, "summary_table.md"))

  eddy_postproc
}


#' Write summary table in markdown format
#'
#' @rdname write_eddy_summary_table_markdown
#' Deprecated
obsolete_write_eddy_glasso_summary_table_markdown <- function(eddy_postproc, output_dir) {
  # create json_dir.  If it already exists, show warnings.
  # shall we remove all the files if the folder exists before moving forward?
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

  eddy_postproc$summary %>%
    select(pathway, n_actual, known_dependency, rewiring, prob.adj, essentiality_mediators_html.l, specificity_mediators_html.l) %>%
    rename(n = n_actual,
           P.value = prob.adj,
           'known dependency (%)' = known_dependency,
           'rewiring (%)' = rewiring,
           'essentiality mediators' = essentiality_mediators_html.l,
           'specificity mediators' = specificity_mediators_html.l) %>%
    mutate(DDN = add_DDN_link(pathway, style = "html")) %>%
    mutate(pathway = add_MSigDB_link(pathway, style = "html")) %>%
    arrange(P.value) %>%
    knitr::kable(format = "markdown") %>%
    gsub("[' ']+,", ",", .) ->      # to get rid of extra spaces between words
    summary_table_markdown

  eddy_postproc[["summary_table_markdown"]] <-
    c(eddy_postproc[["aggregated_markdown"]],
      summary_table_markdown)

  eddy_postproc[["summary_table_markdown"]] %>%
    write(file = file.path(output_dir, "summary_table.md"))

  eddy_postproc
}


# Utility functions -----
#

#' Capitalize words
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' Create DDN link
add_DDN_link <- function(pathway, style = c("markdown", "html")) {
  style <- style[1]

  if (style == "markdown") {
    sprintf('[DDN](ddngraph.html?DDN=%s)', pathway)
  } else {
    sprintf('<a href="ddngraph.html?DDN=%s" target="_blank">DDN</a>', pathway)
  }
}

#' Create MSigDB link
add_MSigDB_link <- function(pathway, style = c("markdown", "html")) {
  style <- style[1]

  if (style == "markdown") {
    sprintf('[%s](https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html)', reactome_cleanup(pathway), pathway)
  } else {
    sprintf('<a href="https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html" target="_blank">%s</a>', pathway, reactome_cleanup(pathway))
  }
}

#' Add GeneCard link
add_GeneCard_link <- function(geneNames, style = c("markdown", "html")) {
  style <- style[1]

  lapply(geneNames,
         function(gene, style) {
           if (style == "markdown") {
             sprintf('[%s](https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s)', gene, gene)
           } else {
             sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank">%s</a>', gene, gene)
           }
         },
         style)
}
