#' Retrieve all nodes in table of DDNs corresponding to pathway_name
#'
#' DDN table should have \code{node_src} and \code{node_dst}, and
#'   Pathway columns.
#'
#' @param ddns_tbl a table of DDNs
#' @param pathway_name Name of pathway to retrieve the nodes for
#' @return array of nodes
#'
#' Obsolete????  Not used in anywhere else.
get_node_set <- function(ddns_tbl, pathway_name) {
  ddns_tbl %>%
    filter(Pathway == pathway_name) %>%
    select(node_src, node_dst) %>%
    unlist() %>% unique()
}


# DDN graph construction -----

#' Create DDN graph from DDN table
#'
#' If there is no condition in DDN table, its considered a single-conditioned
#'   network.  It also adds mediators to the graph if given.
#'
#' @param ddn_tbl DDN table.  DDN table should have \code{node_src} and \code{node_dst}, and
#'   Pathway columns.
#' @param mediators mediators.  See \code{\link{add_mediators_to_DDN_graph}}.
#' @return DDN graph which is a ggraph object of DDN table
#' @seealso [add_mediators_to_DDN_graph()]
as_DDN_graph <- function(ddn_tbl, mediators = NULL) {
  if (!("condition" %in% colnames(ddn_tbl))) {
    ddn_tbl <- bind_cols(ddn_tbl, condition = "condition")
  }

  ddn_tbl %>%
    select(node_src, node_dst, condition, everything()) %>%
    as_tbl_graph %>%
    add_mediators_to_DDN_graph(mediators)
}


# Mediators ----






#' Add mediators to DDN graph
#'
#' @param ddn_graph a ggraph object of DDN table
#' @param mediator_tbl table made of name (mediator) and
#'   additional columns to indicate mediator types
#' @return ddn_graph amended with mediators
add_mediators_to_DDN_graph <- function(ddn_graph, mediator_tbl) {
  if (!is.null(mediator_tbl)) {

    # tidy up mediator table
    mediator_tbl %>%
      filter(mediator != "none") %>%
      distinct(name, mediator) ->
      mediator_tbl

    ddn_graph %>%
      activate(nodes) %>%
      left_join(mediator_tbl, by = "name")
  } else {
    ddn_graph
  }
}


#' Compute if nodes are specificity mediators
#'
#' @param nodes_edges edge list of a specific node (selected row of DDN table)
#' @param p_background Ec / Ec + Es where Ec is the number of condition specific edges
#'   and Es is the number of shared edges.
#' @noRd
compute_node_mediator_specificity_pbinom <- function(node_edges, p_background) {
  if (p_background < 1) {
    pbinom(q = nrow(filter(node_edges, toupper(condition) != "BOTH")),
           size = nrow(node_edges),
           prob = p_background, lower.tail = FALSE)
  } else {
    1  # dummy probability when there is only one condition in DDN
  }
}

#' Compute if nodes are specificity mediators
#'
#' @param nodes_edges edge list of a specific node (selected row of DDN table)
#' @param Ec the number of condition specific edges
#' @param Es the number of shared edges
#' @noRd
compute_node_mediator_specificity_phyper <- function(node_edges, Ec, Es) {
  if (Ec > 0) {  # condition-specific edges exist
    aNode_Ec <- nrow(filter(node_edges, toupper(condition) != "BOTH"))
    aNode_Es <- nrow(node_edges) - aNode_Ec
    phyper(q = aNode_Es, m = Es, n = Ec, k = aNode_Ec + aNode_Es)
  } else {
    1  # dummy probability when there is only one condition in DDN
  }
}



#' Compute specificity mediators for DDN
compute_DDN_mediator_specificity <- function(ddn, p_val.cutoff = 0.05) {
  nodes <- union(ddn$node_src, ddn$node_dst)

  Ec <- nrow(filter(ddn, toupper(condition) != "BOTH")) # no. of condition specific edges
  Es <- nrow(ddn) - Ec  # no. of shared edges

  purrr::map(nodes, function(aNode) {
    ddn %>%
      filter(node_src == aNode | node_dst == aNode) %>%
      compute_node_mediator_specificity_phyper(Ec, Es)
  }) %>% unlist -> phyper_vals

  purrr::map(nodes, function(aNode) {
    ddn %>%
      filter(node_src == aNode | node_dst == aNode) %>%
      compute_node_mediator_specificity_pbinom(p_background = Ec/(Es+Ec))
  }) %>% unlist -> pbinom_vals

  data.frame(name = nodes,
             mediator_specificity_phyper_val = phyper_vals,
             mediator_specificity_pbinom_val = pbinom_vals,
             mediator_specificity_phyper_val.fdr = p.adjust(phyper_vals, method = "fdr"),
             mediator_specificity_pbinom_val.fdr = p.adjust(pbinom_vals, method = "fdr")) %>%
    mutate(mediator_specificity = mediator_specificity_pbinom_val.fdr < p_val.cutoff)
}


#' Compute essentiality mediators for DDN
compute_DDN_mediator_essentiality <- function(ddn, percentile_cutoff = 0.05) {
  conditions <- unique(setdiff(ddn$condition, "BOTH"))

  ddn %>%
    filter(toupper(condition) == "BOTH" | condition == conditions[1]) %>%
    as_DDN_graph %>%
    activate(nodes) %>%
    mutate(betweenness = centrality_betweenness(directed = FALSE)) %>%
    as.data.frame -> nodes_1

  ddn %>%
    filter(toupper(condition) == "BOTH" | condition == conditions[2]) %>%
    as_DDN_graph %>%
    activate(nodes) %>%
    mutate(betweenness = centrality_betweenness(directed = FALSE)) %>%
    as.data.frame -> nodes_2

  if (nrow(nodes_2) == 0) {
    nodes_2 <- nodes_1[1, ]  # setting a dummy condition
    nodes_2$betweenness <- 0
  }

  left_join(nodes_1, nodes_2, by = "name") %>%
    replace_na(list(betweenness.x = 0,
                    betweenness.y = 0,
                    node.shape.x = 0,
                    node.shape.y = 0)) %>%
    mutate(diff = betweenness.x - betweenness.y) %>%
    mutate(abs_diff = abs(diff)) ->
    nodes_mediator_essentiality

  cutoff.n <- quantile(nodes_mediator_essentiality[["diff"]], percentile_cutoff/2)
  cutoff.p <- quantile(nodes_mediator_essentiality[["diff"]], 1-percentile_cutoff/2)

  nodes_mediator_essentiality %>%
    mutate(mediator_essentiality = diff < cutoff.n | diff > cutoff.p)
}




#' Compute mediators for DDN
compute_DDN_mediators <- function(ddn, mediator.type = "both") {
  ddn_specificty <- compute_DDN_mediator_specificity(ddn)
  ddn_essentiality <- compute_DDN_mediator_essentiality(ddn)

  full_join(
    ddn_essentiality %>%
      select(name, mediator_essentiality),
    ddn_specificty %>%
      select(name, mediator_specificity),
    by = "name"
  ) %>%
    replace_na(list(
      mediator_essentiality = FALSE,
      mediator_specificity = FALSE
    )) %>%
    mutate(none = "none") %>%
    mutate(essentiality = c("none", "essentiality")[as.integer(mediator_essentiality) + 1]) %>%
    mutate(specificity = c("none", "specificity")[as.integer(mediator_specificity) + 1]) %>%
    mutate(both = c("none", "specificity", "essentiality", "dual")[as.integer(mediator_essentiality) * 2 +
                                                                     as.integer(mediator_specificity) + 1]) ->
    res
  #    mutate(is_mediator = mediator_essentiality |
  #             mediator_specificity)
  res[["mediator"]] <- res[[mediator.type]]

  # clean up and return
  res %>% select(-none, -essentiality, -specificity, -both)
}


#' Compute mediators for DDN by pathway
compute_DDN_mediators_by_pathway <- function(ddn, ...) {
  split(ddn, ddn$pathway) %>%
    lapply(
      function(a_ddn) {
        compute_DDN_mediators(a_ddn, ...)
      }) %>%
    bind_rows(.id = "pathway")
}



# DDN graph plotting -----

#' Plot DDN graph from DDN table
#'
#' @param ddn_tbl see \code{\link{as_DDN_graph}}
#' @param show_mediators If show_mediators is one of 'none' (no mediator), 'specificitiy', 'essentiality' and 'both',
#'   mediators will be computed and added to DDN graph and shown.  One can also provide a list of mediators.
#' @param show_pathway_group If TRUE, pathway names will be printed around the genes belonging to the same pathway.
#' @param ... Additional parameters passed to \code{\link{plot_DDN_graph}}
#' @return ggraph/ggplot2 object
plot_DDN_tbl_graph <- function(ddn_tbl,
                               show_mediators = "both",
                               show_pathway_group = TRUE,
                               ...) {
  if (!("condition" %in% colnames(ddn_tbl))) {
    ddn_tbl <- bind_cols(ddn_tbl, condition = "C1")
  }

  ddn_tbl %>%
    mutate(condition = toupper(condition)) ->
    ddn_tbl

  ddn_graph <- ddn_tbl %>% as_DDN_graph()

  if (!is.null(show_mediators)) {
    if (is.na(show_mediators)) {
      break
    }

    if (is.vector(show_mediators)) {
      if (length(show_mediators) == 1) {
        if (show_mediators == "none") {
          break
        } else {
          ddn_mediators <- compute_DDN_mediators(ddn_tbl, show_mediators)
        }
      } else {
        ddn_mediators <- data.frame(name = show_mediators,
                                    mediator = "mediator")
      }
    } else {
      ddn_mediators <- show_mediators
    }

    ddn_graph <-
      ddn_graph %>% add_mediators_to_DDN_graph(ddn_mediators)
  }

  pathway_nodeset <- NULL

  if (show_pathway_group & "pathway" %in% colnames(ddn_tbl)) {
    ggraph_data <- ddn_graph$data

    ddn_tbl %>%
      select(pathway, node_src, node_dst) %>%
      pivot_longer(cols = c("node_src", "node_dst"),
                   values_to = "node") %>%
      select(-name) %>%
      distinct() ->
      pathway_nodeset
  } else {
    pathway_nodeset <- NULL
  }

  plot_DDN_graph(ddn_graph = ddn_graph,
                 group2node = pathway_nodeset,
                 ...)
}



#' Plot DDN graph
#'
#' @param ddn_graph see \code{\link{as_DDN_graph}}
#' @param group2node two column table where the first column is group/pathway names, and the second column is node/gene names.
#' @param graph_layout Layout algorithm name, for example, "fr".  See \code{\link{layout_nicely}}.
#' @param show_node_label if TRUE, node/gene names are shown
#' @param show_mediator_label if TRUE, mediator names are shown
#' @return ggraph/ggplot2 object
plot_DDN_graph <- function(ddn_graph,
                           group2node = NULL,
                           graph_layout = "fr",
                           show_node_label = FALSE,
                           show_mediator_label = TRUE) {
  nodes_df <- ddn_graph %>% activate(nodes) %>% data.frame
  edges_df <- ddn_graph %>% activate(edges) %>% data.frame

  if ("mediator" %in% colnames(nodes_df)) {
    ddn_graph %>%
      activate(nodes) %>%
      mutate(name_mediator = ifelse(mediator != "none", name, NA)) ->
      ddn_graph

    nodes_df <- ddn_graph %>% activate(nodes) %>% data.frame
  } else {
    show_mediator_label <- FALSE
  }

  conditions_unique <-
    setdiff(toupper(unique(edges_df$condition)), "BOTH")
  edge_colors.c <-
    c("seagreen", "tomato", "royalblue", "orange")[1:length(conditions_unique)]
  names(edge_colors.c) <- conditions_unique
  edge_colors <-
    c("BOTH" = "gray", edge_colors.c)

  ddn_graph %>%
    ggraph(layout = graph_layout) -> gp

  if ("prior" %in% colnames(edges_df)) {
    gp +
      geom_edge_link(aes(
        color = condition,
        width = condition,
        edge_linetype = prior
      ),
      alpha = 0.75) -> gp
  } else {
    gp +
      geom_edge_link(aes(color = condition,
                         width = condition),
                     alpha = 0.75) -> gp
  }

  gp +
    scale_edge_color_manual(values = edge_colors) +
    scale_edge_linetype_manual(values = c("NONE" = "dotdash", "PRIOR" = "solid")) +
    scale_edge_width_manual(values = c("Both" = 0.5, rep(.75, 10))) -> gp

  if ("mediator" %in% colnames(nodes_df)) {
    # warning: Using size for a discrete variable is not advised.
    gp <-
      gp + geom_node_point(aes(color = mediator, shape = mediator)) +
      scale_color_manual(
        values = c(
          "none" = "gray",
          "specificity" = "yellowgreen",
          "essentiality" = "cornflowerblue",
          "dual" = "tomato"
        )
      )
  }

  if (show_node_label) {
    gp <- gp + geom_node_label(aes(label = name, fill = mediator))
  } else if (show_mediator_label) {
    gp <-
      gp + geom_node_label(aes(label = name_mediator, fill = mediator))
  }

  gp <-
    gp + scale_fill_manual(
      values = c(
        "none" = "white",
        "specificity" = "yellowgreen",
        "essentiality" = "cornflowerblue",
        "dual" = "tomato"
      )
    )

  if (!is.null(group2node)) {
    names(group2node) <- c("pathway", "name")

    group_loc_df <-
      left_join(gp$data, group2node, by = "name") %>%
      group_by(pathway) %>%
      summarise(
        p_x = mean(x),
        p_y = mean(y),
        p_size = max(x = max(x) - min(x), y = max(y) - min(y)),
        n = n()
      ) %>%
      mutate(Pathway = str_wrap(pathway, width = 20)) ->
      pathway_loc_df

    gp <- gp + geom_text_repel(
      data = pathway_loc_df,
      aes(
        x = p_x,
        y = p_y,
        point.size = p_size * sqrt(n) * 2,
        label = Pathway
      ),
      size = 4,
      hjust = 0.5,
      box.padding = 0.3,
      max.overlaps = Inf,
      force = 2,
      force_pull = 1
    )
  }

  gp
}

# END -----
