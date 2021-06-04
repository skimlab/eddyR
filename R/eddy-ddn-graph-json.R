# JSON string/file IO ----


#' Wrap JSON element with either "{ }" or "[ ]" pair
wrap_json_elem <- function(x, wrapper = c("{", "[")) {
  wrapper <- wrapper[1]

  if (wrapper == "{") {
    w_open <- "{"
    w_close <- "}"
  } else {
    w_open <- "["
    w_close <- "]"
  }

  paste(w_open, x, w_close)
}

#' Add tag/name to JSON element
add_elem_name_to_json <- function(x, elem_name = "") {
  paste("\"", elem_name, "\":", x, sep = "")
}


#' Convert DDN table to JSON string and save it to a file
#'
#' @param ddn_tbl data frame: DDN table
#' @param filepath string: path to a file or folder.
#'        If not given, save to a current folder with pathway name as file name.
#' @param return_JSON logical (default FALSE): if true, JSON string is returned
ddn_to_json <- function(ddn_tbl, filepath = NULL, return_JSON = FALSE) {

  if (is.null(filepath)) {
    filepath <- paste(ddn_tbl$pathway[1], ".json", sep = "")
  }

  # It's an output folder, not a filename
  if (!grepl(".json$", tolower(filepath))) {
    filepath <- paste(file.path(filepath, paste(ddn_tbl$pathway[1], ".json", sep = "")))
  }

  # Generate gene to pathway mapping, to group nodes/genes to pathway
  ddn_tbl %>%
    group_by(pathway) %>%
    summarise(genes.l = list(unique(union(node_src, node_dst))),
              n = length(unique(union(node_src, node_dst)))) ->
    x

  y <- x$genes.l
  names(y) <- x$pathway

  z <- stack(y)
  names(z) <- c("gene", "pathway")

  # Map genes to pathway.
  #   If multiple pathways, pick the one with the smallest pathway
  left_join(z, select(x, pathway, n), by = "pathway") %>%
    arrange(n, pathway) %>%
    distinct(gene, .keep_all = T) ->
    gene_to_pathway

  # Create JSON "nodes" elements to connect genes/nodes to pathways
  gene_to_pathway %>%
    distinct(pathway) %>%
    mutate(id = pathway) %>%
    apply(MARGIN = 1,
          function(a_row) {
            a_row %>%
              as.list %>%
              jsonlite::toJSON(auto_unbox = TRUE) %>%
              sprintf("{ \"data\": %s }", .)
          }) %>%
    paste(collapse = ", ") ->
    my_pathway_json

  # convert condition to [C1, C2, and Both]
  condition_names <- c("Both", setdiff(unique(ddn_tbl$condition), c("BOTH", "Both", "both")))
  map_condition <- c("Both", "C1", "C2")[1:length(condition_names)]
  reverse_map_condition <- condition_names
  names(reverse_map_condition) <- map_condition
  names(map_condition) <- condition_names

  # convert to JSON string to be saved in "metadata" slot
  reverse_map_condition %>%
    as.list %>%
    toJSON() %>%
    add_elem_name_to_json(elem_name = "data") %>%
    wrap_json_elem() %>%
    wrap_json_elem("[") ->
    my_metadata_json

  ddn_tbl %>%
    mutate(condition.1 = condition) %>%
    mutate(condition = map_condition[condition.1]) %>%
    mutate(source = node_src, target = node_dst) ->
    my_ddn_tbl

  # Compute mediators and create node table
  #   in preparation to convert it JSON string
  my_ddn_tbl %>%
    compute_DDN_mediators() %>%
    #  compute_DDN_mediators_by_pathway() %>%
    mutate(id = name) %>%
    mutate(mediator = ifelse(mediator == "none", "",
                             paste(capwords(mediator), "mediator"))) %>%
    distinct(id, .keep_all = T) %>%
    #  select(-pathway) %>%
    left_join(gene_to_pathway, by = c("id" = "gene")) %>%
    mutate(parent = pathway) ->
    node_tbl

  # Now conversion to JSON string
  my_nodes <- list()
  for (i in seq(nrow(node_tbl))) {
    node_tbl[i, ] %>%
      as.list %>%
      jsonlite::toJSON(auto_unbox = TRUE) -> x

    my_nodes[[i]] <- sprintf("{\"data\": %s}", x)
  }

  my_nodes_json <- paste(my_nodes, collapse = ",")

  # Combine nodes/genes JSON string and pathway grouping JSON string
  paste(
    my_pathway_json,
    my_nodes_json,
    sep = ", "
  ) %>%
    wrap_json_elem("[") ->
    my_nodes_json

  # Conversion to "edges" JSON string
  eles <- list()
  for (i in seq(nrow(my_ddn_tbl))) {
    my_ddn_tbl[i,] %>%
      as.list %>%
      jsonlite::toJSON(auto_unbox = TRUE) -> x

    eles[[i]] <- sprintf("{\"data\": %s}", x)
  }

  paste(eles, collapse = ",") %>%
    wrap_json_elem("[") ->
    my_eles_json

  # Combine "nodes" and "edges" to a single JSON string
  paste(
    my_metadata_json %>%
      add_elem_name_to_json(elem_name = "metadata"),
    my_nodes_json %>%
      add_elem_name_to_json(elem_name = "nodes"),
    my_eles_json %>%
      add_elem_name_to_json(elem_name = "edges"),
    sep = ", "
  ) %>%
    wrap_json_elem() %>%
    prettify() ->
    my_JSON_string

  # Writing JSON string to a file
  write(my_JSON_string, file = filepath)

  if (return_JSON) {
    return(my_JSON_string)
  }
}


