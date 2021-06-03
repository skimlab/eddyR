# Read GPU-EDDY results ----

#' Read summary table
#'
#' Use this function if EDDY results are already post-processed.
#'
#' @param path path to summary table csv file
#' @param tag tag to indicate the study.  If not given,
#'            the name of folder where summary table csv file is located.
#' @return summary table
#' @export
read_summary_table <- function(path, tag = NA) {
  if (is.na(tag)) {
    tag <- basename(dirname(path))
  }

  readr::read_csv(path, col_types = cols()) %>%
    dplyr::select(-DDN) %>%
    dplyr::rename(
      n.nodes = `# genes`,
      new.dependency = `New Dependency`,
      p.val = `p-val`,
      rewiring = Rewiring,
      essentiality.mediators = `Essentiality Mediators`,
      specificity.mediators = `Specificity Mediators`
    ) %>%
    dplyr::filter(!is.nan(p.val)) %>%
    dplyr::mutate(
      essentiality.mediators = ifelse(essentiality.mediators == "nan", NA, essentiality.mediators),
      specificity.mediators = ifelse(specificity.mediators == "nan", NA, specificity.mediators)
    ) %>%
    dplyr::mutate(
      Essentiality.Mediators = parse_mediators(essentiality.mediators),
      Specificity.Mediators = parse_mediators(specificity.mediators)
    )

  #  dplyr::bind_cols(tag = tag, tbl)
}


#' Read 'results.txt'
#'
#' Use this function if EDDY results are not post-processed yet.
#'   This will be a recommended method in the future.
#'
#' Read 'results.txt' which is a table with all statistics for each gene set run by EDDY.
#'   'results.txt' is space separated file and the format is a bit messy for those FAILED gene set.
#'   So, this function actually reads 'results.txt', remove those FAILED gene set/lines,
#'   reformat them nicely and save them into a tab-delimited file, 'results.tsv'.
#'
#' @param run.dir a path to 'results.txt'
#' @return \code{tibble} table
read_results_txt <- function(run.dir) {
  # read old, ill-formatted 'results.txt'
  x <- readLines(file.path(run.dir, "results.txt"))

  # filter out 'failed' lines
  y <- x[!grepl("failed with", x)]

  # reformat the lines
  y <- gsub("\\s+", "\t", y)

  # and add a missing column headning at the end, 'n'
  y[1] <- paste(y[1], "n", sep = "")

  # replace the first column heading with Pathway
  y[1] <- sub("^.+JS", "Pathway\tJS", y[1])

  # this works with MSigDB originated pathways
  # db.name <- unlist(strsplit(y[2], split = "_"))[1]

  # writing them back into TSV formatted file
  writeLines(y, con = file.path(run.dir, "results.tsv"))

  # now read it back
  readr::read_tsv(file.path(run.dir, "results.tsv"), col_types = readr::cols())
  # bind_cols(
  #   DB = db.name,
  #   read_tsv(file.path(run.dir, "results.tsv"), col_types = cols())
  # ) %>% relocate(
  #   DB, .after = Pathway
  # )
}


#' Read DDN edgelist
#'
#' Read DDN (*_EdgeList.txt) file.
#'
#' @param ddn_filepath a path to *_Edgelist.txt
#' @return \code{tibble} table of DDN edge list.
#'   Table has 'node_src', 'node_dst', 'condition', 'prior', and 'pathway' columns.
read_DDN_EdgeList <- function(ddn_filepath) {
  postfix <- "_EdgeList.txt"
  f <- basename(ddn_filepath)
  dplyr::bind_cols(
    readr::read_tsv(
      ddn_filepath,
      col_names = c("node_src", "node_dst", "condition", "prior"),
      col_types = readr::cols()
    ),
    pathway = sub(postfix, "", f)
  )
}



#' Read all DDNs in a folder
#'
#' @param path path to a folder with all DDNs (*_EdgeList.txt files)
#' @return Combined (\code{bind_rows}) DDN table
read_DDNs <- function(path) {
  postfix <- "_EdgeList.txt"
  files2read <- list.files(path, pattern = postfix)

  lapply(files2read,
         function(f, path) {
           read_DDN_EdgeList(file.path(path, f))
         },
         path) %>% bind_rows
}


# JSON string/file IO ----

#' Capitalize words
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

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
