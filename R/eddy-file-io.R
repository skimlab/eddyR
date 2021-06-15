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
      essentiality_mediators = `Essentiality Mediators`,
      specificity_mediators = `Specificity Mediators`
    ) %>%
    dplyr::filter(!is.nan(p.val)) %>%
    dplyr::mutate(
      essentiality_mediators = ifelse(essentiality_mediators == "nan", NA, essentiality_mediators),
      specificity_mediators = ifelse(specificity_mediators == "nan", NA, specificity_mediators)
    ) %>%
    dplyr::mutate(
      Essentiality_Mediators = parse_mediators(essentiality_mediators),
      Specificity_Mediators = parse_mediators(specificity_mediators)
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
  ) %>%
    mutate(condition = ifelse(condition %in% c("Both", "BOTH", "both"), "common", condition))
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



