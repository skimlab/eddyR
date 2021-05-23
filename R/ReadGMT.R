#' Read a file in GMT format into a list
#'
#' This function reads a GMT file and creates a list of gene sets.
#'
#' @param gmt.file a text file in 'GMT' format.
#'        See 'parse_a_line_in_gmt' function for detail.
#' @param min.size The minimum gene set size to read.  Defaults to 0.
#' @param max.size The maximum gene set size to read.  Defaults to no limit.
#' @param simplify Defaults to TRUE.  See 'output' for detail.
#'        If TRUE, the return structure is simplified to a list of genesets
#'        each of which is a list of genes and the name of each geneset is
#'        the name of gene set from MSigDB file, the first element of each line.
#'        If FALSE, it returns a list of three element list each element of which
#'        is \code{name}, \code{url}, and \code{geneset}.
#' @param sep The character separating fields; name, url, genes.  Defaults to '\t'
#' @return a list of gene set.
#' @export
#' @examples
#' ReadGMT('c2.cp.biocarta.v4.0.symbols.gmt')

ReadGMT <- function(gmt.file, min.size = 0, max.size = -1, simplify = TRUE, sep = "\t") {
  gs.temp <- readLines(con.gmt <- file(gmt.file))
  close(con.gmt)
  gene.sets <- lapply(gs.temp, ParseOneLineInGMT, split = sep)

  list.name <- array("", length(gene.sets))

  # convert to named list -- this will make the access to the list much easier.
  for (ii in 1:length(gene.sets)) {
    list.name[ii] <- gene.sets[[ii]]$name
  }
  names(gene.sets) <- list.name

  # Prune out ones that are too long/short
  for (ii in length(gene.sets):1) {
    if (length(gene.sets[[ii]]$gene.set) < min.size ||
        (max.size > 0 && length(gene.sets[[ii]]$gene.set) > max.size))
      gene.sets[[ii]] <- NULL
  }

  if (simplify) {
    gene.sets <- sapply(gene.sets, "[", "gene.set")
    names(gene.sets) <- sub(".gene.set$", "", names(gene.sets))
  }

  return(gene.sets)
}

#' Parse a line in 'GMT' format
#'
#' This function parses a line of string in 'GMT' format,
#'   for example a line in a MSigDB file.
#'
#' @param one.line a string in 'GMT' format. There are assumed three main elements.
#'        Each element is 'tab' delitmeted.  The first item is a name, the second
#'        item is 'url' or something like that, and the last item is a list of genes,
#'        separated by 'space'.
#' @return a list of three elements: \code{name}, \code{url}, \code{gene.set}
#' @keywords GMT parse
#' @export
#' @examples
#' ParseOneLineInGMT('a_name\tan_url\tgene1 gene2 gene3')
#' @author Seungchan Kim (dolchan@gmail.com)
#'

ParseOneLineInGMT <- function(one.line, split=" ") {
  one.line <- unlist(strsplit(as.character(one.line), split=split))

  return(list(name = one.line[1], url = one.line[2], gene.set = one.line[-c(1:2)]))
}
