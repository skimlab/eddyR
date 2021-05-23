#' Parse mediator columns
#'
#' @param s string made of mediators
#' @param sep separator between mediators, DEFAULT is space/white character.
#' @return a list of mediators or a string if there is only one mediator
#' @export
#'
#' In mediator columns (essentiality or specificity),
#'   mediators are separated by \code{sep},
#'   and convert them into a list of mediators or
#'   a string if there is only one mediator.
#'
parse_mediators <- function(s, sep = "\\s+") {
  if (is.na(s[1]))
    return(NA)

  r <- lapply(s,
              function(s_)
                unlist(str_split(s_, pattern = sep)))
  if (length(r) == 1)
    r[[1]]
  else
    r
}


#' Read summary table from a CSV file
#'
#' @param path path to summary table csv file.
#' @param tag tag to indicate the study.  If not given,
#'            the name of folder where summary table csv file is located.
#' @return summary table.
read_summary_table <- function(path, tag = NA) {
  if (is.na(tag)) {
    tag <- basename(dirname(path))
  }

  read_csv(path, col_types = cols()) %>%
    select(-DDN) %>%
    rename(n.nodes = `# genes`,
           new.dependency = `New Dependency`,
           p.val = `p-val`,
           rewiring = Rewiring,
           essentiality.mediators = `Essentiality Mediators`,
           specificity.mediators = `Specificity Mediators`) %>%
    filter(!is.nan(p.val)) %>%
    mutate(essentiality.mediators = ifelse(essentiality.mediators == "nan", NA, essentiality.mediators),
           specificity.mediators = ifelse(specificity.mediators == "nan", NA, specificity.mediators)) %>%
    mutate(Essentiality.Mediators = parse_mediators(essentiality.mediators),
           Specificity.Mediators = parse_mediators(specificity.mediators))

  #  bind_cols(tag = tag, tbl)
}
