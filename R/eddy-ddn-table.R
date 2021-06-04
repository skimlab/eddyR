#' Clean up Pathway DB name
#'
#' Remove some prefix, like "REACTOME", specified by db_name,
#'   and postfix, like "_PATHWAY" in pathway_name
#'
#' @param pathway_name name for pathway
#' @param db_name prefix indicating pathway database, like "REACTOME"
#' @return new pathway name without prefix and postfix
db_name_cleanup <- function(pathway_name, db_name) {
  pathway_name %>%
    sub(sprintf("^%s_", db_name), "", .) %>%
    sub("_PATHWAY$", "", .) %>%
    gsub("_", " ", .) %>%
    gsub(" +$", "", .) %>%
    gsub("^ +", "", .)
}

#' @rdname db_name_cleanup
reactome_cleanup <- function(pathway_name) {
  pathway_name %>%
    db_name_cleanup(db_name = "REACTOME")
}


#' Expand a pathway line with multiple mediators
#'
#' When there are multiple mediators in a pathway,
#'   the pathway line is replicated in such a way
#'   that each line contains only one mediator and its type.
#'
#' @param p a pathway line in summary table
#' @return expanded summary table
expand_mediators <- function(p) {
  mediator_list <- list(essentiality = (p$Essentiality_Mediators %>% unlist %>% unique()),
                        specificity = (p$Specificity_Mediators %>% unlist %>% unique()))

  right_join(select(p, -c(Essentiality_Mediators, Specificity_Mediators)),
             tibble(Pathway = p$Pathway,
                    mediator = (unlist(mediator_list) %>% unique() %>% as.vector)) %>%
               mutate(essentiality = mediator %in% mediator_list$essentiality,
                      specificity = mediator %in% mediator_list$specificity) %>%
               filter(!is.na(mediator), toupper(mediator) != "NAN"),
             by = "Pathway")
}


#' Flatten expanded summary/mediator table
flatten_DDN_mediators <- function(mediator_tbl) {
  # internal function
  flatten_DDN_mediators_ <- function(mediator_tbl) {
    mediator_tbl %>%
      filter(mediator != "none") %>%
      summarise(
        essentiality_mediators.l = list(name[mediator_essentiality]),
        specificity_mediators.l = list(name[mediator_specificity])
      ) %>%
      mutate(
        essentiality_mediators = paste(unlist(essentiality_mediators.l), collapse = ", "),
        specificity_mediators = paste(unlist(specificity_mediators.l), collapse = ", ")
      ) %>%
      select(ends_with("mediators"), ends_with(".l"))
  }

  # split by pathways and flatten DDNs of each pathway
  split(mediator_tbl, mediator_tbl$pathway) %>%
    lapply(flatten_DDN_mediators_) %>%
    bind_rows(.id = "pathway")
}


#' Add GeneCard to mediators
add_GeneCard_to_mediators <- function(mediator_tbl) {

}

#' Parse mediator string
#'
#' In mediator columns (essentiality or specificity),
#'   mediators are separated by \code{sep},
#'   and convert them into a list of mediators or
#'   a string if there is only one mediator.
#'
#' @param s string made of mediators
#' @param sep separator between mediators, DEFAULT is space/white character.
#' @return a list of mediators or a string if there is only one mediator
#' @export
parse_mediators <- function(s, sep = "\\s+") {
  if (is.null(s))
    return(NA)

  r <- lapply(s,
              function(s_)
                unlist(str_split(s_, pattern = sep)))
  if (length(r) == 1)
    r[[1]]
  else
    r
}



#' Retrieve mediators as an array
#'
#' Mediators are stored as a string.  This function parse them
#'   and return them in an array.  One can retried either specificity
#'   mediators or essentiality mediators or both simultaneously.
#' @param p a pathway line in summary table
#' @param type which type of mediators to retrieve
#' @return array of mediators
get_mediators <- function(p, type = "both") {
  if (type == "both") {
    types = c("specificity", "essentiality")
  } else {
    types = type
  }

  mediators <- c()

  if ("specificity" %in% types)
    mediators <- union(mediators, unlist(parse_mediators(p$specificity_mediators)))

  if ("essentiality" %in% types)
    mediators <- union(mediators, unlist(parse_mediators(p$essentiality_mediators)))

  mm <- unique(mediators[!is.na(mediators)])
  mm[stringr::str_length(mm) > 0]
}


to_DDN_summary <- function(ddn_tbl) {

  geneset_terms <- unique(ddn_tbl$pathway)

  split(ddn_tbl, ddn_tbl$pathway) %>%
    lapply(calc_DDN_score) %>% bind_rows() %>%
    bind_cols(name = geneset_terms, .) %>%
    select(-n) %>%
    mutate(prob.adj = p.adjust(prob)) %>%
    relocate(prob.adj, .after = prob) %>%
    arrange(prob.adj,-rewiring)
}

