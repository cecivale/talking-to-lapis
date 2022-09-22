# ------------------------------------------------------------------------------
#          ---        
#        / o o \    R functions
#        V\ Y /V    Talking to lapis
#    (\   / - \     Created on 5 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA, MIT license
# ------------------------------------------------------------------------------


library(jsonlite)
library(seqinr)
library(tidyverse)

#' Query LAPIS
#' @param database: get samples from open (genbank) or gisaid database (accesskey required for gisaid)
#' @param endpoint: one of aggregated/details/aa-mutations/nuc-mutations/fasta/fasta-aligned/contributors/strain-names/gisaid-epi-isl
#' @param group: fields to group the samples for aggregated endpoint
#' @param filter: attributes list returned by lapis_filter for filtering data
#' @param version: lapis version, v1 default
#' @param accessKey: access key for accessing GISAID data
#' Code modified from lapis documentation by Chaoran Chen
lapis_query <- function(database = c("open", "gisaid"), 
                        endpoint = c("aggregated", "details", "aa-mutations", "nuc-mutations", 
                                     "fasta", "fasta-aligned", "contributors", "strain-names", 
                                     "gisaid-epi-isl"), 
                        group = NULL, filter = NULL, version = "v1", accessKey = NULL) {
  # 1. Build query
  query <- lapis_build_query(match.arg(database),
                             match.arg(endpoint),
                             group, filter, version, accessKey)
  
  # 2. Query the API
  cat("\nTalking to lapis...")
    response <- tryCatch(
      if (endpoint %in% c("fasta", "fasta-aligned")) {
        cat("Downloading secuences...done!") 
        return(seqinr::read.fasta(query))
      } else if (endpoint %in% c("strain-name", "gisaid-epi-isl")) return(read_csv(query, col_names = "sample_id"))
      else fromJSON(URLencode(query)),
      error = function(e){
        message(e) 
        message("\nError from lapis connection, see query in return value.\n")} )
  cat("done!")

  # 3. Check for errors
  errors <- response$errors
  if (length(errors) > 0) {
    stop("Lapis returned some errors")
  }
  
  # 4. Check for deprecation
  deprecationDate <- response$info$deprecationDate
  if (!is.null(deprecationDate)) {
    warning(paste0("This version of the API will be deprecated on ", deprecationDate,
                   ". Message: ", response$info$deprecationInfo))
  }

  # 5. Return the data
  return(list(data = response$data, query = query, errors = response$errors, info = response$info))
}


#' Wrapper query lapis by sequence ids to get one or several endpoints.
#' @param samples: list of sample names 
#' @param sample_id: sample identifier provided, one of gisaidEpiIsl, "genbankAccession", "sraAccession", "strain"
#' @param database: one of open (genbank) or gisaid
#' @endpoint:  one or a list of details/aa-mutations/nuc-mutations/fasta/fasta-aligned/contributors/strain-names/gisaid-epi-isl
lapis_query_by_id <- function(samples, 
                              sample_id,# = c("gisaidEpiIsl", "genbankAccession", "sraAccession", "strain"),  
                              database, endpoint, batch_size = 100, 
                              access_key = Sys.getenv("LAPIS_ACCESS_KEY")) {
  
  data <- lapply(endpoint, function(ep) {
    i = 1
    data_endpoint = list()
  
    while (i < (length(samples) + 1)) {
      attr_list <- list(samples[i:(i + batch_size - 1)])
      names(attr_list) <- sample_id
      
      lapis_data <- lapis_query(database = database,
                                endpoint = ep,
                                filter = attr_list,
                                accessKey = access_key)
      
      if (ep %in% c("fasta", "fasta-aligned")) data_endpoint <- c(data_endpoint, lapis_data)
      else {
        data_endpoint <- bind_rows(data_endpoint, lapis_data$data)
        cat("\nNumber  of samples retrieved: ", nrow(data_endpoint))
        }
      i <- i + batch_size
    }
    return(data_endpoint)
  })
  names(data) <- endpoint
  return(data)
}


#' Helper function to define filters for LAPIS query
#' @param ... All possible filtering attributes in LAPIS queries
#' @return Named list of filters for lapis_query call
lapis_filter <- function(dateFrom = NULL, dateTo = NULL, dateSubmittedFrom = NULL, dateSubmittedTo = NULL,
                         region = NULL, country = NULL, division = NULL, location = NULL,
                         regionExposure = NULL, countryExposure = NULL, divisionExposure = NULL,
                         ageFrom = NULL, ageTo = NULL, sex = NULL, host = NULL, samplingStrategy = NULL,
                         pangoLineage = NULL, nextcladePangoLineage = NULL, 
                         nextstrainClade = NULL, gisaidClade = NULL,
                         submittingLab = NULL, originatingLab = NULL,
                         nucMutations = NULL, aaMutations = NULL, 
                         nextcladeQcOverallScoreFrom = NULL, nextcladeQcOverallScoreTo = NULL,
                         nextcladeQcMissingDataScoreFrom = NULL, nextcladeQcMissingDataScoreTo = NULL,
                         nextcladeQcMixedSitesScoreFrom = NULL, nextcladeQcMixedSitesScoreTo = NULL,
                         nextcladeQcPrivateMutationsScoreFrom = NULL, nextcladeQcPrivateMutationsScoreTo = NULL,
                         nextcladeQcSnpClustersScoreFrom = NULL, nextcladeQcSnpClustersScoreTo = NULL,
                         nextcladeQcFrameShiftsScoreFrom = NULL, nextcladeQcFrameShiftsScoreTo = NULL,
                         nextcladeQcStopCodonsScoreFrom = NULL, nextcladeQcStopCodonsScoreTo = NULL,
                         genbankAccession = NULL, sraAccession = NULL, gisaidEpiIsl = NULL, strain = NULL){
  
  attr_list <- compact(as.list(environment(), all=TRUE))
  return(attr_list)
}


#' Build url query for lapis
#' @inheritParams lapis_query
#' @return LAPIS query url
lapis_build_query <- function(database, endpoint, group, filter, version, accessKey){
  # Transform filter to dataframe if list
  if (!is.null(filter)) {
    if ((is.list(filter) | is.vector(filter)) & !is_tibble(filter)) attr <- enframe(filter)
    else {
      attr <- filter
      if (is_tibble(attr)) stop("The format of filter attributes is not correct, please use lapis_filter function")
    }
    attr <- attr %>% filter(!is.na(value)) %>%
      rowwise() %>%
      mutate_if(is.list, function (x)  paste0(x, collapse= ","))
  }
  
  # Build query
  query <- paste0(
    "https://lapis.cov-spectrum.org/",
    database, "/",
    version, "/",
    "sample", "/",
    endpoint, "?",
    ifelse(!is.null(group), paste0("fields=", paste0(group), collapse = ","), ""),
    "&", ifelse(!is.null(filter), paste0(paste0(attr$name, "=", attr$value), collapse = "&"), ""),
    "&",ifelse(!is.null(accessKey), paste0("accessKey=", accessKey), "")
  )
  return(query)
}


# Test by querying number of sequences in GISAID by country
# test_by_country <- lapis_query(database = "gisaid",
#                                endpoint = "aggregated",
#                                group = "country",
#                                accessKey = Sys.getenv("LAPIS_ACCESS_KEY"))
# 
# 
# 
# t <- lapis_query_by_id(database = "gisaid",
#             endpoint =c("details", "fasta-aligned"),
#             samples =  c("EPI_ISL_7367543", "EPI_ISL_7367544", "EPI_ISL_7367542"),
#             sample_id = "gisaidEpiIsl",
#             batch_size = 1)

# working! :)


lapis_query_by_id(database = "open",
                  endpoint =c("details", "fasta-aligned"),
                  samples =  c("EPI_ISL_7367543", "EPI_ISL_7367544", "EPI_ISL_7367542"),
                  sample_id = "gisaidEpiIsl",
                  batch_size = 1)

