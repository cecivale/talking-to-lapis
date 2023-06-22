# ------------------------------------------------------------------------------
#          ---        
#        / o o \    R functions
#        V\ Y /V    Talking to lapis
#    (\   / - \     Created on 5 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA, MIT license
# ------------------------------------------------------------------------------


library(jsonlite)
library(ape)
library(tidyverse)

#' Query LAPIS
#' 
#' @description Query any of the databases, endpoints, using grouping or filtering through LAPIS. 
#' Code modified from lapis documentation by Chaoran Chen
#' 
#' @param database: get samples from open (genbank) or gisaid database (accesskey required for gisaid)
#' @param endpoint: one of aggregated/details/aa-mutations/nuc-mutations/fasta/fasta-aligned/contributors/strain-names/gisaid-epi-isl
#' @param group: fields to group the samples for aggregated endpoint
#' @param filter: attributes list returned by lapis_filter for filtering data
#' @param version: lapis version, v1 default
#' @param accessKey: access key for accessing GISAID data
#' 
#' @return List with response from LAPIS: data, query, errors and info. If endopoint is fasta, returns DNAbin.
#' @export 
#' 
#' @examples 
#' # Query number of sequences in Genbank by country
#' lapis_query(database = "open",
#'             endpoint = "aggregated",
#'             group = "country")
#' 
#' # Query mutations in Swiss sequences from July to August 2022
#' filter_mut <- lapis_filter(country = "Switzerland", dateFrom = "2022-07-01", dateTo = "2022-08-01")
#' lapis_query(database = "open",
#'             endpoint = "aa-mutation",
#'             filter = filter_mut)
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
        message("\nError from lapis connection, query:\n", query)} )
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
  if (endpoint %in% c("fasta", "fasta-aligned")) {
    print(response$errors) 
    return()
    }
  return(list(data = response$data, query = query, errors = response$errors, info = response$info))
}


#' Query LAPIS by sequence ids
#' 
#' @description Wrapper query lapis by sequence ids to get one or several endpoints for a set of samples.
#' 
#' @param samples: list of sample names 
#' @param sample_id: sample identifier provided, one of gisaidEpiIsl, "genbankAccession", "sraAccession", "strain"
#' @param database: one of open (genbank) or gisaid
#' @endpoint:  one or a list of details/aa-mutations/nuc-mutations/fasta/fasta-aligned/contributors/strain-names/gisaid-epi-isl
#' 
#' @return: Data from LAPIS response, dataframe or list of dataframes, DNAbin for sequence data.
#' @export
#' 
#' @example 
#' lapis_query_by_id(database = "open",
#'                   endpoint =c("details", "fasta-aligned"),
#'                   samples =  c("EPI_ISL_7367543", "EPI_ISL_7367544", "EPI_ISL_7367542"),
#'                   sample_id = "gisaidEpiIsl")
lapis_query_by_id <- function(samples, 
                              sample_id,# = c("gisaidEpiIsl", "genbankAccession", "sraAccession", "strain"),  
                              database, endpoint, 
                              more_filter = NULL,
                              batch_size = 100, 
                              access_key = Sys.getenv("LAPIS_ACCESS_KEY")) {
  
  data <- lapply(endpoint, function(ep) {
    i = 1
    data_endpoint = list()
  
    while (i < (length(samples) + 1)) {
      attr_list <- list(samples[i:(i + batch_size - 1)])
      names(attr_list) <- sample_id
      
      lapis_data <- lapis_query(database = database,
                                endpoint = ep,
                                filter = c(attr_list, more_filter),
                                accessKey = access_key)
        
      if (ep %in% c("fasta", "fasta-aligned")) data_endpoint <- c(data_endpoint, lapis_data)
      else {
        if (ep %in% c("nuc-mutations", "aa-mutations") & batch_size == 1) data <- lapis_data$data %>% mutate(sample = attr_list)
        else data <- lapis_data$data
        
        data_endpoint <- bind_rows(data_endpoint, data)
        cat("\nNumber  of samples retrieved: ", nrow(data_endpoint))
        }
      i <- i + batch_size
    }
    return(data_endpoint)
  })
  names(data) <- endpoint
  if (length(data) == 1) return(data[[1]])
  return(data)
}

# Query lapis
# @param samples: list of sample names to get the acknowledment table
# @param seqid: sample identifier provided, default gisaidEpiIsl
# @param output: output format of acknowledgement table: "all" or "group", default grouped by lab
lapis_acknowledgment_table <- function(samples, sample_id = "gisaidEpiIsl", batch_size = 20, output = NA) {
  
  contributors <- lapis_query_by_id(samples = samples, 
                           sample_id = sample_id,  
                           database = "gisaid", 
                           endpoint = "contributors", 
                           batch_size = batch_size) %>%
    select(gisaidEpiIsl)
    # select(-genbankAccession, -sraAccession)
  open_contributors <- lapis_query_by_id(samples = contributors$gisaidEpiIsl, 
                                sample_id = "gisaidEpiIsl",  
                                database = "open", 
                                endpoint = "contributors", 
                                batch_size = batch_size) #%>%
    # select(gisaidEpiIsl, genbankAccession, sraAccession)
    data <- left_join(contributors, open_contributors, by = "gisaidEpiIsl")
    cat("\nNumber  of samples retrieved: ", nrow(data))
    cat("\nGenbank ID available for: ", sum(!is.na(data$genbankAccession)), " samples of ", nrow(data))
  
  if (!is.na(output) & output == "grouped") {
    return(data %>% group_by(originatingLab, submittingLab, authors) %>%
                                  summarise(all_gisaid = paste(gisaidEpiIsl, collapse = ","),
                                            all_genbank = paste(genbankAccession, collapse = ","),
                                            all_ena = paste(sraAccession, collapse = ","), .groups = "drop"))
  }
  return(data)
}

#' LAPIS filter
#' 
#' @description Helper function to define filters for LAPIS query
#' @param ... All possible filtering attributes in LAPIS queries
#' @return Named list of filters for lapis_query call
#' @export
#' 
#' @example 
#' lapis_filter(country = "Switzerland", dateFrom = "2022-07-01")
lapis_filter <- function(dateFrom = NULL, dateTo = NULL, dateSubmittedFrom = NULL, dateSubmittedTo = NULL,
                         region = NULL, country = NULL, division = NULL, location = NULL,
                         regionExposure = NULL, countryExposure = NULL, divisionExposure = NULL,
                         ageFrom = NULL, ageTo = NULL, sex = NULL, host = NULL, samplingStrategy = NULL,
                         pangoLineage = NULL, nextcladePangoLineage = NULL, 
                         nextstrainClade = NULL, gisaidClade = NULL,
                         submittingLab = NULL, originatingLab = NULL,
                         nucMutations = NULL, aaMutations = NULL, variantQuery = NULL, 
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

lapis_query_sequence_reference <- function(database = c("open", "gisaid"), 
                                           endpoint = c("aa-mutations", "nuc-mutations"), 
                                           filter = NULL, version = "v1", accessKey = NULL,
                                           min_proportion = 0.75,
                                           reference_genome,
                                           output_fasta = NULL) {
  mutations <- lapis_query(database = database, 
                           endpoint = endpoint, 
                           filter = filter,
                           version = version, 
                           accessKey = accessKey)$data
  mutations_filtered <- mutations %>% filter(proportion >= min_proportion) %>%
    mutate(ref = tolower(str_sub(mutation, 1, 1)),
           pos = str_sub(mutation, 2, -2),
           alt = tolower(str_sub(mutation, -1, -1)))
  
  reference <- ape::as.character.DNAbin(ape::read.FASTA(reference_genome))
  reference_mutated <- reference
  
  for (row in 1:nrow(mutations_filtered)) { 
    if (reference_mutated[[1]][as.numeric(mutations_filtered[row, "pos"])] == mutations_filtered[row, "ref"]) {
      reference_mutated[[1]][as.numeric(mutations_filtered[row, "pos"])] <- mutations_filtered[row, "alt"]
    } else { 
      warning("Error: Reference sequence different to lapis reference.")
      break
      }
  }
  
  reference_mutated_dnabin <- ape::as.DNAbin(reference_mutated)
  names(reference_mutated_dnabin) <- paste0("lapis_", paste0(filter, collapse = "_"))
    
  if (!is.null(output_fasta)) 
    write.FASTA(reference_mutated_dnabin, output_fasta)
  
  return(reference_mutated_dnabin)
  }
