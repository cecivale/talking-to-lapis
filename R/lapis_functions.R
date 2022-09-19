# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Util functions
#        V\ Y /V    Talking to lapis
#    (\   / - \     5 July 2022, updated 20 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(jsonlite)
library(tidyverse)

# Function to get empty list of lapis attributes
lapis_attributes <- function(){
  names_attr_list <- read_tsv("~/Tools/utils/lapis_attributes_list.tsv", col_names = "attribute")
  attr_list <- rep(NA, nrow(names_attr_list))
  names(attr_list) <- names_attr_list$attribute
  return(as.list(attr_list))
}

# Build query for lapis
lapis_build_query <- function(database, endpoint, fields, attributes, version, accessKey){
  # Transform attributes to dataframe if list
  if (!is.null(attributes)) {
    if ((is.list(attributes) | is.vector(attributes)) & !is_tibble(attributes)) attr <- enframe(attributes)
    else {
      attr <- attributes
      if (is_tibble(attr)) stop("The format of attributes is not correct, please use lapis_attributes or by lapis_attributes_table functions")
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
    ifelse(!is.null(fields), paste0("fields=", paste0(fields), collapse = ","), ""),
    "&", ifelse(!is.null(attributes), paste0(paste0(attr$name, "=", attr$value), collapse = "&"), ""),
    "&",ifelse(!is.null(accessKey), paste0("accessKey=", accessKey), "")
  )
  return(query)
}

# Query lapis
# @param database: get samples from open (genbank) or gisaid database (accesskey required for gisaid)
# @param endpoint: one of aggregated/details/aa-mutations/nuc-mutations/fasta/fasta-aligned/contributors/strain-names/gisaid-epi-isl
# @param fields: fields to group the samples for aggregated endpoint
# @param attributes: attributes list returned by lapis_attributes or by lapis_attributes_table to filter data
# @param version: lapis version, v1 default
# @param accessKey: access key for accessing GISAID data
# Code from lapis documentation @Chaoran Chen
lapis_query <- function(database, endpoint, fields = NULL, attributes = NULL, version = "v1", accessKey = NULL) {
  # 1. Build query
  query <- lapis_build_query(database, endpoint, fields, attributes, version, accessKey)
  
  # 2. Query the API
  cat("\nTalking to lapis...")
    response <- tryCatch(
      if (endpoint %in% c("fasta", "fasta-aligned")) return(ape::read.FASTA(query))
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


# Query lapis
# @param samples: list of sample names to get the acknowledment table
# @param seqid: sample identifier provided, default gisaidEpiIsl
# @param output: output format of acknowledgement table: "all" or "group", default grouped by lab
lapis_acknowledgment_table <- function(samples, seqid = "gisaidEpiIsl", output = "group") {
  i = 1
  data = tibble()
  
  while (i < (length(samples) + 1)) {
    # Contributors from gisaid database
    attr_list <- list(samples[i:(i+40)])
    names(attr_list) <- seqid
    contributors <- lapis_query(database = "gisaid",
                                endpoint = "contributors",
                                attributes = attr_list,
                                accessKey = Sys.getenv("LAPIS_ACCESS_KEY"))$data %>%
      select(-genbankAccession, -sraAccession)
    # Get genbank ids from open database
    attr_list <- list(gisaidEpiIsl = contributors$gisaidEpiIsl)
    open_contributors <- lapis_query(database = "open",
                                     endpoint = "contributors",
                                     attributes = attr_list)$data %>%
      select(gisaidEpiIsl, genbankAccession, sraAccession)
    
    data <- bind_rows(data, left_join(contributors, open_contributors, by = "gisaidEpiIsl"))
    cat("\nNumber  of samples retrieved: ", nrow(data))
    cat("\nGenbank ID available for: ", sum(!is.na(data$genbankAccession)), " samples of ", nrow(data))
    i <- i + 41
  }
  
  if (output == "group") return(data %>% group_by(originatingLab, submittingLab, authors) %>%
                                    summarise(all_gisaid = paste(gisaidEpiIsl, collapse = ","),
                                              all_genbank = paste(genbankAccession, collapse = ","),
                                              all_ena = paste(sraAccession, collapse = ","), .groups = "drop"))
  return(data)
}


# Test by querying number of sequences in GISAID by country
# test_by_country <- lapis_query(database = "gisaid", 
#                                endpoint = "aggregated", 
#                                fields = "country",
#                                accessKey = Sys.getenv("LAPIS_ACCESS_KEY"))
# 
# 
# attr_list <- lapis_attributes()
# attr_list$gisaidEpiIsl <- "EPI_ISL_7367543"
# t <- lapis_query(database = "gisaid", 
#             endpoint = "details", 
#             attributes = attr_list,
#             accessKey = Sys.getenv("LAPIS_ACCESS_KEY"))

# working! :)




