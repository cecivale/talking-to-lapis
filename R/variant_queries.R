

# Variant queries examples, mutation or not mutation (without Ns)

library(tidyverse)
library(lubridate)
source("~/Tools/talking-to-LAPIS/R/lapis_functions.R")

df_mut <- lapis_query_by_id(database = "gisaid",
                            endpoint = "fasta-aligned",
                            samples = df_join$gisaid_id,
                            sample_id = "gisaidEpiIsl",
                            more_filter = lapis_filter(variantQuery = "29557T"),
                            access_key =  Sys.getenv("LAPIS_ACCESS_KEY"))

df_wt <- lapis_query_by_id(database = "gisaid",
                           endpoint = "fasta-aligned",
                           samples = df_join$gisaid_id,
                           sample_id = "gisaidEpiIsl",
                           more_filter = lapis_filter(variantQuery = "!maybe(29557T)"),
                           access_key =  Sys.getenv("LAPIS_ACCESS_KEY"))