#'
#' BEAT-PD Scoring Function:
#' Metric is MSE, scored within individual,
#' and weighted by either log or sqrt
#'
#' Input: Rscript {this_script} --help
#' Output: a json file
#'
#' If the input to --submission_file meets the validation requirements
#' (link to main project wiki), a json with the format
#'
#' {"validation_and_scoring_error": [false]
#'  "sqrt_weighted_mse": [float],
#'  "log_weighted_mse": [float]}
#'
#' is returned.
#'
#' Otherwise, a json with the format
#'
#' {"validation_and_scoring_error": [true],
#'  "message": [character]}
#'
#'  is returned.
#'
#'  In the special case that we cannot even read the submission file as a csv,
#'  An additional property `problems` is returned, which is the output of readr::problems
#'
#' {"validation_and_scoring_error": [true],
#'  "message": [character],
#'  "problems": [json_formatted_character]}

library(tidyverse)
library(synapser)
library(optparse)
library(jsonlite)

TEMPLATES <- list(
  on_off = "syn21344933",
  dyskinesia = "syn21344934",
  tremor = "syn21344949")

read_args <- function() {
  option_list <- list(
    make_option("--submission_file", type = "character",
                help = "Path to submission (prediction) file."),
    make_option("--entity_type", type = "character",
                default = "org.sagebionetworks.repo.model.FileEntity",
                help = "The entity type (must be org.sagebionetworks.repo.model.FileEntity)"),
    make_option("--phenotype", type = "character",
                help = "One of 'on_off', 'tremor', 'dyskinesia'."),
    make_option("--synapse_config", type = "character",
                help = "Path to Synapse config containing credentials",
                default = "~"),
    make_option("--output_file", type = "character",
                help = "Where to output results."))
  opt <- parse_args(OptionParser(option_list = option_list))
  return(opt)
}

validate_submission <- function(submission_file, entity_type, trait) {
  result <- list()
  # Did the participant submit something unscoreable like a project?
  if (entity_type != "org.sagebionetworks.repo.model.FileEntity") {
    result$validation_and_scoring_error <- TRUE
    result$message <- paste("The submission is expected to be a Synapse file,",
                            "but instead we found a submission of type",
                            entity_type, ".")
    if (entity_type == "org.sagebionetworks.repo.model.Project") {
      result$message <- paste(result$message,
                              "Did you accidentally submit your project?")
    }
    return(result)
  }
  parsing_error_text <- "There were problems reading the submission file."
  # Can we read this file as a CSV?
  df <- tryCatch({
    df <- read_csv(submission_file)
    if (nrow(problems(df))) {
      stop(parsing_error_text)
    }
    df
  }, error = function(e) {
    result$validation_and_scoring_error <- TRUE
    error_message <- gettext(e)
    result$message <- error_message
    if (str_detect(error_message, parsing_error_text)) {
      result$message <- parsing_error_text
      result$problems <- problems(df)
    }
    return(result)
  })
  if (is.list(df) && has_name(df, "validation_and_scoring_error") && df$validation_and_scoring_error) {
    return(df) # actually 'result', our object containing error info
  }
  # Does this file have all the required columns?
  if (!("measurement_id" %in% names(df))) {
    result$validation_and_scoring_error <- TRUE
    result$message <- "Did not find the column 'measurement_id'."
    return(result)
  } else if (!("prediction" %in% names(df))) {
    result$validation_and_scoring_error <- TRUE
    result$message <- "Did not find the column 'prediction'."
    return(result)
  }
  # Is measurement_id a character string?
  if (!is.character(df$measurement_id)) {
    result$validation_and_scoring_error <- TRUE
    result$message <- "Column 'measurement_id' must be a character string."
    return(result)
  }
  # Is prediction numeric?
  if (!is.numeric(df$prediction)) {
    result$validation_and_scoring_error <- TRUE
    result$message <- "Column 'prediction' must be unquoted and numeric"
    return(result)
  }
  # Are there NA values in column `prediction`?
  if (any(is.na(df$prediction))) {
    result$validation_and_scoring_error <- TRUE
    result$message <- "Column 'prediction' must not contain any NA or missing values."
    return(result)
  }
  # Do we have all measurement_ids for this trait?
  template <- read_csv(synGet(TEMPLATES[[trait]])$path)
  missing_ids <- template %>%
    anti_join(df, by = "measurement_id")
  if (nrow(missing_ids)) {
    missing_ids_str <- str_c(missing_ids$measurement_id, collapse = "\n")
    result$validation_and_scoring_error <- TRUE
    result$message <- paste0(paste(
      "Not all required measurement_id values are present",
      "for phenotype", paste0(trait, "."),
      "The following measurement_id values are missing:\n\n"),
      missing_ids_str)
    return(result)
  }
  result$validation_and_scoring_error <- FALSE
  return(result)
}

weightedMSE<-function(filename, trait){
  # Get predictions and truth
  pred<-read_csv(filename, col_names=T)
  dat<-getTruth(trait)

  dat$pred<-pred$prediction[match(dat$measurement_id, pred$measurement_id)]

  mseframe<- dat %>% group_by(subject_id) %>% summarize(n=length(truth), mse=getMSE(truth, pred))
  res<-list(sqrt_weighted_mse=weighted.mean(mseframe$mse, sqrt(mseframe$n)), log_weighted_mse=weighted.mean(mseframe$mse, log(mseframe$n)))
  return(res)
}

getTruth<-function(trait){
  realsynid<-"syn21292051"
  cissynid<-"syn21291582"

  realsyn<-synGet(realsynid)
  real<-read.csv(realsyn$path, header=T, as.is=T)
  cissyn<-synGet(cissynid)
  cis<-read.csv(cissyn$path, header=T, as.is=T)

  truth<-rbind(cis, real)
  truth<-truth[, c("measurement_id", "subject_id", trait)]
  names(truth)[3]<-"truth"
  truth<-truth[!is.na(truth$truth),]
  return(truth)
}

getMSE<-function(x, y){
  return(mean((x-y)^2))
}

main <- function() {
  args <- read_args()
  # hacky method to login, waiting on Jira issue SYNR-1007
  synapse_config_home <- file.path(path.expand("~"), ".synapseConfig")
  if (path.expand(args$synapse_config) != synapse_config_home) {
    file.copy(args$synapse_config, synapse_config_home, overwrite=TRUE)
  }
  synLogin()
  validation <- validate_submission(submission_file = args$submission_file,
                                    entity_type = args$entity_type,
                                    trait = args$phenotype)
  if (validation$validation_and_scoring_error) {
    validation$sqrt_weighted_mse <- NA
    validation$log_weighted_mse <- NA
    write_json(validation, args$output_file, auto_unbox = TRUE)
    return()
  }
  result <- weightedMSE(args$submission_file, args$phenotype)
  result$validation_and_scoring_error <- FALSE
  write_json(result, args$output_file, auto_unbox = TRUE)
}

main()
