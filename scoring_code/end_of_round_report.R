library(synapser)
library(tidyverse)
library(optparse)

SUBMISSION_QUEUES <- list(
  "on_off" = 9614407, # 9614375,
  "dyskinesia" = 9614409, # 9614377,
  "tremor" = 9614411) # 9614376)
NULL_MODEL_SUBMISSIONS <- list(
  "on_off" = 9696769,
  "dyskinesia" = 9696770,
  "tremor" = 9696771)
ROUND_END_TIMES <- list(
  "1" = "2020-04-18 17:00:00",
  "2" = "2020-04-25 17:00:00",
  "3" = "2020-05-08 17:00:00",
  "4" = "2020-05-22 17:00:00") # All PT
ROUND_END_TIMES <- purrr::map(ROUND_END_TIMES,
                              lubridate::as_datetime, tz = "America/Los_Angeles")


#' Parse command line arguments
read_args <- function() {
  option_list <- list(
    make_option("--current-round", type = "integer",
                help = "The current round as an integer."))
  opt <- parse_args2(OptionParser(option_list = option_list))
  return(opt)
}

#' Get the scores for the null models
#'
#' @param null_model_submissions A list with phenotypes for names and
#' submission IDs for values.
#' @return A list with phenotypes for names and scores for values.
get_null_model_scores <- function(null_model_submissions) {
  null_model_scores <- purrr::map(null_model_submissions, function(submission_id) {
    sub_status <- synGetSubmissionStatus(submission_id)
    double_annotations <- reshape_annotation_list(
      sub_status$annotations["doubleAnnos"])
    return(double_annotations$sqrt_weighted_mse)
  })
  return(null_model_scores)
}

#' Get the most recent ACCEPTED submission for each submitting entity
#'
#' @param phenotype The phenotype that corresponds to the evaluation queue of interest
#' @param round_end The end date of the round as a string. Submissions
#' made after the end of the round are not considered.
#' @return A data frame containing submission information about the most
#' recently submitted submission for that team/individual (the submitting entity)
get_most_recent_submissions <- function(phenotype, round_end) {
  evaluation_name <- synGetEvaluation(SUBMISSION_QUEUES[[phenotype]])$name
  subs <- synGetSubmissionBundles(
    as.character(SUBMISSION_QUEUES[[phenotype]]), status = "ACCEPTED")$asList()
  if (length(subs) == 0) {
    return()
  }
  subs_df <- purrr::map_dfr(subs, function(sub) {
      status <- sub[[2]]
      sub <- sub[[1]]
      string_annotations <- reshape_annotation_list(status$annotations["stringAnnos"])
      if (string_annotations$validation_and_scoring_error == "false") {
        double_annotations <- reshape_annotation_list(status$annotations["doubleAnnos"])
        tibble(
          createdOn = lubridate::as_datetime(sub$createdOn,
                                             tz = "America/Los_Angeles"),
          entityId = sub$entityId,
          evaluationId = sub$evaluationId,
          evaluationName = evaluation_name,
          id = sub$id,
          name = sub$name,
          userId = sub$userId,
          score = as.double(double_annotations$sqrt_weighted_mse)) %>%
          mutate(teamId = ifelse(hasName(sub, "teamId"), sub$teamId, NA))
      } else {
        return(tibble())
      }
    }) %>%
    mutate(submittingEntity = ifelse(is.na(teamId), userId, teamId)) %>%
    filter(createdOn < round_end) %>%
    arrange(desc(createdOn)) %>%
    distinct(submittingEntity, .keep_all = T) %>%
    mutate(rank = row_number(score))
  return(subs_df)
}

#' Reshape a list of annotations a more convenient format
#'
#' @param annotations A list of lists where each inner list has a \code{key}
#' named element and a \code{value} named element.
#' @return a list where \code{key} maps to \code{value}
reshape_annotation_list <- function(annotations) {
      annotation_values <- annotations %>%
        purrr::map(~ .$value)
      annotation_names <- annotations %>%
        purrr::map(~ .$key)
      names(annotation_values) <- annotation_names
      return(annotation_values)
}

#' Produce the email body to send to the submitter
#'
#' @param score The score (sqrt_weighted_mse)
#' @param rank The rank of this participant relative to the other participants
#' as an integer.
#' @param rank_out_of The number of submissions that were scored as an integer.
#' @param phenotype The phenotype that corresponds to the evaluation queue submitted to
#' as a character string.
#' @param null_model_score The score of the null model for this \code{phenotype}.
#' @param submission_id The ID of the submission.
#' @param submission_name The name of the submission.
#' @param evaluation_queue The name of the submission queue.
#' @param current_round The round number, as an integer.
#' @param recipient The name of the submitting entitity (team or individual).
#' @return A character string containing the email body.
draft_message <- function(score, rank, rank_out_of, phenotype, null_model_score,
                          submission_id, submission_name, evaluation_queue,
                          current_round, recipient) {
  next_round_end_date <- ROUND_END_TIMES[[as.character(current_round + 1)]]
  submission_name <- ifelse(is.na(submission_name), "unnamed", submission_name)
  if (!is.null(next_round_end_date)) {
    date_as_character <- strftime(next_round_end_date, "%B %d, %R %Z")
  }
  if (score > null_model_score) {
    null_model_message <- glue::glue(
      "Unfortunately, your submission performed worse than the \\
      null model. ")
    if (!is.null(next_round_end_date)) {
      null_model_message <- glue::glue(
        null_model_message, "But there's still time to resubmit!")
    }
  } else if (score < null_model_score) {
    null_model_message <- glue::glue(
      "Your submission performed better than the null model. \\
       Way to go!")
  } else if (score == null_model_score) {
    null_model_message <- glue::glue(
      "Your submission performed equally as well as the null model. \\
      That's actually pretty impressive. \\
      I didn't think this was even possible! \\
      Did you find our null model somewhere and \\
      submit it to the queue?",
      .sep = " ")
  }
  if (!is.null(next_round_end_date)) {
    next_round_message <- glue::glue(
    "The next round ends { date_as_character }. \\
    Until then — happy submitting!

    - BEAT-PD Challenge Administrator")
  } else {
    next_round_message <- "- BEAT-PD Challenge Administrator"
  }
  message <- glue::glue("Hello { recipient },

                        The rank of your most recently scored submission \\
                        ({ submission_id }, { submission_name }) \\
                        to { evaluation_queue } is { rank } \\
                        out of { rank_out_of } other participants.

                        { null_model_message }

                        { next_round_message }")
  return(message)
}


#' Sends an email to each submitting team in a given round
#'
#' @param most_recent_submissions A dataframe with columns \code{submittingEntity}
#' and \code{message}.
email_submitting_entities <- function(most_recent_submissions, phenotype,
                                      null_model_scores, current_round) {
  rank_out_of <- nrow(most_recent_submissions)
  most_recent_submissions %>%
    select(score, rank, submission_id = id, submission_name = name,
           evaluation_queue = evaluationName, teamId, submittingEntity) %>%
    purrr::pmap(function(score, rank, submission_id, submission_name,
                         evaluation_queue, teamId, submittingEntity) {
      if (!is.na(teamId)) {
        recipient <- synGetTeam(submittingEntity)$name
      } else {
        recipient <- synGetUserProfile(submittingEntity)$firstName
      }
      message_body <- draft_message(
        score = score,
        rank = rank,
        rank_out_of = rank_out_of,
        phenotype = phenotype,
        null_model_score = null_model_scores[[phenotype]],
        submission_id = submission_id,
        submission_name = submission_name,
        evaluation_queue = evaluation_queue,
        current_round = current_round,
        recipient = recipient)
      synSendMessage(userIds = list(submittingEntity),
                     messageSubject = as.character(
                       glue::glue("BEAT-PD Round { current_round } { evaluation_queue } Results")),
                     messageBody = message_body,
                     contentType = "text")
    })
}

main <- function() {
  args <- read_args()
  if (is.null(args$options$current_round)) {
    stop("Specifying the current round is required.")
  }
  synLogin()
  null_model_scores <- get_null_model_scores(NULL_MODEL_SUBMISSIONS)
  purrr::map(names(SUBMISSION_QUEUES), function(phenotype) {
    most_recent_submissions <- get_most_recent_submissions(
      phenotype = phenotype,
      round_end = ROUND_END_TIMES[[as.character(args$options$current_round)]])
    if (is.null(most_recent_submissions)) { # no one has submitted a valid submission yet
      return()
    }
    email_submitting_entities(
      most_recent_submissions = most_recent_submissions,
      phenotype = phenotype,
      null_model_scores = null_model_scores,
      current_round = args$options$current_round)
  })
}

main()
