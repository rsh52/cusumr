#' Cumulative sum analysis on binary dataset.
#'
#' @description
#' The cusumr command takes in a binary vector of data and returns a dataframe of rolling cusum scores per event. The dataframe is indicative of compliance with the dataset.
#' Users can choose to include a learning phase and reset options based on default values for upper and lower decision limits as well as acceptable rates of  success and error.
#'
#' For plotting see the \link[cusumr]{cusumr_plot} command.
#'
#' @param events_outcomes vector = [n_samples of 0s and 1s]. Binary Integer values. If intergers are not 0 or 1, then the FALSE alarms should be explicitly given.
#' @param acceptable_rate numeric, default=0.2. Number between 0 and 1. Acceptable success rate of the process being monitored.
#' @param type1_error_rate numeric, default=0.1. Number between 0 and 1. A false positive or type 1 error rate
#' @param type2_error_rate numeric, default=0.2 Number between 0 and 1. A false negative or type 2 error rate
#' @param learning boolean, optional (default=True). Whether to start from learning phase or from monitoring phase.
#' @param reset boolean, optional (default=True). Whether to reset score when the score hits the decision limit in monitoring phase. If yes, the cusum score will start at zero again and restart monitoring.
#'
#' @examples
#' df <- c(0,0,0,1,1,0,1,0,1,1,1,0)
#' cusumr(df)
#' cusumr(events_outcomes = df, learning = FALSE, reset = TRUE)
#'
#'
#' @references
#' The use of the Cusum technique in the assessment of trainee competence in new procedures. Int J Qual Health Care. 2000 Oct;12(5):433-8.
#'
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/11079224}
#'
#' Cumulative sum (CUSUM) assessment and medical education: a square peg in a round hole. Anaesthesia, 2011, 66, pages 243-254.
#'
#' \url{https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-2044.2011.06692.x}
#'
#' An application of the learning curve cumulative summation test to evaluate training for endotracheal intubation in emergency medicine. Emerg Med J, 2015;32:291?..294.
#'
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/24154942}


# The 'cusumr' function
# Compute cumulative summation score in binary outcome measures.
# The cusum score started from 0.
# The score will fall with each success and rise with each failure.

# The 'cusumr_graph' function
# Draw a cusum graph.The cart is a plot of the cusumative sum of a successful or adverse outcome.
# Desired outcomes result in downward steps, and upward steps are produced when adverse outcomes are encountered.

# Parameters:
# events_outcomes: vector = [n_samples of 0s and 1s]
#   Binary Integer values. If intergers are not 0 or 1, then the FALSE alarms should be explicitly given.
# acceptable_rate: numeric, default=0.2
#   Number between 0 and 1. Acceptable success rate of the process being monitored.
# unacceptable_rate: numeric, default=0.4
#   Number between 0 and 1. Unacceptable success rate of the process being monitored.
#   If the acceptable rate is equal or greater with the unacceptable rate, then a warning of
#   "Accepted failure rate should be less than unacceptable failure rate." should be explicitly given.
# type1_error_rate: numeric, default=0.1
#   Number between 0 and 1. A false positive or type 1 error rate
# type2_error_rate: numeric, default=0.2
#   Number between 0 and 1. A false negative or type 2 error rate
# learning: boolean, optional (default=True)
#   Whether to start from learning phase or from monitoring phase.
# reset: boolean, optional (default=True)
#   Whether to reset score when the score hits the decision limit in monitoring phase.
#   If yes, the cusum score will start at zero again and restart monitoring.

# Returns:
# cusum_matrix: matrix, shape = [n_samples, 5]
# The matrix have "cusum_score", "learning_signal" and "reset_signal" for every successive events.

# Note:
# The learning curve cumulative summation (LC-CUSUM) test allows for quantitative assessments of the learning process.
# The acceptable and unacceptable ratios are determined by prior research or by expert consent for each process.
# The type 1 and type 2 error rates were usually set by 0.1 to monitor acceptable and unacceptable events at same time.

# Reference
# The use of the Cusum technique in the assessment of trainee competence in new procedures. Int J Qual Health Care. 2000 Oct;12(5):433-8.
# https://www.ncbi.nlm.nih.gov/pubmed/11079224
# Cumulative sum (CUSUM) assessment and medical education: a square peg in a round hole. Anaesthesia, 2011, 66, pages 243-254.
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-2044.2011.06692.x
# An application of the learning curve cumulative summation test to evaluate training for endotracheal intubation in emergency medicine. Emerg Med J, 2015;32:291?..294.
# https://www.ncbi.nlm.nih.gov/pubmed/24154942

# The 'cusumr' function
#' @export
cusumr <- function(events_outcomes, acceptable_rate = 0.2, unacceptable_rate = 0.4,
                     type1_error_rate = 0.1, type2_error_rate = 0.1, learning = TRUE, reset = TRUE )
{
  # Start of assertion
  # Review the parameters entered and give error message explicitly in case by using checkmate/assert calls
  events_outcomes <- as.integer(events_outcomes)
  checkmate::assert_integer(events_outcomes, lower = 0, upper = 1, any.missing = FALSE,
                            min.len = 1)
  checkmate::assert_numeric(acceptable_rate, lower = 0, upper = 1,
                            finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(unacceptable_rate, lower = 0, upper = 1,
                 finite = TRUE, any.missing = FALSE, len = 1)
  if (acceptable_rate >= unacceptable_rate) {
    stop("Accepted failure rate should be less than unacceptable failure rate.")
  }
  checkmate::assert_numeric(type1_error_rate, lower = 0, upper = 1,
                 finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(type2_error_rate, lower = 0, upper = 1,
                 finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_logical(learning, any.missing = FALSE, len = 1)
  checkmate::assert_logical(reset, any.missing = FALSE, len = 1)

  # calculation of cusum symbols and formulae from features entered.
  a <- log((1-type2_error_rate)/type1_error_rate)
  b <- log((1-type1_error_rate)/type2_error_rate)
  P <- log(unacceptable_rate/acceptable_rate)
  Q <- log((1-acceptable_rate)/(1-unacceptable_rate))
  s <- Q/(P+Q)
  h0 <- b/(P+Q)
  h1 <- a/(P+Q)

  # calculation of a vector of the score_change for events_outcome
  # cusum value decrease with each success by s, increase with each failure by (1-s)
  # score_change is a vector with a mixture of increments (with each failure) and decrements (with each sucess)
  number_events <- length(events_outcomes)
  score_change <- rep(0,number_events)
  score_change[events_outcomes == 1] <- -s
  score_change[events_outcomes == 0] <- (1 - s)

  # calculation cusum score
  # the score start from zero
  cusum_score <- 0
  cusum_matrix <- matrix(0, nrow = number_events, ncol = 5)
  cusum_matrix[, 2] <- events_outcomes

  # learning signal on-off
  # if "learning == TRUE", the learning signal starts from 1
  # it will be assumed that the process is not performed successfully at the begining
  # if "learning == FALSE", the learning signal starts from 0
  # we will skip the learning phase and start from the monitoring phase
  if (learning == TRUE) {
    cusum_matrix[, 4] <- 1
  } else {
    cusum_matrix[, 4] <- 0
  }

  for (i in 1:number_events) {
    # learning phase start
    # if learnning signal is on, the score change between holding barrier (0) and decision limit (-h0).
    # holding barrier prevents the score from drifting too far away from the decision limit.
    if (cusum_matrix[i, 4] == 1) {
      cusum_score <- min(0, cusum_score + score_change[i])
      cusum_matrix[i, 3] <- cusum_score
      # when the score fall below the decision limit, the process can be considered to be competently performed.
      # then the learnig signal will be off, the monitoring phase will start.
      if (cusum_score <= -h0) {
        cusum_matrix[(i+1):number_events, 4] <- 0
      }
    }
    # monitoring phase start
    # if learnig signal is off, the score change between holding barrier (0) and decision limit (h1).
    # holding barrier prevents the score from drifting too far away from the decision limit.
    if (cusum_matrix[i, 4] == 0) {
      cusum_score <- max(0, cusum_score + score_change[i])
      cusum_matrix[i, 3] <- cusum_score
      if (cusum_score >= h1) {
      # when the score rise above the decision limit, a signal is given for the unacceptable adverse outcome rate.
      # when this happens, the monitoring process should be stopped to allow for appropriate action to be taken.
      # After this, monitoring will re-start
        if (reset == TRUE) {
          cusum_matrix[i, 5] <- 1
          cusum_score <- 0
        }
      }
    }
    cusum_matrix[i, 1] <- i
  }
  cusum_matrix <- as.data.frame(cusum_matrix)
  names(cusum_matrix) <- c("number_events", "events_outcomes", "cusum_score", "learning_signal",
                           "reset_signal")
  class(cusum_matrix) <- c("cusum_matrix", "data.frame")
  return(cusum_matrix)
}
