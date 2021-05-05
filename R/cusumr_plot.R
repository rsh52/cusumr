#' Graphic display of cumulative sum analysis on binary dataset
#'
#' @description
#' The cusumr_plot() command builds on the same architecture of the \link[cusumr]{cusumr} command, while automating the production of a graphic visual. Visuals are developed using the ggplot2 package.
#'
#' @param events_outcomes vector = [n_samples of 0s and 1s]. Binary Integer values. If intergers are not 0 or 1, then the FALSE alarms should be explicitly given.
#' @param acceptable_rate numeric, default=0.2. Number between 0 and 1. Acceptable success rate of the process being monitored.
#' @param unacceptable_rate numeric, default=0.4. Number between 0 and 1.
#' @param type1_error_rate numeric, default=0.1. Number between 0 and 1. A false positive or type 1 error rate
#' @param type2_error_rate numeric, default=0.2 Number between 0 and 1. A false negative or type 2 error rate
#' @param learning boolean, optional (default=True). Whether to start from learning phase or from monitoring phase.
#' @param reset boolean, optional (default=True). Whether to reset score when the score hits the decision limit in monitoring phase. If yes, the cusum score will start at zero again and restart monitoring.
#'
#' @examples
#' df <- df <- c(0,0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
#' cusumr_plot(df)
#' cusumr_plot(events_outcomes = df, learning = FALSE, reset = TRUE)
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

# Perform cumulative sum analysis on binary dataset and plot
# The 'cusumr_plot' function
#' @export
cusumr_plot <- function(events_outcomes, acceptable_rate = 0.2, unacceptable_rate = 0.4,
                        type1_error_rate = 0.1, type2_error_rate = 0.1, learning = TRUE, reset = TRUE){
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

  # call cusum_matrix from cusumr fx
  # cusum_matrix <- cusumr(events_outcomes, acceptable_rate, unacceptable_rate,
  #                        type1_error_rate, type2_error_rate, learning, reset)

  # make range_tibble where the learning_ and reset_signals start & end
  # the tibble looks like (learning_phase=c(start,end),monitoring_1=c(start,end),last_phase=c(start,end)
  # first, let's find end points for each phase
  if (learning == TRUE){
    # if learning is never over, let's draw just one graph
    if (cusum_matrix[-(1:(length(cusum_matrix$number_events)-1)),4] == 1){
      range_end <- c(length(cusum_matrix$number_events))
    } else if (cusum_matrix[-(1:(length(cusum_matrix$number_events)-1)),4] == 0){
      # if learning is over, let's draw one learning graph + several monitoring graph
      range_end <- c(sum(cusum_matrix$learning_signal),cusum_matrix$number_events[cusum_matrix$reset_signal == 1],length(cusum_matrix$number_events))
    }
  } else if  (learning == FALSE){ # set defalut of no_graph, when learning == Not True
    range_end <- c(cusum_matrix$number_events[cusum_matrix$reset_signal == 1],length(cusum_matrix$number_events))
  }
  # let's find staring points and assemble range_tibble
  range_start <- c(1,c(range_end+1)[1:length(range_end)-1])
  range_tibble <- tibble::tibble(range_start,range_end)
  # next function will find cusum score for each phase of cusum chart
  # the matrix returned hold events numbers and each cusum scores to feed your ggplot()

  #' @export
  find_score <- function(range_tibble,i){
    starting_point <- c(as.integer(range_tibble[i,1])-1,0)
    names(starting_point) <- c("event", "event_score")
    event <- c(as.integer(range_tibble[i,1]):as.integer(range_tibble[i,2]))
    event_score <- cusum_matrix$cusum_score[event]
    event_matrix <- tibble::tibble(event, event_score)
    range_matrix <- rbind(starting_point, event_matrix)
    return(range_matrix)
  }
  # Let's plot graph
  # configuring x and y axis
  # the range of the x-axis will increase according to the number of events.
  # if the number of events is a two-digit number, x axis tick value is 10.
  x_range <- length(cusum_matrix$number_events)
  x_ticks <- 10^(floor(log10(x_range)))
  # first graph
  graph_first <- as.data.frame(find_score(range_tibble,1))
  p <- ggplot2::ggplot(data=graph_first, ggplot2::aes(x=event, y=event_score))+
    ggplot2::scale_x_continuous(limits = c(0, x_range), breaks=seq(0, x_range, x_ticks))+
    ggplot2::scale_y_continuous(limits = c(-h0 * 1.25, h1 * 1.2 ), breaks=seq(-10,10,1)) +
    ggplot2::geom_line(color="black", size=0.5) +
    ggplot2::geom_point(shape=124, size=3) +
    ggplot2::labs(x="Number of Events",y="LC and CUSUM Score",tag="",title="",subtitle="",caption="") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    #draw upper and lower decision limits
    ggplot2::geom_abline(ggplot2::aes(intercept=-h0,slope=0), size=0.5, linetype = "dotted", color = 'black') +
    ggplot2::geom_abline(ggplot2::aes(intercept=0,slope=0), size=0.5, linetype = "solid", color = 'black') +
    ggplot2::geom_abline(ggplot2::aes(intercept=h1,slope=0), size=0.5, linetype = "dotted", color = 'black') +
    ggplot2::annotate("text", x = x_range * 0.95, y = -h0 - 0.25, label = "LDL") +
    ggplot2::annotate("text", x = x_range * 0.95, y = h1 + 0.25, label = "UDL")
  # add more graph
  if (length(range_tibble$range_end) >1){
    for (i in 2:length(range_tibble$range_end)){
      graph_next <- as.data.frame(find_score(range_tibble,i))
      p <- p + ggplot2::geom_line(data=graph_next, ggplot2::aes(x=event, y=event_score), color="black", size=0.5)+
        ggplot2::geom_point(data=graph_next, ggplot2::aes(x=event, y=event_score), shape=124, size=2)
    }
  }
  # if learning == true and learning ends during observation, mark the learning point as a vertical line with the point.
  if (( learning == TRUE) & (length(range_tibble$range_end) >1)){
    p <- p +
      ggplot2::geom_vline(ggplot2::aes(xintercept=range_end[1]),
                          color='darkgoldenrod1', size=0.5, linetype = "dashed")+
      ggplot2::geom_point(ggplot2::aes(range_end[1], cusum_matrix$cusum_score[range_end[1]]), shape=19, size=2, color="darkgoldenrod1")
  }
  # if reset happen during monitor, mark the reset point as the red point
  if (sum(cusum_matrix$reset_signal)>0){
    reset_events <- cusum_matrix$number_events[cusum_matrix$reset_signal == 1]
    reset_scores <- cusum_matrix$cusum_score[cusum_matrix$reset_signal == 1]
    reset_matrix <- as.data.frame(tibble::tibble(reset_events,reset_scores))
    p <- p +
      ggplot2::geom_point(data = reset_matrix, ggplot2::aes(reset_matrix[,1], reset_matrix[,2]), shape=19, size=2, color="red")
  }
  return(p)
}
# end of function
