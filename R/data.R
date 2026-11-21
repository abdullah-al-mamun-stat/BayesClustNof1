#' Simulated series of N-of-1 trials with two responders
#'
#' A series of 10 N-of-1 trials simulated with the data-generating model of
#' Simulation study 2 in the accompanying paper: 5 blocks, each with one
#' treatment period and one control period in random order, and 6 daily
#' outcomes per period (600 observations). Patients 3 and 7 have a treatment
#' effect of 1.5 and all other patients have no effect. Errors follow an AR(1)
#' process with \eqn{\rho = 0.4} that restarts in every period, and the time,
#' block and period effects are non-zero.
#'
#' @format A data frame with 600 rows, sorted by patient and time, and columns
#' \describe{
#'   \item{ID}{patient (1 to 10)}
#'   \item{Y}{outcome}
#'   \item{A}{treatment indicator (1 = active, 0 = control)}
#'   \item{t}{day within the treatment period (0 to 5)}
#'   \item{B}{block (1 to 5)}
#'   \item{P}{position of the period within its block (1 or 2)}
#'   \item{beta_true}{the true treatment effect of the patient}
#' }
#' @examples
#' head(nof1_example)
#' tapply(nof1_example$beta_true, nof1_example$ID, unique)
"nof1_example"

#' Raynaud phenomenon series of N-of-1 trials: placebo vs sildenafil 40 mg
#'
#' Daily Raynaud Condition Scores from the series of N-of-1 trials of
#' on-demand sildenafil by Roustit et al. (2018), restricted to the placebo
#' and sildenafil 40 mg weeks, as analysed in the accompanying paper. Each of
#' the 38 patients completed two to five blocks; each block contained one week
#' of placebo, one of 40 mg and one of 80 mg in random order, and only the
#' placebo and 40 mg weeks are kept here.
#'
#' @format A data frame with 1,537 rows, sorted by patient, block, week and day:
#' \describe{
#'   \item{ID}{patient (\code{id_total} in the original data)}
#'   \item{Y}{Raynaud Condition Score, 0 to 10; lower values mean fewer symptoms (\code{rcs})}
#'   \item{A}{treatment: 1 = sildenafil 40 mg, 0 = placebo}
#'   \item{t}{day of the week, 1 to 7 (\code{jour_semaine})}
#'   \item{B}{block (\code{cycle_patient})}
#'   \item{P}{position of the week among the two retained weeks of its block (1 = earlier)}
#'   \item{week}{week number within the patient's trial (\code{semaine_patient})}
#' }
#' @source Roustit M, Giai J, Gaget O, et al. (2018). Data from: On-demand
#'   sildenafil as a treatment for Raynaud phenomenon: a series of n-of-1 trials.
#'   Dryad Digital Repository, \doi{10.5061/dryad.c670tq2} (CC0 1.0).
#' @references Roustit M, Giai J, Gaget O, et al. (2018). On-demand sildenafil
#'   as a treatment for Raynaud phenomenon: a series of n-of-1 trials.
#'   \emph{Annals of Internal Medicine}, 169, 694--703.
#' @examples
#' head(raynaud)
#' length(unique(raynaud$ID))
"raynaud"
