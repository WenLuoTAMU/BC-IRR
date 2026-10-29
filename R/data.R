#' Task Engagement and Math Performance (Ota & DuPaul, 2002)
#'
#' Data from a study examining the effects of supplemental computer instruction
#' on math performance for three students with ADHD. This dataset follows an
#' A-B replication design across subjects.
#'
#' @format A data frame with 46 rows and 4 variables:
#' \describe{
#'   \item{Subject}{Identifier for the student (1, 2, or 3)}
#'   \item{Session}{Session number (time variable)}
#'   \item{Y}{Outcome variable (Math performance/score)}
#'   \item{Phase}{Experimental phase indicator (0 = Baseline, 1 = Intervention/Computer Instruction)}
#' }
#' @source Ota, K. R., & DuPaul, G. J. (2002). Task-engagement and mathematics
#' performance in children with attention-deficit hyperactivity disorder: Effects
#' of supplemental computer instruction. \emph{School Psychology Quarterly}, 17(3), 242–257.
"OtaDuPaul2002"

#' Task Engagement in Language Arts (Cihak et al., 2010)
#'
#' Data from a study examining the effects of self-modeling static-picture
#' prompts delivered via a handheld computer on task engagement for three
#' students with autism. This dataset includes each student's initial
#' baseline and intervention phases in the language arts setting.
#'
#' @format A data frame with 46 rows and 5 variables:
#' \describe{
#'   \item{Subject}{Student identifier:
#'     1 = Adam, 2 = Jordan, 3 = Richard.}
#'   \item{Session}{Original session number within student.}
#'   \item{Y}{Task engagement, expressed as the percentage of observation
#'     intervals with task engagement (0–100).}
#'   \item{Phase}{Experimental phase indicator
#'     (0 = Baseline, 1 = Handheld Computer Intervention).}
#'   \item{Total}{Total number of observation intervals per session (60).}
#' }
#'
#' @details
#' Outcome values were extracted from Figures 1–3 and are approximate.
#' The dataset includes Adam's sessions 1–10, Jordan's sessions 1–23,
#' and Richard's sessions 1–13. Withdrawal and subsequent intervention
#' phases are excluded. Eight outcomes are missing for Jordan,
#' representing sessions without recorded observations.
#'
#' @source Cihak, D. F., Wright, R., & Ayres, K. M. (2010).
#' Use of self-modeling static-picture prompts via a handheld computer
#' to facilitate self-monitoring in the general education classroom.
#' \emph{Education and Training in Autism and Developmental Disabilities},
#' 45(1), 136–149.
"Cihak2010"
