#' Two-by-two genetic correlation (2x2r) analysis using HDL
#'
#' Estimates six genetic correlations from four GWAS summary statistics
#' (two traits X, Y \eqn{\times} two samples) using HDL, tests for
#' cross-sample asymmetry, and infers the likely causal direction between
#' the two traits.
#'
#' @param gwasX1.df A data frame of GWAS summary statistics for trait X in
#'   sample 1. Must contain columns \code{SNP}, \code{A1}, \code{A2},
#'   \code{N}, and either \code{Z} or both \code{b} and \code{se}.
#' @param gwasY1.df A data frame of GWAS summary statistics for trait Y in
#'   sample 1. Same format as \code{gwasX1.df}.
#' @param gwasX2.df A data frame of GWAS summary statistics for trait X in
#'   sample 2. Same format as \code{gwasX1.df}.
#' @param gwasY2.df A data frame of GWAS summary statistics for trait Y in
#'   sample 2. Same format as \code{gwasX1.df}.
#' @param LD.path Path to the directory containing the LD reference panel
#'   (eigenvalues and eigenvectors). See \code{\link{HDL.rg}} for details.
#' @param Nref Sample size of the LD reference panel. Default is 335,265
#'   (UK Biobank).
#' @param N0.11 Number of overlapping individuals between the \code{gwasX1}
#'   and \code{gwasY1} cohorts (sample overlap within sample 1). Set to
#'   \code{0} if the two traits were measured in non-overlapping individuals.
#'   Default: \code{NULL}, which internally uses
#'   \code{min(gwasX1.df$N, gwasY1.df$N)} (the HDL default).
#' @param N0.22 Number of overlapping individuals between the \code{gwasX2}
#'   and \code{gwasY2} cohorts (sample overlap within sample 2). Same
#'   semantics as \code{N0.11}. Default: \code{NULL}.
#' @param output.file Path to a file where the summary results should be
#'   written. If \code{""} (default), results are printed to the console only.
#' @param eigen.cut Eigenvalue threshold passed to \code{\link{HDL.rg}} for
#'   each of the six calls. Either \code{"automatic"} (default) or a numeric
#'   value between 0 and 1.
#' @param fill.missing.N How to handle SNPs with missing sample size. One of
#'   \code{NULL} (default, remove such SNPs), \code{"min"}, \code{"max"}, or
#'   \code{"median"}. Passed to \code{\link{HDL.rg}}.
#' @param lim Optimisation tolerance passed to \code{\link{HDL.rg}}.
#'   Default: \code{exp(-18)}.
#' @param verbose Logical, \code{FALSE} by default. Whether to print the
#'   genetic covariance optimisation process on the console. Passed to
#'   \code{\link{HDL.rg}}.
#'
#' @details
#' The 2\eqn{\times}2r design estimates six pairwise genetic correlations
#' from the four input GWAS by calling \code{\link{HDL.rg}} on every
#' relevant pair:
#'
#' \tabular{lll}{
#'   \strong{Symbol} \tab \strong{Pair} \tab \strong{N0 assumption} \cr
#'   \eqn{r_G^{11}} \tab X(s1) vs Y(s1) \tab user-specified (\code{N0.11}) \cr
#'   \eqn{r_G^{22}} \tab X(s2) vs Y(s2) \tab user-specified (\code{N0.22}) \cr
#'   \eqn{r_G^{12}} \tab X(s1) vs Y(s2) \tab 0 (cross-sample) \cr
#'   \eqn{r_G^{21}} \tab X(s2) vs Y(s1) \tab 0 (cross-sample) \cr
#'   \eqn{r_G^X}    \tab X(s1) vs X(s2) \tab 0 (cross-sample) \cr
#'   \eqn{r_G^Y}    \tab Y(s1) vs Y(s2) \tab 0 (cross-sample) \cr
#' }
#'
#' \eqn{r_G^X} and \eqn{r_G^Y} quantify the genetic similarity of each
#' trait across the two samples and can be used to assess cohort
#' heterogeneity before interpreting asymmetry.
#'
#' \strong{Asymmetry test.}
#' The cross-sample directional asymmetry is
#' \deqn{\Delta = r_G^{12} - r_G^{21}.}
#' Under \eqn{H_0: \Delta = 0} the test statistic
#' \deqn{z_\Delta = \Delta \,/\, \sqrt{\mathrm{se}_{12}^2 + \mathrm{se}_{21}^2}}
#' follows a standard normal distribution (two-sided test). Independence of
#' the two estimates is assumed, which is valid when sample 1 and sample 2
#' are non-overlapping.
#'
#' \strong{Causal direction.}
#' When the asymmetry is significant, the causal direction is inferred by
#' comparing the sign of \eqn{\Delta} with the sign of the within-sample
#' difference \eqn{\delta_\mathrm{within} = r_G^{22} - r_G^{11}}:
#' \itemize{
#'   \item Same sign \eqn{\Rightarrow} X \eqn{\rightarrow} Y.
#'   \item Opposite signs \eqn{\Rightarrow} Y \eqn{\rightarrow} X.
#' }
#'
#' @return An invisible list with the following elements:
#' \describe{
#'   \item{\code{estimates.df}}{A data frame with columns \code{rG},
#'     \code{Estimate}, \code{SE}, and \code{P} for all six genetic
#'     correlations.}
#'   \item{\code{asymmetry}}{A named list: \code{Delta} (point estimate),
#'     \code{Delta.se} (SE), \code{z} (test statistic), \code{P}
#'     (two-sided p-value).}
#'   \item{\code{causal.direction}}{A named list: \code{delta.within}
#'     (\eqn{r_G^{22} - r_G^{11}}), \code{direction} (one of
#'     \code{"X->Y"}, \code{"Y->X"}, or
#'     \code{"Asymmetry not significant"}), \code{note}.}
#'   \item{\code{rg.list}}{Named list of the six raw \code{HDL.rg} result
#'     objects (\code{rG.11}, \code{rG.22}, \code{rG.12}, \code{rG.21},
#'     \code{rG.X}, \code{rG.Y}).}
#' }
#'
#' @author Yi Zheng, Xia Shen
#'
#' @references
#' Ning Z, Pawitan Y, Shen X (2020). High-definition likelihood inference of
#' genetic correlations across human complex traits. \emph{Nature Genetics}.
#'
#' @seealso \code{\link{HDL.rg}}
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' res <- HDL.2x2r(
#'   gwasX1.df = gwasX1, gwasY1.df = gwasY1,
#'   gwasX2.df = gwasX2, gwasY2.df = gwasY2,
#'   LD.path   = "/path/to/LD/reference",
#'   N0.11 = 0, N0.22 = 0   # traits measured independently within each sample
#' )
#' res$estimates.df
#' res$asymmetry
#' res$causal.direction
#' }
#' @export

HDL.2x2r <- function(
    gwasX1.df, gwasY1.df, gwasX2.df, gwasY2.df,
    LD.path,
    Nref          = 335265,
    N0.11         = NULL,
    N0.22         = NULL,
    output.file   = "",
    eigen.cut     = "automatic",
    fill.missing.N = NULL,
    lim           = exp(-18),
    verbose       = FALSE
) {

  ## ---- initialise output file -------------------------------------------- ##
  time.start <- date()
  .log <- function(..., append = TRUE) {
    msg <- paste0(...)
    cat(msg)
    if (output.file != "") cat(msg, file = output.file, append = append)
  }

  .log("HDL.2x2r analysis starts on ", time.start, "\n", append = FALSE)

  ## ---- resolve N0 defaults ------------------------------------------------ ##
  if (is.null(N0.11)) N0.11 <- min(c(gwasX1.df$N, gwasY1.df$N), na.rm = TRUE)
  if (is.null(N0.22)) N0.22 <- min(c(gwasX2.df$N, gwasY2.df$N), na.rm = TRUE)

  ## ---- helper: run one HDL.rg call ---------------------------------------- ##
  .run_hdl <- function(g1, g2, N0, label) {
    .log("\n--- Estimating ", label, " ---\n")
    tryCatch(
      HDL.rg(
        gwas1.df       = g1,
        gwas2.df       = g2,
        LD.path        = LD.path,
        Nref           = Nref,
        N0             = N0,
        output.file    = "",        # managed here, not inside HDL.rg
        eigen.cut      = eigen.cut,
        fill.missing.N = fill.missing.N,
        lim            = lim,
        verbose        = verbose
      ),
      error = function(e) {
        .log("WARNING: HDL.rg failed for ", label, ": ", conditionMessage(e), "\n")
        NULL
      }
    )
  }

  ## ---- run all six pairs -------------------------------------------------- ##
  res.11 <- .run_hdl(gwasX1.df, gwasY1.df, N0 = N0.11,
                     label = "r_G^{11} [X(s1) vs Y(s1)]")
  res.22 <- .run_hdl(gwasX2.df, gwasY2.df, N0 = N0.22,
                     label = "r_G^{22} [X(s2) vs Y(s2)]")
  res.12 <- .run_hdl(gwasX1.df, gwasY2.df, N0 = 0,
                     label = "r_G^{12} [X(s1) vs Y(s2)]")
  res.21 <- .run_hdl(gwasX2.df, gwasY1.df, N0 = 0,
                     label = "r_G^{21} [X(s2) vs Y(s1)]")
  res.X  <- .run_hdl(gwasX1.df, gwasX2.df, N0 = 0,
                     label = "r_G^X  [X(s1) vs X(s2)]")
  res.Y  <- .run_hdl(gwasY1.df, gwasY2.df, N0 = 0,
                     label = "r_G^Y  [Y(s1) vs Y(s2)]")

  ## ---- compile estimates table -------------------------------------------- ##
  .extract <- function(res, name) {
    if (is.null(res))
      return(data.frame(rG = name, Estimate = NA_real_, SE = NA_real_,
                        P = NA_real_, stringsAsFactors = FALSE))
    data.frame(rG = name, Estimate = res$rg, SE = res$rg.se, P = res$P,
               stringsAsFactors = FALSE)
  }

  estimates.df <- rbind(
    .extract(res.11, "r_G^{11}"),
    .extract(res.22, "r_G^{22}"),
    .extract(res.12, "r_G^{12}"),
    .extract(res.21, "r_G^{21}"),
    .extract(res.X,  "r_G^X"),
    .extract(res.Y,  "r_G^Y")
  )
  rownames(estimates.df) <- NULL

  ## ---- asymmetry test ----------------------------------------------------- ##
  rg12   <- estimates.df$Estimate[estimates.df$rG == "r_G^{12}"]
  rg21   <- estimates.df$Estimate[estimates.df$rG == "r_G^{21}"]
  se12   <- estimates.df$SE[estimates.df$rG == "r_G^{12}"]
  se21   <- estimates.df$SE[estimates.df$rG == "r_G^{21}"]

  Delta    <- rg12 - rg21
  Delta.se <- sqrt(se12^2 + se21^2)
  z.asym   <- Delta / Delta.se
  P.asym   <- 2 * pnorm(-abs(z.asym))

  asymmetry <- list(
    Delta    = Delta,
    Delta.se = Delta.se,
    z        = z.asym,
    P        = P.asym
  )

  ## ---- causal direction --------------------------------------------------- ##
  rg11         <- estimates.df$Estimate[estimates.df$rG == "r_G^{11}"]
  rg22         <- estimates.df$Estimate[estimates.df$rG == "r_G^{22}"]
  delta.within <- rg22 - rg11

  if (!is.na(P.asym) && P.asym < 0.05) {
    if (!is.na(delta.within) && sign(Delta) == sign(delta.within)) {
      direction <- "X->Y"
      note <- paste0(
        "Delta (r_G^{12} - r_G^{21}) = ", round(Delta, 4),
        " and delta.within (r_G^{22} - r_G^{11}) = ", round(delta.within, 4),
        " share the same sign, supporting X -> Y."
      )
    } else {
      direction <- "Y->X"
      note <- paste0(
        "Delta (r_G^{12} - r_G^{21}) = ", round(Delta, 4),
        " and delta.within (r_G^{22} - r_G^{11}) = ", round(delta.within, 4),
        " have opposite signs, supporting Y -> X."
      )
    }
  } else {
    direction <- "Asymmetry not significant"
    note <- paste0(
      "Asymmetry test P = ", formatC(P.asym, format = "e", digits = 2),
      " (>= 0.05). Causal direction cannot be determined from asymmetry alone."
    )
  }

  causal.direction <- list(
    delta.within = delta.within,
    direction    = direction,
    note         = note
  )

  ## ---- print summary ------------------------------------------------------ ##
  .fmt <- function(x) {
    if (is.na(x)) return("NA")
    if (abs(x) < 1e-4) formatC(x, format = "e", digits = 2) else round(x, 4)
  }

  .log("\n\n========== HDL.2x2r Results ==========\n")
  .log("\nSix genetic correlations:\n")
  .log(sprintf("  %-12s  %8s  %8s  %12s\n", "rG", "Estimate", "SE", "P"))
  for (i in seq_len(nrow(estimates.df))) {
    .log(sprintf("  %-12s  %8s  %8s  %12s\n",
                 estimates.df$rG[i],
                 .fmt(estimates.df$Estimate[i]),
                 .fmt(estimates.df$SE[i]),
                 .fmt(estimates.df$P[i])))
  }

  .log("\nAsymmetry test (H0: Delta = r_G^{12} - r_G^{21} = 0):\n")
  .log("  Delta    = ", .fmt(Delta),    "\n")
  .log("  Delta.se = ", .fmt(Delta.se), "\n")
  .log("  z        = ", .fmt(z.asym),   "\n")
  .log("  P        = ", formatC(P.asym, format = "e", digits = 2), "\n")

  .log("\nCausal direction: ", direction, "\n")
  .log("  ", note, "\n")
  .log("======================================\n\n")

  end.time <- date()
  .log("Analysis finished at ", end.time, "\n")

  ## ---- return ------------------------------------------------------------- ##
  invisible(list(
    estimates.df     = estimates.df,
    asymmetry        = asymmetry,
    causal.direction = causal.direction,
    rg.list = list(
      rG.11 = res.11,
      rG.22 = res.22,
      rG.12 = res.12,
      rG.21 = res.21,
      rG.X  = res.X,
      rG.Y  = res.Y
    )
  ))
}
