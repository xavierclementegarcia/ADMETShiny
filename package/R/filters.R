# ---------------------------------------------------------------------------
# filters.R
#
# Drug-likeness filters based on the classical literature rules.
#
# Each *Filter() combines two layers of evidence:
#   (a) the individual physicochemical property limits (adjustable from the UI),
#       taken from the original publication of each rule;
#   (b) the "#violations" column pre-computed by computeViolationColumns(),
#       as an additional safety layer.
#
# The filters use the generic `LogP` column (not SwissADME-specific
# MLOGP/WLOGP/XLOGP3). The thresholds are from the original publications:
#   Lipinski (1997): LogP <= 5 (original Rule-of-Five)
#   Ghose (1999):    LogP -0.4 to 5.6
#   Veber (2002):    RB <= 10 AND (TPSA <= 140 OR HBA+HBD <= 12)
#   Egan (2000):     LogP <= 5.88, TPSA <= 131.6
#   Muegge (2001):   LogP -2 to 5
#
# References:
#   Lipinski et al., 1997 (Rule of Five)
#   Ghose et al., 1999    (Ghose filter / qualifying range)
#   Veber et al., 2002    (oral bioavailability rule)
#   Egan et al., 2000     (Egan / PSA-LogP rule)
#   Muegge et al., 2001   (pharmacophore point filter)
# ---------------------------------------------------------------------------

## ---------------------------- helpers --------------------------------------

#' Check that required columns exist before filtering
#'
#' Internal helper used by every drug-likeness filter to fail early with an
#' informative message when an expected column is missing from the input data.
#'
#' @param data A data.frame to check.
#' @param required_columns Character vector of column names that must be
#'   present in \code{data}.
#' @param filter_name Character scalar; name of the calling filter, used in the
#'   error message.
#' @return Invisible \code{NULL}; called for its side-effect of stopping when
#'   columns are missing.
#' @keywords internal
.checkColumns <- function(data, required_columns, filter_name) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "[%s] Missing required columns: %s",
        filter_name, paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

## ---------------------------- Lipinski ------------------------------------
## Rule of Five (Lipinski et al., 1997): MW <= 500, LogP <= 5, HBA <= 10,
## HBD <= 5. Up to 1 violation is usually considered acceptable, but it is
## left configurable.

#' Lipinski Rule-of-Five filter
#'
#' Filters compounds according to the Rule of Five (Lipinski et al., 1997):
#' molecular weight, LogP, number of H-bond acceptors and donors, plus the
#' pre-computed number of Lipinski violations.
#'
#' Uses the generic \code{LogP} column with the original threshold of 5 (not
#' the SwissADME-specific MLOGP <= 4.15).
#'
#' @param data A data.frame with the standard column schema (processed by
#'   \code{\link{fixSwissADME}}, \code{\link{fixADMETlab}},
#'   \code{\link{fixDeepPK}}, or \code{\link{mapCDKDescriptors}}).
#' @param mw Numeric. Maximum molecular weight. Default 500.
#' @param logp Numeric. Maximum LogP. Default 5.
#' @param hba Numeric. Maximum number of H-bond acceptors. Default 10.
#' @param hbd Numeric. Maximum number of H-bond donors. Default 5.
#' @param violations Integer. Maximum tolerated Lipinski violations. Default 0.
#' @return A data.frame with the rows of \code{data} that satisfy all
#'   thresholds.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = 300, LogP = 2, "#H-bond acceptors" = 4,
#'   "#H-bond donors" = 2, "Lipinski #violations" = 0, check.names = FALSE)
#' lipinskiFilter(d)
#' }
#' @export
#' @references Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J.
#'   (1997). Experimental and computational approaches to estimate solubility
#'   and permeability in drug discovery and development settings.
#'   \emph{Advanced Drug Delivery Reviews}, 23(1-3), 3-25.
lipinskiFilter <- function(
    data,
    mw = 500,
    logp = 5,
    hba = 10,
    hbd = 5,
    violations = 0) {

  required_columns <- c(
    "MW", "LogP", "#H-bond acceptors", "#H-bond donors",
    "Lipinski #violations"
  )
  .checkColumns(data, required_columns, "Lipinski")

  filter(
    data,
    MW <= mw,
    LogP <= logp,
    `#H-bond acceptors` <= hba,
    `#H-bond donors` <= hbd,
    `Lipinski #violations` <= violations
  )
}

## ----------------------------- Ghose ---------------------------------------
## Ghose et al., 1999 -- "qualifying range" (two-sided range):
##   MW: 160-480 | MR: 40-130 | LogP: -0.4 to 5.6 | #heavy atoms: 20-70

#' Ghose drug-likeness filter
#'
#' Filters compounds according to the Ghose qualifying range (Ghose et al.,
#' 1999): molecular weight, molar refractivity, LogP and number of heavy
#' atoms.
#'
#' @param data A data.frame with the standard column schema.
#' @param g_mw_min Numeric. Minimum molecular weight. Default 160.
#' @param g_mw_max Numeric. Maximum molecular weight. Default 480.
#' @param g_mr_min Numeric. Minimum molar refractivity. Default 40.
#' @param g_mr_max Numeric. Maximum molar refractivity. Default 130.
#' @param g_logp_min Numeric. Minimum LogP. Default -0.4.
#' @param g_logp_max Numeric. Maximum LogP. Default 5.6.
#' @param g_ha_min Numeric. Minimum number of heavy atoms. Default 20.
#' @param g_ha_max Numeric. Maximum number of heavy atoms. Default 70.
#' @param violations Integer. Maximum tolerated Ghose violations. Default 0.
#' @return A data.frame with the rows of \code{data} that satisfy all
#'   thresholds.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = 300, MR = 90, LogP = 2, "#Heavy atoms" = 22,
#'   "Ghose #violations" = 0, check.names = FALSE)
#' ghoseFilter(d)
#' }
#' @export
#' @references Ghose, A. K., Viswanadhan, V. N., & Wendoloski, J. J. (1999).
#'   A knowledge-based approach in designing combinatorial or medicinal
#'   chemistry libraries for drug discovery. 1. A qualitative and
#'   quantitative characterization of known drug databases.
#'   \emph{Journal of Combinatorial Chemistry}, 1(1), 55-68.
ghoseFilter <- function(
    data,
    g_mw_min = 160,
    g_mw_max = 480,
    g_mr_min = 40,
    g_mr_max = 130,
    g_logp_min = -0.4,
    g_logp_max = 5.6,
    g_ha_min = 20,
    g_ha_max = 70,
    violations = 0) {

  required_columns <- c(
    "MW", "MR", "LogP", "#Heavy atoms", "Ghose #violations"
  )
  .checkColumns(data, required_columns, "Ghose")

  filter(
    data,
    MW >= g_mw_min, MW <= g_mw_max,
    MR >= g_mr_min, MR <= g_mr_max,
    LogP >= g_logp_min, LogP <= g_logp_max,
    `#Heavy atoms` >= g_ha_min, `#Heavy atoms` <= g_ha_max,
    `Ghose #violations` <= violations
  )
}

## ------------------------------ Egan ----------------------------------------
## Egan et al., 2000: TPSA <= 131.6 and LogP <= 5.88 (PSA-LogP / BBB model).

#' Egan drug-likeness filter
#'
#' Filters compounds according to the Egan rule (Egan et al., 2000) based on
#' TPSA and LogP.
#'
#' @param data A data.frame with the standard column schema.
#' @param e_tpsa Numeric. Maximum TPSA. Default 131.6.
#' @param e_logp Numeric. Maximum LogP. Default 5.88.
#' @param violations Integer. Maximum tolerated Egan violations. Default 0.
#' @return A data.frame with the rows of \code{data} that satisfy all
#'   thresholds.
#' @examples
#' \dontrun{
#' d <- data.frame(TPSA = 80, LogP = 2, "Egan #violations" = 0,
#'   check.names = FALSE)
#' eganFilter(d)
#' }
#' @export
#' @references Egan, W. J., Merz, K. M., & Baldwin, J. J. (2000). Prediction
#'   of drug absorption using multivariate statistics. \emph{Journal of
#'   Medicinal Chemistry}, 43(21), 3867-3877.
eganFilter <- function(
    data,
    e_tpsa = 131.6,
    e_logp = 5.88,
    violations = 0) {

  required_columns <- c("TPSA", "LogP", "Egan #violations")
  .checkColumns(data, required_columns, "Egan")

  filter(
    data,
    LogP <= e_logp,
    TPSA <= e_tpsa,
    `Egan #violations` <= violations
  )
}

## ----------------------------- Muegge ---------------------------------------
## Muegge et al., 2001: MW 200-600 | LogP -2 to 5 | TPSA <= 150 |
##   HBA <= 10 | HBD <= 5 | rotatable bonds <= 15.
## Note: the original Muegge rule also includes "pharmacophore points >= 4"
## (HBA + HBD + rings), but rings are not directly available from
## SwissADME/ADMETlab/Deep-PK/CDK data. The aromatic heavy atoms
## criterion was removed because it is NOT part of the original Muegge rule
## and the threshold of <= 7 is far too restrictive (a single benzene ring
## has 6 aromatic heavy atoms; two rings have 10-12).

#' Muegge drug-likeness filter
#'
#' Filters compounds according to the Muegge pharmacophore-point filter
#' (Muegge et al., 2001): MW, LogP, HBA, HBD, TPSA and rotatable bonds.
#'
#' @param data A data.frame with the standard column schema.
#' @param m_mw_min Numeric. Minimum molecular weight. Default 200.
#' @param m_mw_max Numeric. Maximum molecular weight. Default 600.
#' @param m_logp_min Numeric. Minimum LogP. Default -2.
#' @param m_logp_max Numeric. Maximum LogP. Default 5.
#' @param m_hba Numeric. Maximum H-bond acceptors. Default 10.
#' @param m_hbd Numeric. Maximum H-bond donors. Default 5.
#' @param m_rb Numeric. Maximum rotatable bonds. Default 15.
#' @param m_tpsa Numeric. Maximum TPSA. Default 150.
#' @param violations Integer. Maximum tolerated Muegge violations. Default 0.
#' @return A data.frame with the rows of \code{data} that satisfy all
#'   thresholds.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = 300, LogP = 2, TPSA = 80,
#'   "#H-bond acceptors" = 4, "#H-bond donors" = 2,
#'   "#Rotatable bonds" = 3, "Muegge #violations" = 0, check.names = FALSE)
#' mueggeFilter(d)
#' }
#' @export
#' @references Muegge, I., Heald, S. L., & Brittelli, D. (2001). Simple
#'   selection criteria for drug-like chemical matter. \emph{Journal of
#'   Medicinal Chemistry}, 44(12), 1841-1846.
mueggeFilter <- function(
    data,
    m_mw_min = 200,
    m_mw_max = 600,
    m_logp_min = -2,
    m_logp_max = 5,
    m_hba = 10,
    m_hbd = 5,
    m_rb = 15,
    m_tpsa = 150,
    violations = 0) {

  required_columns <- c(
    "#Rotatable bonds", "MW", "LogP", "TPSA",
    "#H-bond acceptors", "#H-bond donors",
    "Muegge #violations"
  )
  .checkColumns(data, required_columns, "Muegge")

  filter(
    data,
    `#Rotatable bonds` <= m_rb,
    MW >= m_mw_min, MW <= m_mw_max,
    LogP >= m_logp_min, LogP <= m_logp_max,
    TPSA <= m_tpsa,
    `#H-bond acceptors` <= m_hba,
    `#H-bond donors` <= m_hbd,
    `Muegge #violations` <= violations
  )
}

## ----------------------------- Veber ---------------------------------------
## Veber et al., 2002: rotatable bonds <= 10 AND (TPSA <= 140 OR HBA+HBD <= 12).
## The two polarity conditions are alternative equivalents in the original
## paper, not cumulative; an OR is used to respect the rule.

#' Veber drug-likeness filter
#'
#' Filters compounds according to the Veber oral-bioavailability rule (Veber
#' et al., 2002): rotatable bonds and a polarity condition (TPSA or
#' HBA+HBD).
#'
#' @param data A data.frame with the standard column schema.
#' @param v_rb Numeric. Maximum rotatable bonds. Default 10.
#' @param v_tpsa Numeric. Maximum TPSA. Default 140.
#' @param v_hb_sum Numeric. Maximum H-bond acceptors + donors. Default 12.
#' @param violations Integer. Maximum tolerated Veber violations. Default 0.
#' @return A data.frame with the rows of \code{data} that satisfy all
#'   thresholds.
#' @examples
#' \dontrun{
#' d <- data.frame(TPSA = 80, "#Rotatable bonds" = 3,
#'   "#H-bond acceptors" = 4, "#H-bond donors" = 2,
#'   "Veber #violations" = 0, check.names = FALSE)
#' veberFilter(d)
#' }
#' @export
#' @references Veber, D. F., Johnson, S. R., Cheng, H. Y., Smith, B. R.,
#'   Ward, K. W., & Kopple, K. D. (2002). Molecular properties that influence
#'   the oral bioavailability of drug candidates. \emph{Journal of Medicinal
#'   Chemistry}, 45(12), 2615-2623.
veberFilter <- function(
    data,
    v_rb = 10,
    v_tpsa = 140,
    v_hb_sum = 12,
    violations = 0) {

  required_columns <- c(
    "#Rotatable bonds", "TPSA",
    "#H-bond acceptors", "#H-bond donors", "Veber #violations"
  )
  .checkColumns(data, required_columns, "Veber")

  data %>%
    mutate(.hb_sum = `#H-bond acceptors` + `#H-bond donors`) %>%
    filter(
      `#Rotatable bonds` <= v_rb,
      (TPSA <= v_tpsa | .hb_sum <= v_hb_sum),
      `Veber #violations` <= violations
    ) %>%
    select(-.hb_sum)
}

## --------------------------- Dispatcher --------------------------------------

#' Apply drug-likeness filters in sequence
#'
#' Applies the selected drug-likeness filters to a normalized data.frame. The
#' filters are applied in the order Lipinski, Veber, Ghose, Egan, Muegge; only
#' those named in \code{filters} are executed.
#'
#' @param data A data.frame already normalized (e.g. by \code{fixSwissADME},
#'   \code{fixADMETlab}, \code{fixDeepPK}, or \code{mapCDKDescriptors}).
#' @param filters Character vector with the names of the filters to apply. Any
#'   subset of \code{c("Lipinski", "Veber", "Ghose", "Egan", "Muegge")}.
#' @param lipinski,veber,ghose,egan,muegge Named lists of parameters forwarded
#'   to the corresponding filter function.
#' @return A data.frame with the rows of \code{data} that pass every selected
#'   filter.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = 300, LogP = 2, TPSA = 80, MR = 90,
#'   "#Heavy atoms" = 22, "#Rotatable bonds" = 3,
#'   "#H-bond acceptors" = 4, "#H-bond donors" = 2,
#'   "Lipinski #violations" = 0, "Ghose #violations" = 0,
#'   "Veber #violations" = 0, "Egan #violations" = 0,
#'   "Muegge #violations" = 0, check.names = FALSE)
#' applyFilters(d, filters = c("Lipinski", "Veber", "Egan"))
#' }
#' @export
#' @seealso \code{\link{lipinskiFilter}}, \code{\link{veberFilter}},
#'   \code{\link{ghoseFilter}}, \code{\link{eganFilter}},
#'   \code{\link{mueggeFilter}}.
applyFilters <- function(
    data,
    filters,
    lipinski = list(),
    veber = list(),
    ghose = list(),
    egan = list(),
    muegge = list()) {

  resultado <- data

  if ("Lipinski" %in% filters) {
    resultado <- do.call(lipinskiFilter, c(list(data = resultado), lipinski))
  }

  if ("Veber" %in% filters) {
    resultado <- do.call(veberFilter, c(list(data = resultado), veber))
  }

  if ("Ghose" %in% filters) {
    resultado <- do.call(ghoseFilter, c(list(data = resultado), ghose))
  }

  if ("Egan" %in% filters) {
    resultado <- do.call(eganFilter, c(list(data = resultado), egan))
  }

  if ("Muegge" %in% filters) {
    resultado <- do.call(mueggeFilter, c(list(data = resultado), muegge))
  }

  resultado
}
