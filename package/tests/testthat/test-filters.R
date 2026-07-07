# Tests for the drug-likeness filters and violation-column calculator.
# These exercise the pure-R logic and do not require any suggested package.

## ---- helper: a small synthetic dataset in the standard schema -----------
.make_test_data <- function() {
  data.frame(
    Name = c("druglike", "violator"),
    SMILES = c("CCO", "CCCCCCCCCCCCCCCCCCCCCCCCCC"),
    MW = c(300, 650),
    LogP = c(2, 7),
    TPSA = c(40, 160),
    MR = c(70, 150),
    "#H-bond acceptors" = c(4, 12),
    "#H-bond donors" = c(2, 7),
    "#Rotatable bonds" = c(3, 14),
    "#Heavy atoms" = c(20, 80),
    "#Aromatic heavy atoms" = c(6, 9),
    "Lipinski #violations" = c(0, 4),
    "Ghose #violations" = c(0, 4),
    "Veber #violations" = c(0, 1),
    "Egan #violations" = c(0, 2),
    "Muegge #violations" = c(0, 5),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

## ---- computeViolationColumns ---------------------------------------------
test_that("computeViolationColumns computes the five violation columns", {
  d <- .make_test_data()
  d <- computeViolationColumns(d)

  expect_true("Lipinski #violations" %in% names(d))
  expect_true("Ghose #violations" %in% names(d))
  expect_true("Veber #violations" %in% names(d))
  expect_true("Egan #violations" %in% names(d))
  expect_true("Muegge #violations" %in% names(d))

  # druglike compound has zero violations across the board
  expect_equal(as.numeric(d[1, "Lipinski #violations"]), 0)
  expect_true(as.numeric(d[2, "Lipinski #violations"]) > 0)
  expect_equal(as.numeric(d[1, "Egan #violations"]), 0)
  expect_true(as.numeric(d[2, "Egan #violations"]) > 0)
})

test_that("computeViolationColumns handles missing columns gracefully", {
  d <- data.frame(MW = c(300, 650))
  d <- computeViolationColumns(d)
  expect_true("Lipinski #violations" %in% names(d))
  expect_true(all(is.na(d[["Lipinski #violations"]])))
})

test_that("Veber violations include rotatable bonds (RB > 10)", {
  # A compound with RB = 15 but TPSA = 50 (low polarity) should still
  # have a Veber violation because of the rotatable bonds.
  d <- data.frame(
    MW = 300, LogP = 2, TPSA = 50,
    "#H-bond acceptors" = 4, "#H-bond donors" = 2,
    "#Rotatable bonds" = 15,
    check.names = FALSE
  )
  d <- computeViolationColumns(d)
  expect_equal(as.numeric(d[1, "Veber #violations"]), 1)

  # A compound with RB = 5 and TPSA = 50 should have 0 Veber violations
  d2 <- data.frame(
    MW = 300, LogP = 2, TPSA = 50,
    "#H-bond acceptors" = 4, "#H-bond donors" = 2,
    "#Rotatable bonds" = 5,
    check.names = FALSE
  )
  d2 <- computeViolationColumns(d2)
  expect_equal(as.numeric(d2[1, "Veber #violations"]), 0)
})

## ---- lipinskiFilter -------------------------------------------------------
test_that("lipinskiFilter keeps drug-like compounds", {
  d <- .make_test_data()
  out <- lipinskiFilter(d)
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
})

test_that("lipinskiFilter errors on missing columns", {
  d <- data.frame(MW = 300)
  expect_error(lipinskiFilter(d), "Lipinski")
})

## ---- ghoseFilter ----------------------------------------------------------
test_that("ghoseFilter respects the qualifying range", {
  d <- .make_test_data()
  out <- ghoseFilter(d)
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
})

## ---- eganFilter -----------------------------------------------------------
test_that("eganFilter respects TPSA and LogP limits", {
  d <- .make_test_data()
  out <- eganFilter(d)
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
})

## ---- mueggeFilter ---------------------------------------------------------
test_that("mueggeFilter respects the pharmacophore-point filter", {
  d <- .make_test_data()
  out <- mueggeFilter(d)
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
})

## ---- veberFilter ----------------------------------------------------------
test_that("veberFilter respects rotatable bonds and polarity", {
  d <- .make_test_data()
  out <- veberFilter(d)
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
  # the temporary .hb_sum column must be removed
  expect_false(".hb_sum" %in% names(out))
})

## ---- applyFilters ---------------------------------------------------------
test_that("applyFilters chains the selected filters", {
  d <- .make_test_data()
  out <- applyFilters(d, filters = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"))
  expect_true("druglike" %in% out$Name)
  expect_false("violator" %in% out$Name)
})

test_that("applyFilters with empty filters returns the input unchanged", {
  d <- .make_test_data()
  out <- applyFilters(d, filters = character(0))
  expect_equal(nrow(out), nrow(d))
})
