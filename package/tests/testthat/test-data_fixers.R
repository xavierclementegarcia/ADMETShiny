# Tests for fixSwissADME and computeADMETProperties (pure-R, no suggested deps).

test_that("fixSwissADME coerces numeric and character columns", {
  d <- data.frame(
    Name = c("a", "b"),
    MW = c("300", "650"),
    "Consensus Log P" = c("2", "6"),
    "#H-bond acceptors" = c("4", "12"),
    "#H-bond donors" = c("2", "7"),
    "GI absorption" = c("High", "Low"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- fixSwissADME(d)
  expect_true(is.numeric(out$MW))
  expect_true(is.numeric(out[["Consensus Log P"]]))
  expect_true(is.character(out[["GI absorption"]]))
  # LogP column should be created from Consensus Log P
  expect_true("LogP" %in% names(out))
  expect_equal(out$LogP, out[["Consensus Log P"]])
})

test_that("fixSwissADME warns on missing required columns", {
  d <- data.frame(Name = "a")
  expect_warning(fixSwissADME(d), "Missing expected columns")
})

test_that("fixSwissADME converts empty strings to NA", {
  d <- data.frame(
    MW = c("300", ""),
    "Consensus Log P" = c("2", "6"),
    "#H-bond acceptors" = c("4", "12"),
    "#H-bond donors" = c("2", "7"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- fixSwissADME(d)
  expect_true(is.na(out$MW[2]))
})

test_that("computeADMETProperties adds GI, BBB and Pgp columns", {
  d <- data.frame(
    LogP = c(2.5, 6),
    TPSA = c(70, 220),
    MW = c(300, 650),
    "#H-bond donors" = c(2, 7),
    "#H-bond acceptors" = c(4, 12),
    check.names = FALSE
  )
  out <- computeADMETProperties(d)
  expect_true("GI absorption" %in% names(out))
  expect_true("BBB permeant" %in% names(out))
  expect_true("Pgp substrate" %in% names(out))
  # compound inside the HIA ellipse -> High
  expect_equal(out[1, "GI absorption"], "High")
  # compound outside the HIA ellipse (TPSA = 220 is beyond the ellipse) -> Low
  expect_equal(out[2, "GI absorption"], "Low")
})
