# Tests for mapADMETColumns and computeADMETProperties (pure-R, no suggested deps).

test_that("mapADMETColumns renames and converts columns", {
  d <- data.frame(
    name = c("a", "b"),
    mol_weight = c("300", "650"),
    lipophilicity = c("2", "6"),
    smiles = c("CCO", "CCCCCCCC"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mapping <- c(
    name = "Name",
    mol_weight = "MW",
    lipophilicity = "LogP",
    smiles = "SMILES"
  )
  out <- mapADMETColumns(d, mapping, calculate_cdk = FALSE)
  expect_true("MW" %in% names(out))
  expect_true("LogP" %in% names(out))
  expect_true(is.numeric(out$MW))
  expect_true(is.numeric(out$LogP))
})

test_that("mapADMETColumns handles empty mapping", {
  d <- data.frame(a = 1, b = 2)
  out <- mapADMETColumns(d, c(a = "None", b = "None"), calculate_cdk = FALSE)
  expect_true(is.data.frame(out))
})

test_that("mapADMETColumns converts ADMET numeric to categorical", {
  d <- data.frame(
    SMILES = c("CCO", "CC(=O)O"),
    hia_prob = c(0.8, 0.3),
    bbb_prob = c(0.6, 0.2),
    pgp_prob = c(0.9, 0.1),
    check.names = FALSE
  )
  mapping <- c(
    SMILES = "SMILES",
    hia_prob = "GI Absorption_num",
    bbb_prob = "BBB Permeant_num",
    pgp_prob = "Pgp Substrate_num"
  )
  out <- mapADMETColumns(d, mapping, calculate_cdk = FALSE)
  expect_true("GI absorption" %in% names(out))
  expect_true("BBB permeant" %in% names(out))
  expect_true("Pgp substrate" %in% names(out))
  expect_equal(out[1, "GI absorption"], "High")
  expect_equal(out[2, "GI absorption"], "Low")
  expect_equal(out[1, "BBB permeant"], "Yes")
  expect_equal(out[2, "BBB permeant"], "No")
})

test_that("mapADMETColumns skips None-mapped columns", {
  d <- data.frame(
    SMILES = c("CCO"),
    extra_col = c("unrelated"),
    check.names = FALSE
  )
  mapping <- c(
    SMILES = "SMILES",
    extra_col = "None"
  )
  out <- mapADMETColumns(d, mapping, calculate_cdk = FALSE)
  expect_true("SMILES" %in% names(out))
  expect_true("extra_col" %in% names(out))
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
