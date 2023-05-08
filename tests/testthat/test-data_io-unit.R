# Unit test load_demohaz_data in data_io.R

library(testthat)

test_that("Testing load_demohaz_data() function", {
  # (1) Ache Lifetable
  expect_error(
    ache_lifetable <- load_demohaz_data("ache"),
    NA
  )
  
  expect_equal(
    dim(ache_lifetable),
    c(78,9)
  )
  
  # (2) Hadza Lifetable
  expect_error(
    hadza_lifetable <- load_demohaz_data("hadza"),
    NA
  )
  
  expect_equal(
    dim(hadza_lifetable),
    c(18,9)
  )
  
  # (3) Hiwi Lifetable
  expect_error(
    hiwi_lifetable <- load_demohaz_data("hiwi"),
    NA
  )
  
  expect_equal(
    dim(hiwi_lifetable),
    c(16,9)
  )
  
  # (4) New Tsimane Lifetable
  expect_error(
    tsimane_lifetable <- load_demohaz_data("tsimane"),
    NA
  )
  
  expect_equal(
    dim(tsimane_lifetable),
    c(17,9)
  )

  # Test for an invalid dataset name
  expect_error(load_demohaz_data("invalid_dataset_name"))
})

